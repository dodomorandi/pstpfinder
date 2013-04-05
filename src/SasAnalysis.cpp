/*
 *  This file is part of PSTP-finder, an user friendly tool to analyze GROMACS
 *  molecular dynamics and find transient pockets on the surface of proteins.
 *  Copyright (C) 2011 Edoardo Morandi.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "SasAnalysis.h"

#include <sys/sysinfo.h>

#include <iterator>
#include <string>
#include <iostream>
#include <sstream>
#include <type_traits>

namespace PstpFinder
{
  template<typename T>
  class SasAnalysis<T,
      typename enable_if<
        not is_base_of<base_stream(basic_istream, T), T>::value and
        is_base_of<base_stream(basic_ostream, T), T>::value>::type> :
      public SasAnalysis_Write<T>
  {
    public:
      SasAnalysis(unsigned int nAtoms, const Gromacs& gromacs,
                  Session<T>& sessionFile) :
          SasAnalysis_Write<T>(nAtoms, gromacs, sessionFile) {}
      SasAnalysis(const Gromacs& gromacs, Session<T>& sessionFile) :
          SasAnalysis_Write<T>(gromacs, sessionFile) {}
      SasAnalysis(const Gromacs& gromacs, const string& sessionFileName) :
          SasAnalysis_Write<T>(gromacs, sessionFileName) {}
  };

  template<typename T>
  class SasAnalysis<T,
      typename enable_if<
        is_base_of<base_stream(basic_istream, T), T>::value and
        is_base_of<base_stream(basic_ostream, T), T>::value>::type> :
      public SasAnalysis_Write<T>
  {
    public:
      SasAnalysis(unsigned int nAtoms, const Gromacs& gromacs,
                  Session<T>& sessionFile) :
          SasAnalysis_Write<T>(nAtoms, gromacs, sessionFile) { init(); }
      SasAnalysis(const Gromacs& gromacs, Session<T>& sessionFile) :
          SasAnalysis_Write<T>(gromacs, sessionFile) { init(); }
      SasAnalysis(const Gromacs& gromacs, const string& sessionFileName) :
          SasAnalysis_Write<T>(gromacs, sessionFileName) { init(); }

    private:
      typedef SasAnalysis_Write<T> Base;

      void
      init()
      {
          size_t sasAtomSize(Base::serializer->getSerializedSize(SasAtom()));
          Base::sasMetaStream.seekg(0);
          Base::changeable = false;
          Base::serializer = new Serializer<MetaStream<T>>(Base::sasMetaStream);

          unsigned int validatedChunks(0);
          unsigned int chunkSize;
          unsigned long totalFrames(0);
          streampos backupPosition;

          *Base::serializer >> chunkSize;
          backupPosition = Base::sasMetaStream.tellg();

          while(not Base::sasMetaStream.eof())
          {
            Base::sasMetaStream.seekg(chunkSize * sasAtomSize * Base::nAtoms,
                                      ios_base::cur);

            if(Base::sasMetaStream.eof())
            {
              Base::sasMetaStream.seekg(backupPosition);
              break;
            }

            totalFrames += chunkSize;
            validatedChunks++;
            backupPosition = Base::sasMetaStream.tellg();
            *Base::serializer >> chunkSize;
          }

          Base::sasMetaStream.clear();
          if(totalFrames == 0)
            Base::sasMetaStream.seekg(0);
          else
            Base::sasMetaStream.seekg(0, ios_base::end);

          Base::readFrames = totalFrames;
          Base::sasMetaStream.seekp(Base::sasMetaStream.tellg());
          Base::analysisThread = new SasAnalysisThread<T>(*this);
      }
  };

  template<typename T>
  SasAnalysis_Base<T>::SasAnalysis_Base(unsigned int nAtoms,
                                        const Gromacs& gromacs,
                                        Session<T>& session) :
      sasMetaStream(session.getSasStream())
  {
    this->nAtoms = nAtoms;
    this->gromacs = &gromacs;
    init();
  }

  template<typename T>
  SasAnalysis_Base<T>::SasAnalysis_Base(const Gromacs& gromacs,
                                        Session<T>& session) :
      sasMetaStream(session.getSasStream())
  {
    nAtoms = gromacs.getGroup("Protein").size();
    this->gromacs = &gromacs;
    init();
  }

  template<typename T>
  SasAnalysis_Base<T>::SasAnalysis_Base(const Gromacs& gromacs,
                                        const string& sessionFileName) :
      rawSession(sessionFileName), sasMetaStream(rawSession.getSasStream())
  {
    nAtoms = gromacs.getGroup("Protein").size();
    this->gromacs = &gromacs;
    init();
  }

  template<typename T>
  void
  SasAnalysis_Base<T>::init()
  {
    changeable = true;

    struct sysinfo info;
    if(sysinfo(&info) == 0)
    {
      maxBytes = info.freeram * 0.8; // 80% of free ram
      if(maxBytes > 1073741824) // 1 GB
        maxBytes = 1073741824;
    }
    else
      maxBytes = 134217728;

    maxChunk = 8388608;
  }

  template<typename T>
  SasAnalysis_Read<T>::~SasAnalysis_Read()
  {
    if(Base::changeable)
      return;

    Base::analysisThread->stop();
    delete Base::analysisThread;

    delete Base::serializer;
  }

  template<typename T>
  SasAnalysis_Write<T>::~SasAnalysis_Write()
  {
    if(Base::changeable)
      return;

    if(Base::frames.size() != 0)
      flush();
    Base::analysisThread->waitForFlush();
    Base::analysisThread->stop();
    delete Base::analysisThread;

    delete Base::serializer;
    Base::sasMetaStream.close();
  }

  // TODO: Add a push_back(vector<SasAtom>&&)
  template<typename T>
  const SasAnalysis_Write<T>&
  SasAnalysis_Write<T>::push_back(const vector<SasAtom>& sasAtoms)
  {
    if(Base::changeable)
    {
      Base::changeable = false;
      Base::serializer = new Serializer<MetaStream<T>>(Base::sasMetaStream);

      Base::analysisThread = new SasAnalysisThreadType(*this);
    }

    Base::bufferMutex.lock();
    Base::frames.push_back(sasAtoms);

    bool cond = (Base::frames.size() == Base::maxFrames);
    Base::bufferMutex.unlock();

    if(cond)
      flush();
    return *this;
  }

  template<typename T>
  SasAnalysis_Read<T>&
  SasAnalysis_Read<T>::removeMe(const vector<SasAtom>*& sasAtom)
  {
    static std::vector<vector<SasAtom>>::const_iterator i;

    if(Base::changeable)
    {
      Base::changeable = false;
      Base::serializer = new Serializer<MetaStream<T>>(Base::sasMetaStream);

      Base::analysisThread = new SasAnalysisThreadType(*this);
    }

    if(Base::frames.size() == 0)
    {
      Base::bufferMutex.lock();
      while(Base::bufferCount == 0)
      {
        unique_lock<mutex> bufferCountLock(Base::bufferCountMutex);
        Base::bufferMutex.unlock();
        Base::bufferCountCondition.wait(bufferCountLock);
        Base::bufferMutex.lock();
      }
      Base::frames = Base::chunks.front();
      i = Base::frames.begin();

      Base::bufferCount--;
      Base::bufferMutex.unlock();
    }
    else if(i == Base::frames.end())
    {
      Base::bufferMutex.lock();

      Base::frames.clear();
      Base::chunks.pop_front();

      if(Base::bufferCount == 0)
      {
        Base::bufferMutex.unlock();
        unique_lock<mutex> bufferCountLock(Base::bufferCountMutex);
        Base::bufferCountCondition.wait(bufferCountLock);
        Base::bufferMutex.lock();
      }

      Base::bufferCount--;

      if(Base::chunks.size() == 0)
      {
        Base::bufferMutex.unlock();
        sasAtom = nullptr;
        return *this;
      }

      Base::frames = Base::chunks.front();
      i = Base::frames.begin();

      Base::analysisThread->wakeUp();
      Base::bufferMutex.unlock();
    }

    sasAtom = &(*i++);

    return *this;
  }

  template<typename T>
  typename SasAnalysis_Read<T>::const_iterator
  SasAnalysis_Read<T>::begin()
  {
    if(beginTriggered)
      throw;

    beginTriggered = true;

    if(Base::changeable)
    {
      Base::changeable = false;
      Base::serializer = new Serializer<MetaStream<T>>(Base::sasMetaStream);

      Base::analysisThread = new SasAnalysisThreadType(*this);
    }

    return const_iterator(this);
  }

  template<typename T>
  typename SasAnalysis_Read<T>::const_iterator
  SasAnalysis_Read<T>::end()
  {
    return const_iterator(this, std::end(Base::frames));
  }

  template<typename T>
  void
  SasAnalysis_Read<T>::const_iterator::handleEmptyFrames()
  {
    parent->bufferMutex.lock();
    while(parent->bufferCount == 0)
    {
      unique_lock<mutex> bufferCountLock(parent->bufferCountMutex);
      parent->bufferMutex.unlock();
      parent->bufferCountCondition.wait(bufferCountLock);
      parent->bufferMutex.lock();
    }
    parent->frames = parent->chunks.front();
    iter = std::begin(parent->frames);

    parent->bufferCount--;
    parent->bufferMutex.unlock();
  }

  template<typename T>
  typename SasAnalysis_Read<T>::const_iterator&
  SasAnalysis_Read<T>::const_iterator::operator++()
  {
    if(parent->frames.size() == 0)
      handleEmptyFrames();
    else if(++iter == std::end(parent->frames))
    {
      parent->bufferMutex.lock();

      parent->frames.clear();
      parent->chunks.pop_front();

      if(parent->bufferCount == 0)
      {
        parent->bufferMutex.unlock();
        unique_lock<mutex> bufferCountLock(parent->bufferCountMutex);
        parent->bufferCountCondition.wait(bufferCountLock);
        parent->bufferMutex.lock();
      }

      parent->bufferCount--;

      if(parent->chunks.size() == 0)
      {
        parent->bufferMutex.unlock();
        iter = std::end(parent->frames);
        return *this;
      }

      parent->frames = parent->chunks.front();
      iter = std::begin(parent->frames);

      parent->analysisThread->wakeUp();
      parent->bufferMutex.unlock();
    }

    return *this;
  }

  template<typename T>
  typename SasAnalysis_Read<T>::const_iterator::reference
  SasAnalysis_Read<T>::const_iterator::operator*() const
  {
    return iter.operator*();
  }

  template<typename T>
  typename SasAnalysis_Read<T>::const_iterator::pointer
  SasAnalysis_Read<T>::const_iterator::operator->() const
  {
    return iter.operator->();
  }

  template<typename T>
  bool
  SasAnalysis_Read<T>::const_iterator::operator!=(
      const SasAnalysis_Read<T>::const_iterator& other) const
  {
    return parent != other.parent or iter != other.iter;
  }

  template<typename T>
  void
  SasAnalysis_Write<T>::flush()
  {
    if(Base::changeable)
    {
      Base::changeable = false;
      Base::serializer = new Serializer<MetaStream<T>>(Base::sasMetaStream);

      Base::analysisThread = new SasAnalysisThreadType(*this);
    }

    Base::bufferMutex.lock();
    while(Base::bufferCount == 0)
    {
      unique_lock<mutex> bufferCountLock(Base::bufferCountMutex);
      Base::bufferMutex.unlock();
      Base::bufferCountCondition.wait(bufferCountLock);
      Base::bufferMutex.lock();
    }
    Base::bufferCount--;
    Base::chunks.push_back(move(Base::frames));

    Base::frames = vector<vector<SasAtom>>();
    Base::analysisThread->wakeUp();
    Base::bufferMutex.unlock();
  }

  template<typename T>
  bool
  SasAnalysis_Write<T>::save()
  {
    if(Base::chunks.empty())
      return false;

    dumpChunk(Base::chunks.front(), *Base::serializer);

    return true;
  }

  template<typename T>
  bool
  SasAnalysis_Read<T>::open()
  {
    Base::sasMetaStream.peek();
    if(Base::sasMetaStream.eof())
      return false;

    vector<vector<SasAtom>> chunk = loadChunk(*Base::serializer);
    unsigned long chunkSize = chunk.capacity()
                              * (sizeof(vector<SasAtom>)
                                 + sizeof(SasAtom) * Base::nAtoms)
                              + sizeof(vector<vector<SasAtom>> );
    if(chunkSize > Base::maxChunk)
    {
      // bufferMutex already lock_upgrade from ThreadOpen
      while(Base::bufferCount != 0)
      {
        unique_lock<mutex> lock(Base::bufferCountMutex);
        Base::bufferMutex.unlock();
        Base::bufferCountCondition.wait(lock);
        Base::bufferMutex.lock();
      }
      Base::maxChunk = chunkSize;
      updateChunks();
    }
    Base::chunks.push_back(move(chunk));

    return true;
  }

  template<typename T>
  void
  SasAnalysis_Write<T>::dumpChunk(const vector<vector<SasAtom>>& chunk,
                         Serializer<MetaStream<T>>& out) const
  {
    unsigned int size = chunk.size();
    out << size;

    for(const vector<SasAtom>& frame : chunk)
    {
      if(Base::gromacs and Base::gromacs->isAborting())
        break;
      for(const SasAtom& sasAtom : frame)
        out << sasAtom;
    }
  }

  template<typename T>
  std::vector<vector<SasAtom>>
  SasAnalysis_Read<T>::loadChunk(Serializer<MetaStream<T>>& in)
  {
    unsigned int size;
    std::vector<vector<SasAtom>> chunk;

    in >> size;
    chunk.reserve(size);

    for(unsigned int i = 0; i < size; i++)
    {
      vector<SasAtom> atoms;
      atoms.reserve(Base::nAtoms);

      auto iter = std::begin(atoms);
      for(unsigned int k = 0; k < Base::nAtoms; k++, iter++)
      {
        if(Base::gromacs and Base::gromacs->isAborting())
          return chunk;
        try
        {
          SasAtom atom;
          in >> atom;
          atoms.push_back(move(atom));
        }
        catch(int e)
        {
          std::cerr << "Warning: inconsistent binary file." << std::endl;
          break;
        }
      }

      chunk.push_back(move(atoms));
    }

    return chunk;
  }

  template<typename T>
  bool
  SasAnalysis_Base<T>::setMaxBytes(unsigned long bytes)
  {
    if(not changeable)
      return false;

    maxBytes = bytes;
    updateChunks();

    return true;
  }

  template<typename T>
  unsigned long
  SasAnalysis_Base<T>::getMaxBytes()
  {
    return maxBytes;
  }

  template<typename T>
  bool
  SasAnalysis_Base<T>::setMaxChunkSize(unsigned long bytes)
  {
    if(not changeable)
      return false;

    maxChunk = bytes;
    updateChunks();

    return true;
  }

  template<typename T>
  unsigned long
  SasAnalysis_Base<T>::getMaxChunkSize()
  {
    return maxChunk;
  }

  template<typename T>
  void
  SasAnalysis_Base<T>::updateChunks()
  {
    // SasAtom serialization produces real * 4
    maxFrames = maxChunk / (nAtoms * sizeof(real) * 4);
    vector<vector<SasAtom>> testVector(maxFrames);
    // Let's check the true allocation
    unsigned int vectorSize = testVector.capacity()
        * (sizeof(vector<SasAtom> ) + sizeof(SasAtom) * nAtoms)
                              + sizeof(vector<vector<SasAtom>> );

    bufferMax = maxBytes / vectorSize;
    if(bufferMax == 0)
    {
      cerr << "Can't allocate memory... Strange. Is your RAM full?" << endl;
      throw bad_alloc();
    }

    chunks = boost::circular_buffer<std::vector<vector<SasAtom>>>(bufferMax);
  }

  template<typename T>
  void
  SasAnalysis_Read<T>::updateChunks()
  {
    Base::updateChunks();
    Base::bufferCount = 0;
  }

  template<typename T>
  void
  SasAnalysis_Write<T>::updateChunks()
  {
    Base::updateChunks();
    Base::bufferCount = Base::bufferMax - 1;
  }

  template<typename T>
  unsigned int
  SasAnalysis_Write<T>::getReadFrames() const
  {
    return readFrames;
  }
}

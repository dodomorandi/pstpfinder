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

#ifndef _SASANALYSIS_H
#define _SASANALYSIS_H

namespace PstpFinder
{
  template<typename T> class SasAnalysis_Base;
  template<typename T> class SasAnalysis_Read;
  template<typename T> class SasAnalysis_Write;
  template<typename T, typename = void> class SasAnalysis;
}

#include "SasAtom.h"
#include "Gromacs.h"
#include "Session.h"
#include "SasAnalysisThread.h"
#include "utils.h"
#include "Serializer.h"

#include <thread>
#include <sys/sysinfo.h>
#include <boost/circular_buffer.hpp>

#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <type_traits>

namespace PstpFinder
{
  template<typename T, typename>
  class SasAnalysisThread;

  template<typename T>
  class SasAnalysis_Base
  {
    public:
      SasAnalysis_Base(unsigned int nAtoms, const Gromacs&, Session<T>&);
      SasAnalysis_Base(const Gromacs& gromacs,
                       const std::string& sessionFileName);
      SasAnalysis_Base(const Gromacs&, Session<T>&);
      virtual bool setMaxBytes(unsigned long bytes);
      virtual unsigned long getMaxBytes();
      virtual bool setMaxChunkSize(unsigned long bytes);
      virtual unsigned long getMaxChunkSize();

    protected:
      boost::circular_buffer<std::vector<SasAtom*>> chunks;
      std::vector<SasAtom*> frames;
      unsigned int nAtoms;
      Session<T> rawSession;
      MetaStream<T>& sasMetaStream;
      Serializer<MetaStream<T>>* serializer;
      std::streampos fileStreamEnd;
      const Gromacs* gromacs;
      unsigned long maxFrames, maxBytes, maxChunk;
      unsigned int bufferCount;
      unsigned int bufferMax;
      mutable std::condition_variable bufferCountCondition;
      mutable std::mutex bufferMutex, bufferCountMutex;
      bool changeable;

      template<typename, typename> friend class SasAnalysisThread_Base;
      template<typename, typename> friend class SasAnalysisThread;
      typedef SasAnalysisThread<T> SasAnalysisThreadType;
      SasAnalysisThreadType* analysisThread;

      virtual void init();
      virtual void updateChunks();
  };

  template<typename T>
  class SasAnalysis_Read : public SasAnalysis_Base<T>
  {
    public:
      typedef SasAnalysisThread<T> SasAnalysisThreadType;
      SasAnalysis_Read(unsigned int nAtoms, const Gromacs& gromacs,
                        Session<T>& sessionFile) :
          Base(nAtoms, gromacs, sessionFile) { updateChunks(); }
      SasAnalysis_Read(const Gromacs& gromacs,
                       const std::string& sessionFileName) :
          Base(gromacs, sessionFileName) { updateChunks(); }
      SasAnalysis_Read(const Gromacs& gromacs, Session<T>& sessionFile) :
          Base(gromacs, sessionFile) { updateChunks(); }
      virtual ~SasAnalysis_Read();
      virtual bool read(std::vector<SasAtom>& sasAtom);

    private:
      typedef SasAnalysis_Base<T> Base;
      template<typename, typename> friend class SasAnalysisThread_Base;
      template<typename, typename> friend class SasAnalysisThread;

      virtual std::vector<SasAtom*>
        loadChunk(Serializer<MetaStream<T>>& in);
      bool open();
      virtual void updateChunks();
  };


  template<typename T>
  class SasAnalysis_Write : public SasAnalysis_Base<T>
  {
    public:
      typedef SasAnalysisThread<T> SasAnalysisThreadType;
      SasAnalysis_Write(unsigned int nAtoms, const Gromacs& gromacs,
                        Session<T>& sessionFile) :
          Base(nAtoms, gromacs, sessionFile), readFrames(0) { updateChunks(); }
      SasAnalysis_Write(const Gromacs& gromacs,
                        const std::string& sessionFileName) :
          Base(gromacs, sessionFileName), readFrames(0) { updateChunks(); }
      SasAnalysis_Write(const Gromacs& gromacs, Session<T>& sessionFile) :
          Base(gromacs, sessionFile), readFrames(0) { updateChunks(); }
      virtual ~SasAnalysis_Write();
      virtual void write(const std::vector<SasAtom>& sasAtoms);
      unsigned int getReadFrames() const;

    protected:
      unsigned int readFrames;

    private:
      typedef SasAnalysis_Base<T> Base;
      template<typename, typename> friend class SasAnalysisThread_Base;
      template<typename, typename> friend class SasAnalysisThread;

      virtual void dumpChunk(const std::vector<SasAtom*>& chunk,
                Serializer<MetaStream<T>>& out) const;
      virtual void flush();
      virtual bool save();
      virtual void updateChunks();
  };

  template<typename T>
  class SasAnalysis<T,
      typename std::enable_if<
          is_stream_base_of<std::basic_istream, T>::value
          and not is_stream_base_of<std::basic_ostream, T>
        ::value>::type> :
      public SasAnalysis_Read<T>
  {
    public:
      SasAnalysis(unsigned int nAtoms, const Gromacs& gromacs,
                  Session<T>& sessionFile) :
          SasAnalysis_Read<T>(nAtoms, gromacs, sessionFile) {}
      SasAnalysis(const Gromacs& gromacs, Session<T>& sessionFile) :
          SasAnalysis_Read<T>(gromacs, sessionFile) {}
      SasAnalysis(const Gromacs& gromacs,
                  const std::string& sessionFileName) :
          SasAnalysis_Read<T>(gromacs, sessionFileName) {}
  };

  template<typename T>
  class SasAnalysis<T,
      typename std::enable_if<
        not is_stream_base_of<std::basic_istream, T>::value and
        is_stream_base_of<std::basic_ostream, T>::value>::type> :
      public SasAnalysis_Write<T>
  {
    public:
      SasAnalysis(unsigned int nAtoms, const Gromacs& gromacs,
                  Session<T>& sessionFile) :
          SasAnalysis_Write<T>(nAtoms, gromacs, sessionFile) {}
      SasAnalysis(const Gromacs& gromacs, Session<T>& sessionFile) :
          SasAnalysis_Write<T>(gromacs, sessionFile) {}
      SasAnalysis(const Gromacs& gromacs,
                  const std::string& sessionFileName) :
          SasAnalysis_Write<T>(gromacs, sessionFileName) {}
  };

  template<typename T>
  class SasAnalysis<T,
      typename std::enable_if<
      is_stream_base_of<std::basic_istream, T>::value and
      is_stream_base_of<std::basic_ostream, T>::value>::type> :
      public SasAnalysis_Write<T>
  {
    public:
      SasAnalysis(unsigned int nAtoms, const Gromacs& gromacs,
                  Session<T>& sessionFile) :
          SasAnalysis_Write<T>(nAtoms, gromacs, sessionFile) { init(); }
      SasAnalysis(const Gromacs& gromacs, Session<T>& sessionFile) :
          SasAnalysis_Write<T>(gromacs, sessionFile) { init(); }
      SasAnalysis(const Gromacs& gromacs,
                  const std::string& sessionFileName) :
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
          std::streampos backupPosition;

          *Base::serializer >> chunkSize;
          backupPosition = Base::sasMetaStream.tellg();

          while(not Base::sasMetaStream.eof())
          {
            Base::sasMetaStream.seekg(chunkSize * sasAtomSize * Base::nAtoms,
                                      std::ios_base::cur);

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
            Base::sasMetaStream.seekg(0, std::ios_base::end);

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
                                        const std::string& sessionFileName) :
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

    for(auto& frame : Base::frames)
      delete[] frame;

    delete Base::serializer;
  }

  template<typename T>
  SasAnalysis_Write<T>::~SasAnalysis_Write()
  {
    if(Base::changeable)
      return;

    /* Do not save the last chunk if we are aborting! */
    if(not Base::gromacs or Base::gromacs->isAborting())
    {
      Base::bufferMutex.lock();
      Base::frames.clear();
      Base::bufferMutex.unlock();
    }
    else if(Base::frames.size() != 0)
      flush();

    Base::bufferMutex.lock();
    while(Base::bufferCount != Base::bufferMax - 1)
    {
      std::unique_lock<std::mutex> bufferCountLock(Base::bufferCountMutex);
      Base::bufferMutex.unlock();
      Base::bufferCountCondition.wait(bufferCountLock);
      Base::bufferMutex.lock();
    }
    Base::bufferMutex.unlock();

    Base::analysisThread->stop();
    delete Base::analysisThread;

    for(auto& chunk : Base::chunks)
    {
      for(auto& frame : chunk)
        delete[] frame;
    }

    delete Base::serializer;
    Base::sasMetaStream.close();
  }

  template<typename T>
  void
  SasAnalysis_Write<T>::write(const std::vector<SasAtom>& sasAtoms)
  {
    assert(sasAtoms.size() == Base::nAtoms);
    SasAtom* tmpFrame = new SasAtom[Base::nAtoms];
    if(Base::changeable)
    {
      Base::changeable = false;
      Base::serializer = new Serializer<MetaStream<T>>(Base::sasMetaStream);

      Base::analysisThread = new SasAnalysisThreadType(*this);
    }

    std::copy(std::begin(sasAtoms), std::end(sasAtoms), tmpFrame);
    Base::bufferMutex.lock();
    Base::frames.push_back(tmpFrame);

    bool cond = (Base::frames.size() == Base::maxFrames);
    Base::bufferMutex.unlock();

    if(cond)
      flush();
  }

  template<typename T>
  bool
  SasAnalysis_Read<T>::read(std::vector<SasAtom>& sasAtom)
  {
    static std::vector<SasAtom*>::const_iterator currentFrameIter;

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
        std::unique_lock<std::mutex> bufferCountLock(Base::bufferCountMutex);
        Base::bufferMutex.unlock();
        Base::bufferCountCondition.wait(bufferCountLock);
        Base::bufferMutex.lock();
      }
      Base::frames = Base::chunks.front();
      currentFrameIter = Base::frames.begin();

      Base::bufferCount--;
      Base::bufferMutex.unlock();
    }
    else if(currentFrameIter == Base::frames.end())
    {
      Base::bufferMutex.lock();

      for(SasAtom* frame : Base::frames)
        delete[] frame;
      Base::frames.clear();

      Base::chunks.pop_front();

      if(Base::bufferCount == 0)
      {
        Base::bufferMutex.unlock();
        std::unique_lock<std::mutex> bufferCountLock(Base::bufferCountMutex);
        Base::bufferCountCondition.wait(bufferCountLock);
        Base::bufferMutex.lock();
      }

      Base::bufferCount--;

      if(Base::chunks.size() == 0)
      {
        Base::bufferMutex.unlock();
        sasAtom.clear();
        return false;
      }

      Base::frames = Base::chunks.front();
      currentFrameIter = Base::frames.begin();

      Base::analysisThread->wakeUp();
      Base::bufferMutex.unlock();
    }

    if(sasAtom.size() != Base::nAtoms)
    {
      sasAtom.clear();
      sasAtom.resize(Base::nAtoms);
    } 
    std::copy_n(*currentFrameIter, Base::nAtoms, std::begin(sasAtom));
    ++currentFrameIter;

    return true;
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
      std::unique_lock<std::mutex> bufferCountLock(Base::bufferCountMutex);
      Base::bufferMutex.unlock();
      Base::bufferCountCondition.wait(bufferCountLock);
      Base::bufferMutex.lock();
    }
    Base::bufferCount--;
    Base::chunks.push_back(Base::frames);

    Base::frames.clear();
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

    std::vector<SasAtom*> chunk = loadChunk(*Base::serializer);
    unsigned long chunkSize = chunk.capacity()
                              * (sizeof(SasAtom*)
                                 + sizeof(SasAtom) * Base::nAtoms)
                              + sizeof(std::vector<SasAtom*> );
    if(chunkSize > Base::maxChunk)
    {
      // bufferMutex already lock_upgrade from ThreadOpen
      while(Base::bufferCount != 0)
      {
        std::unique_lock<std::mutex> lock(Base::bufferCountMutex);
        Base::bufferMutex.unlock();
        Base::bufferCountCondition.wait(lock);
        Base::bufferMutex.lock();
      }
      Base::maxChunk = chunkSize;
      updateChunks();
    }
    Base::chunks.push_back(chunk);

    return true;
  }

  template<typename T>
  void
  SasAnalysis_Write<T>::dumpChunk(const std::vector<SasAtom*>& chunk,
                         Serializer<MetaStream<T>>& out) const
  {
    unsigned int size = chunk.size();
    out << size;

    for(const SasAtom* frame : chunk)
    {
      if(Base::gromacs and Base::gromacs->isAborting())
        break;
      const SasAtom* end = frame + Base::nAtoms;
      for(const SasAtom* atom = frame; atom < end; ++atom)
        out << *atom;
    }
  }

  template<typename T>
  std::vector<SasAtom*>
  SasAnalysis_Read<T>::loadChunk(Serializer<MetaStream<T>>& in)
  {
    unsigned int size;
    std::vector<SasAtom*> chunk;
    SasAtom* atoms;
    SasAtom* atom;

    in >> size;
    chunk.reserve(size);

    for(unsigned int i = 0; i < size; i++)
    {
      atoms = new SasAtom[Base::nAtoms];
      atom = atoms;

      for(unsigned int k = 0; k < Base::nAtoms; k++, atom++)
      {
        if(Base::gromacs and Base::gromacs->isAborting())
          return chunk;
        try
        {
          in >> *atom;
        }
        catch(int e)
        {
          std::cerr << "Warning: inconsistent binary file." << std::endl;
          SasAtom* tmpatoms = new SasAtom[k];
          std::copy(atoms, atoms + k, tmpatoms);
          delete[] atoms;
          atoms = tmpatoms;
          break;
        }
      }

      chunk.push_back(atoms);
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
    std::vector<SasAtom*> testVector(maxFrames);
    unsigned int vectorSize = testVector.capacity()
                              * (sizeof(SasAtom*) + sizeof(SasAtom) * nAtoms)
                              + sizeof(std::vector<SasAtom*> );

    bufferMax = maxBytes / vectorSize;
    if(bufferMax == 0)
    {
      std::cerr << "Can't allocate memory... Strange. Is your RAM full?" <<
          std::endl;
      throw std::bad_alloc();
    }

    chunks = boost::circular_buffer<std::vector<SasAtom*> >(bufferMax);
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
#endif


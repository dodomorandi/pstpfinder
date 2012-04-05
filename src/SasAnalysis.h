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

#include "SasAtom.h"
#include "Gromacs.h"
#include "Session.h"
#include "SasAnalysisThread.h"

#include <thread>
#include <sys/sysinfo.h>
#include <boost/circular_buffer.hpp>
#include <boost/filesystem.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <vector>
#include <string>
#include <iostream>
#include <sstream>

namespace archive = boost::archive;
namespace fs = boost::filesystem;

namespace PstpFinder
{
  template<typename T, typename>
  class SasAnalysisThread;

  template<typename T>
  class SasAnalysis
  {
    public:
      SasAnalysis(unsigned int nAtoms, const Gromacs& gromacs,
                  Session<T>& sessionFile,
                  bool savingMode = true);

      SasAnalysis(const Gromacs& gromacs, Session<T>& sessionFile,
                  bool savingMode = true);

      SasAnalysis(const Gromacs& gromacs, const string& sessionFileName,
                  bool savingMode = true);

      ~SasAnalysis();

      const SasAnalysis&
      operator <<(SasAtom* sasAtoms);

      SasAnalysis&
      operator >>(SasAtom*& sasAtom);

      bool
      setMaxBytes(unsigned long bytes);

      unsigned long
      getMaxBytes();

      bool
      setMaxChunkSize(unsigned long bytes);

      unsigned long
      getMaxChunkSize();

    private:
      boost::circular_buffer<std::vector<SasAtom*> > chunks;
      std::vector<SasAtom*> frames;
      unsigned int nAtoms;
      Session<T> rawSession;
      MetaStream<T>& sasMetaStream;
      std::streampos fileStreamEnd;
      const Gromacs* gromacs;
      unsigned long maxFrames, maxBytes, maxChunk;
      unsigned int bufferCount;
      unsigned int bufferMax;
      mutable condition_variable bufferCountCondition;
      mutable mutex bufferMutex, bufferCountMutex;
      bool changeable;
      boost::archive::binary_iarchive* inArchive;
      boost::archive::binary_oarchive* outArchive;

      enum
      {
        MODE_OPEN,
        MODE_SAVE
      } mode;

      template<typename> friend class SasAnalysisThread_Base;
      template<typename, typename> friend class SasAnalysisThread;
      typedef SasAnalysisThread<T> SasAnalysisThreadType;
      SasAnalysisThreadType* analysisThread;

      void
      init(bool savingMode = true);

      void
      dumpChunk(const std::vector<SasAtom*>& chunk,
                boost::archive::binary_oarchive& out) const;

      std::vector<SasAtom*>
      loadChunk(boost::archive::binary_iarchive& in);

      void
      flush();

      bool
      save(const std::string& filename);

      bool
      save();

      bool
      open();

      void
      updateChunks();
  };

  template<typename T>
  SasAnalysis<T>::SasAnalysis(unsigned int nAtoms, const Gromacs& gromacs,
                           Session<T>& session, bool savingMode) :
      sasMetaStream(session.getSasStream())
  {
    this->nAtoms = nAtoms;
    this->gromacs = &gromacs;
    init(savingMode);
  }

  template<typename T>
  SasAnalysis<T>::SasAnalysis(const Gromacs& gromacs, Session<T>& session,
                           bool savingMode) :
      sasMetaStream(session.getSasStream())
  {
    nAtoms = gromacs.getGroup("Protein").size();
    this->gromacs = &gromacs;
    init(savingMode);
  }

  template<typename T>
  SasAnalysis<T>::SasAnalysis(const Gromacs& gromacs,
                              const string& sessionFileName, bool savingMode) :
      rawSession(sessionFileName), sasMetaStream(rawSession.getSasStream())
  {
    nAtoms = gromacs.getGroup("Protein").size();
    this->gromacs = &gromacs;
    init(savingMode);
  }

  template<typename T>
  void
  SasAnalysis<T>::init(bool savingMode)
  {
    changeable = true;
    if(savingMode)
      mode = MODE_SAVE;
    else
      mode = MODE_OPEN;

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

    updateChunks();
  }

  template<typename T>
  SasAnalysis<T>::~SasAnalysis()
  {
    if(changeable)
      return;

    if(mode == MODE_SAVE)
    {
      if(frames.size() != 0)
        flush();
      analysisThread->stop();
      delete analysisThread;

      delete outArchive;
    }
    else
    {
      analysisThread->stop();
      delete analysisThread;

      delete inArchive;
    }
  }

  template<typename T>
  const SasAnalysis<T>&
  SasAnalysis<T>::operator <<(SasAtom* sasAtoms)
  {
    SasAtom* tmpFrame = new SasAtom[nAtoms];
    if(changeable)
    {
      changeable = false;
      outArchive = new archive::binary_oarchive(sasMetaStream);

      analysisThread = new SasAnalysisThreadType(*this);
    }

    std::copy(sasAtoms, sasAtoms + nAtoms, tmpFrame);
    bufferMutex.lock();
    frames.push_back(tmpFrame);

    bool cond = (frames.size() == maxFrames);
    bufferMutex.unlock();

    if(cond)
      flush();
    return *this;
  }

  template<typename T>
  SasAnalysis<T>&
  SasAnalysis<T>::operator >>(SasAtom*& sasAtom)
  {
    static std::vector<SasAtom*>::const_iterator i;

    if(changeable)
    {
      changeable = false;
      inArchive = new archive::binary_iarchive(sasMetaStream);

      analysisThread = new SasAnalysisThreadType(*this);
    }

    if(frames.size() == 0)
    {
      bufferMutex.lock();
      while(bufferCount == 0)
      {
        unique_lock<mutex> bufferCountLock(bufferCountMutex);
        bufferMutex.unlock();
        bufferCountCondition.wait(bufferCountLock);
        bufferMutex.lock();
      }
      frames = chunks.front();
      i = frames.begin();

      bufferCount--;
      bufferMutex.unlock();
    }
    else if(i == frames.end())
    {
      bufferMutex.lock();

      for(i = frames.begin(); i < frames.end(); i++)
        delete[] *i;
      frames.clear();

      chunks.pop_front();

      if(bufferCount == 0)
      {
        bufferMutex.unlock();
        unique_lock<mutex> bufferCountLock(bufferCountMutex);
        bufferCountCondition.wait(bufferCountLock);
        bufferMutex.lock();
      }

      bufferCount--;

      if(chunks.size() == 0)
      {
        bufferMutex.unlock();
        sasAtom = 0;
        return *this;
      }

      frames = chunks.front();
      i = frames.begin();

      analysisThread->wakeUp();
      bufferMutex.unlock();
    }

    sasAtom = *(i++);

    return *this;
  }

  template<typename T>
  void
  SasAnalysis<T>::flush()
  {
    if(mode != MODE_SAVE)
      return;

    if(changeable)
    {
      changeable = false;
      outArchive = new archive::binary_oarchive(sasMetaStream);

      analysisThread = new SasAnalysisThreadType(*this);
    }

    bufferMutex.lock();
    while(bufferCount == 0)
    {
      unique_lock<mutex> bufferCountLock(bufferCountMutex);
      bufferMutex.unlock();
      bufferCountCondition.wait(bufferCountLock);
      bufferMutex.lock();
    }
    bufferCount--;
    chunks.push_back(frames);

    frames.clear();
    analysisThread->wakeUp();
    bufferMutex.unlock();

  }

  template<typename T>
  bool
  SasAnalysis<T>::save()
  {
    if(mode != MODE_SAVE)
      return false;

    if(chunks.empty())
      return false;

    dumpChunk(chunks.front(), *outArchive);

    return true;
  }

  template<typename T>
  bool
  SasAnalysis<T>::open()
  {
    if(mode != MODE_OPEN)
      return false;

    sasMetaStream.peek();
    if(sasMetaStream.eof())
      return false;

    vector<SasAtom*> chunk = loadChunk(*inArchive);
    unsigned long chunkSize = chunk.capacity()
                              * (sizeof(SasAtom*) + sizeof(SasAtom) * nAtoms)
                              + sizeof(vector<SasAtom*>);
    if(chunkSize > maxChunk)
    {
      // bufferMutex already lock_upgrade from ThreadOpen
      while(bufferCount != 0)
      {
        unique_lock<mutex> lock(bufferCountMutex);
        bufferMutex.unlock();
        bufferCountCondition.wait(lock);
        bufferMutex.lock();
      }
      maxChunk = chunkSize;
      updateChunks();
    }
    chunks.push_back(chunk);

    return true;
  }

  template<typename T>
  void
  SasAnalysis<T>::dumpChunk(const std::vector<SasAtom*>& chunk,
                         archive::binary_oarchive& out) const
  {
    unsigned int size = chunk.size();
    out << size;

    for(std::vector<SasAtom*>::const_iterator i = chunk.begin();
        i < chunk.end(); i++)
    {
      if(gromacs and gromacs->isAborting())
        break;
      const SasAtom* end = *i + nAtoms;
      for(SasAtom* j = *i; j < end; j++)
        out << *j;
    }
  }

  template<typename T>
  std::vector<SasAtom*>
  SasAnalysis<T>::loadChunk(archive::binary_iarchive& in)
  {
    unsigned int size;
    std::vector<SasAtom*> chunk;
    SasAtom* atoms;
    SasAtom* atom;

    in >> size;
    chunk.reserve(size);

    for(unsigned int i = 0; i < size; i++)
    {
      atoms = new SasAtom[nAtoms];
      atom = atoms;

      for(unsigned int k = 0; k < nAtoms; k++, atom++)
      {
        if(gromacs and gromacs->isAborting())
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
  SasAnalysis<T>::setMaxBytes(unsigned long bytes)
  {
    if(not changeable)
      return false;

    maxBytes = bytes;
    updateChunks();

    return true;
  }

  template<typename T>
  unsigned long
  SasAnalysis<T>::getMaxBytes()
  {
    return maxBytes;
  }

  template<typename T>
  bool
  SasAnalysis<T>::setMaxChunkSize(unsigned long bytes)
  {
    if(not changeable)
      return false;

    maxChunk = bytes;
    updateChunks();

    return true;
  }

  template<typename T>
  unsigned long
  SasAnalysis<T>::getMaxChunkSize()
  {
    return maxChunk;
  }

  template<typename T>
  void
  SasAnalysis<T>::updateChunks()
  {
    // SasAtom serialization produces real * 4
    maxFrames = maxChunk / (nAtoms * sizeof(real) * 4);
    vector<SasAtom*> testVector(maxFrames);
    unsigned int vectorSize = testVector.capacity()
                              * (sizeof(SasAtom*) + sizeof(SasAtom) * nAtoms)
                              + sizeof(vector<SasAtom*> );

    bufferMax = maxBytes / vectorSize;

    if(mode == MODE_SAVE)
      bufferCount = bufferMax - 1;
    else
      bufferCount = 0;

    chunks = boost::circular_buffer<std::vector<SasAtom*> >(bufferMax);
  }
}
#endif


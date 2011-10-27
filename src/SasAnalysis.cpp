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

namespace archive = boost::archive;
namespace io = boost::iostreams;
namespace fs = boost::filesystem;

namespace PstpFinder
{

  SasAnalysis::SasAnalysis(unsigned int nAtoms, const Gromacs& gromacs,
                           std::string filename, bool savingMode)
  {
    this->nAtoms = nAtoms;
    this->gromacs = &gromacs;
    init(filename, savingMode);
  }

  SasAnalysis::SasAnalysis(const Gromacs& gromacs, std::string filename,
                           bool savingMode)
  {
    nAtoms = gromacs.getGroup("Protein").size();
    this->gromacs = &gromacs;
    init(filename, savingMode);
  }

  void
  SasAnalysis::init(string filename, bool savingMode)
  {
    changeable = true;
    this->filename = filename;
    if(savingMode)
    {
      fileIO = io::file_descriptor(filename, BOOST_IOS::trunc |
                                             BOOST_IOS::out |
                                             BOOST_IOS::binary);
      mode = MODE_SAVE;
    }
    else
    {
      sessionFile = Session(filename);
      mode = MODE_OPEN;
    }

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

  SasAnalysis::~SasAnalysis()
  {
    if(changeable)
      return;

    if(mode == MODE_SAVE)
    {
      if(frames.size() != 0)
        flush();
      operationThread->stop();
      delete operationThread;

      delete outArchive;
    }
    else
    {
      operationThread->stop();
      delete operationThread;

      delete inArchive;
    }
  }

  const SasAnalysis&
  SasAnalysis::operator <<(SasAtom* sasAtoms)
  {
    SasAtom* tmpFrame = new SasAtom[nAtoms];
    if(changeable)
    {
      changeable = false;

      //outFilter.push(io::zlib_compressor());
      outFilter.push(fileIO);
      outArchive = new archive::binary_oarchive(outFilter);

      operationThread = new OperationThread(*this);
    }

    std::copy(sasAtoms, sasAtoms + nAtoms, tmpFrame);
    bufferMutex.lock();
    frames.push_back(tmpFrame);
    bufferMutex.unlock_and_lock_shared();

    bool cond = (frames.size() == maxFrames);
    bufferMutex.unlock_shared();

    if(cond)
      flush();
    return *this;
  }

  SasAnalysis&
  SasAnalysis::operator >>(SasAtom*& sasAtom)
  {
    static std::vector<SasAtom*>::const_iterator i;

    if(changeable)
    {
      changeable = false;

      //inFilter.push(io::zlib_decompressor());
      inFilter.push(sessionFile.getSasStream());
      inArchive = new archive::binary_iarchive(inFilter);

      operationThread = new OperationThread(*this);
    }

    if(frames.size() == 0)
    {
      bufferMutex.lock_upgrade();
      while(bufferCount == 0)
      {
        boost::unique_lock<boost::mutex> bufferCountLock(bufferCountMutex);
        bufferMutex.unlock_upgrade();
        bufferCountCondition.wait(bufferCountLock);
        bufferMutex.lock_upgrade();
      }
      bufferMutex.unlock_upgrade_and_lock();
      frames = chunks.front();
      i = frames.begin();

      bufferCount--;
      bufferMutex.unlock();
    }
    else if(i == frames.end())
    {
      bufferMutex.lock_upgrade();

      for(i = frames.begin(); i < frames.end(); i++)
        delete[] *i;
      frames.clear();

      bufferMutex.unlock_upgrade_and_lock();
      chunks.pop_front();

      if(bufferCount == 0)
      {
        bufferMutex.unlock();
        boost::unique_lock<boost::mutex> bufferCountLock(bufferCountMutex);
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

      operationThread->wakeUp();
      bufferMutex.unlock();
    }

    sasAtom = *(i++);

    return *this;
  }

  void
  SasAnalysis::flush()
  {
    if(mode != MODE_SAVE)
      return;

    if(changeable)
    {
      changeable = false;

      //outFilter.push(io::zlib_compressor());
      outFilter.push(fileIO);
      outArchive = new archive::binary_oarchive(outFilter);

      operationThread = new OperationThread(*this);
    }

    bufferMutex.lock_upgrade();
    while(bufferCount == 0)
    {
      boost::unique_lock<boost::mutex> bufferCountLock(bufferCountMutex);
      bufferMutex.unlock_upgrade();
      bufferCountCondition.wait(bufferCountLock);
      bufferMutex.lock_upgrade();
    }
    bufferMutex.unlock_upgrade_and_lock();
    bufferCount--;
    chunks.push_back(frames);

    frames.clear();
    operationThread->wakeUp();
    bufferMutex.unlock();

  }

  bool
  SasAnalysis::save(const std::string& filename)
  {
    if(mode != MODE_SAVE)
      return false;

    if(this->filename == filename)
      return true;

    std::string old_filename = this->filename;
    fs::copy_file(old_filename, filename);
    if(not fs::exists(filename))
      return false;

    this->filename = filename;
    fileIO = io::file_descriptor(filename, BOOST_IOS::trunc | BOOST_IOS::out |
    BOOST_IOS::binary);

    if(not save())
    {
      fileIO.close();
      boost::filesystem::remove(filename);
      this->filename = old_filename;
      fileIO = io::file_descriptor(old_filename, BOOST_IOS::app | BOOST_IOS::out |
      BOOST_IOS::binary);
      return false;
    }

    return true;
  }

  bool
  SasAnalysis::save()
  {
    if(mode != MODE_SAVE)
      return false;

    if(chunks.empty())
      return false;

    dumpChunk(chunks.front(), *outArchive);

    return true;
  }

  bool
  SasAnalysis::open()
  {
    if(mode != MODE_OPEN)
      return false;

    inFilter.peek();
    if(inFilter.eof())
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
        boost::unique_lock<boost::mutex> lock(bufferCountMutex);
        bufferMutex.unlock_upgrade();
        bufferCountCondition.wait(lock);
        bufferMutex.lock_upgrade();
      }
      bufferMutex.unlock_upgrade_and_lock();
      maxChunk = chunkSize;
      updateChunks();
      bufferMutex.unlock_and_lock_upgrade();
    }
    chunks.push_back(chunk);

    return true;
  }

  void
  SasAnalysis::dumpChunk(const std::vector<SasAtom*>& chunk,
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

  std::vector<SasAtom*>
  SasAnalysis::loadChunk(archive::binary_iarchive& in)
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

  SasAnalysis::OperationThread::OperationThread(SasAnalysis& parent)
  {
    this->parent = &parent;
    isStopped = false;
    if(parent.mode == SasAnalysis::MODE_SAVE)
      thread = boost::thread(
          boost::bind(&SasAnalysis::OperationThread::threadSave,
                      boost::ref(*this)));
    else
      thread = boost::thread(
          boost::bind(&SasAnalysis::OperationThread::threadOpen,
                      boost::ref(*this)));
  }

  SasAnalysis::OperationThread::~OperationThread()
  {
    stop();
  }

  void
  SasAnalysis::OperationThread::threadSave()
  {
    while(not isStopped)
    {
      parent->bufferMutex.lock_upgrade();

      while(parent->bufferCount == parent->bufferMax - 1 and not isStopped)
      {
        boost::unique_lock<boost::mutex> lock(wakeMutex);
        parent->bufferMutex.unlock_upgrade();
        wakeCondition.wait(lock);
        parent->bufferMutex.lock_upgrade();
      }

      if(isStopped or (parent->gromacs and parent->gromacs->isAborting()))
      {
        parent->bufferMutex.unlock_upgrade();
        break;
      }

      if(parent->save())
      {
        parent->bufferMutex.unlock_upgrade_and_lock();
        vector<SasAtom*>& curChunk = parent->chunks.front();
        for(std::vector<SasAtom*>::iterator i = curChunk.begin();
            i < curChunk.end(); i++)
          delete[] *i;
        curChunk.clear();

        parent->chunks.pop_front();
        parent->bufferCount++;
        parent->bufferCountCondition.notify_all();

        parent->bufferMutex.unlock_and_lock_upgrade();
      }

      parent->bufferMutex.unlock_upgrade();
    }

    if(parent->gromacs and parent->gromacs->isAborting())
      return;

    parent->bufferMutex.lock_upgrade();
    while(parent->bufferCount < parent->bufferMax - 1)
    {
      parent->bufferMutex.unlock_upgrade_and_lock();
      parent->save();
      parent->chunks.pop_front();
      parent->bufferCount++;
      parent->bufferCountCondition.notify_all();
      parent->bufferMutex.unlock_and_lock_upgrade();
    }
    parent->bufferMutex.unlock_upgrade();
  }

  void
  SasAnalysis::OperationThread::threadOpen()
  {
    while(not isStopped)
    {
      parent->bufferMutex.lock_upgrade();

      while(parent->bufferCount == parent->bufferMax - 1 and not isStopped)
      {
        boost::unique_lock<boost::mutex> lock(wakeMutex);
        parent->bufferMutex.unlock_upgrade();
        wakeCondition.wait(lock);
        parent->bufferMutex.lock_upgrade();
      }

      if(isStopped)
      {
        while(parent->bufferCount < parent->bufferMax - 1)
        {
          parent->bufferMutex.unlock_upgrade_and_lock();
          parent->bufferCount++;
          parent->bufferCountCondition.notify_all();
          parent->bufferMutex.unlock_and_lock_upgrade();
        }
        parent->bufferMutex.unlock_upgrade();
        break;
      }

      if(parent->gromacs and parent->gromacs->isAborting())
      {
        parent->bufferMutex.unlock_upgrade();
        break;
      }

      if(not parent->open())
        isStopped = true;

      parent->bufferMutex.unlock_upgrade_and_lock();
      parent->bufferCount++;
      parent->bufferCountCondition.notify_all();
      parent->bufferMutex.unlock();
    }
  }

  void
  SasAnalysis::OperationThread::wakeUp()
  {
    wakeCondition.notify_one();
  }

  void
  SasAnalysis::OperationThread::stop()
  {
    if(isStopped)
    {
      thread.join();
      return;
    }

    parent->bufferMutex.lock();
    isStopped = true;
    wakeCondition.notify_one();
    parent->bufferMutex.unlock();

    thread.join();
  }

  bool
  SasAnalysis::setMaxBytes(unsigned long bytes)
  {
    if(not changeable)
      return false;

    maxBytes = bytes;
    updateChunks();

    return true;
  }

  unsigned long
  SasAnalysis::getMaxBytes()
  {
    return maxBytes;
  }

  bool
  SasAnalysis::setMaxChunkSize(unsigned long bytes)
  {
    if(not changeable)
      return false;

    maxChunk = bytes;
    updateChunks();

    return true;
  }

  unsigned long
  SasAnalysis::getMaxChunkSize()
  {
    return maxChunk;
  }

  void
  SasAnalysis::updateChunks()
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

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

using namespace Gromacs;
namespace archive = boost::archive;
namespace io = boost::iostreams;
namespace ip = boost::interprocess;
namespace fs = boost::filesystem;

SasAnalysis::SasAnalysis(unsigned int nAtoms,
                         std::string filename,
                         bool savingMode)
{
  this->nAtoms = nAtoms;
  gromacs = 0;
  init(filename, savingMode);
}

SasAnalysis::SasAnalysis(const Gromacs& gromacs,
                         std::string filename,
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
    fileIO = io::file_descriptor(filename, BOOST_IOS::trunc | BOOST_IOS::out |
                                           BOOST_IOS::binary);
    mode = MODE_SAVE;
  }
  else
  {
    fileIO = io::file_descriptor(filename, BOOST_IOS::in | BOOST_IOS::binary);
    mode = MODE_OPEN;
  }

  maxBytes = 134217728;
  maxChunk = 16777216;

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
    
  delete bufferSemaphore;
}

const SasAnalysis&
SasAnalysis::operator <<(SasAtom* sasAtoms)
{
  SasAtom* tmpFrame = new SasAtom[nAtoms];
  if(changeable)
  {
    changeable = false;
    bufferSemaphore = new ip::interprocess_semaphore(bufferSemaphoreMax);

    outFilter.push(io::zlib_compressor());
    outFilter.push(fileIO);
    outArchive = new archive::binary_oarchive(outFilter);

    operationThread = new OperationThread(*this);
  }

  std::copy(sasAtoms, sasAtoms + nAtoms, tmpFrame);
  
  frames.push_back(tmpFrame);
  if(frames.size() == maxFrames)
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
    bufferSemaphore = new ip::interprocess_semaphore(0);

    inFilter.push(io::zlib_decompressor());
    inFilter.push(fileIO);
    inArchive = new archive::binary_iarchive(inFilter);

    operationThread = new OperationThread(*this);
  }

  if(frames.size() == 0)
  {
    bufferSemaphore->wait();
    bufferMutex.lock();
    frames = chunks.front();
    i = frames.begin();

    bufferSemaphoreCount--;
    bufferMutex.unlock();
  }
  else if(i == frames.end())
  {
    bufferMutex.lock();

    for(i = frames.begin(); i < frames.end(); i++)
      delete[] *i;
    frames.clear();
    chunks.pop_front();

    bufferMutex.unlock();
    bufferSemaphore->wait();
    bufferMutex.lock();
    bufferSemaphoreCount--;

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
    bufferSemaphore = new ip::interprocess_semaphore(bufferSemaphoreMax);

    outFilter.push(io::zlib_compressor());
    outFilter.push(fileIO);
    outArchive = new archive::binary_oarchive(outFilter);

    operationThread = new OperationThread(*this);
  }

  bufferSemaphore->wait();
  bufferMutex.lock();
  bufferSemaphoreCount--;
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

  if (this->filename == filename)
    return true;
  
  std::string old_filename = this->filename;
  fs::copy_file(old_filename, filename);
  if(not fs::exists(filename))
    return false;

  this->filename = filename;
  fileIO = io::file_descriptor(filename, BOOST_IOS::app | BOOST_IOS::out |
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

  chunks.push_back(loadChunk(*inArchive));

  return true;
}

void
SasAnalysis::dumpChunk(const std::vector<SasAtom*>& chunk,
                       archive::binary_oarchive& out) const
{
  unsigned int size = chunk.size();
  out << size;

  for
  (
    std::vector<SasAtom*>::const_iterator i = chunk.begin();
    i < chunk.end();
    i++
  )
  {
    const SasAtom* end = *i + nAtoms;
    for
    (
      SasAtom* j = *i;
      j < end;
      j++
    )
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
      in >> *atom;

    chunk.push_back(atoms);
  }

  return chunk;
}

SasAnalysis::OperationThread::OperationThread(SasAnalysis& parent)
{
  this->parent = &parent;
  isStopped = false;
  if(parent.mode == SasAnalysis::MODE_SAVE)
    thread = boost::thread(boost::bind(
                &SasAnalysis::OperationThread::threadSave, boost::ref(*this)));
  else
    thread = boost::thread(boost::bind(
                &SasAnalysis::OperationThread::threadOpen, boost::ref(*this)));
}

SasAnalysis::OperationThread::~OperationThread()
{
  stop();
}

void
SasAnalysis::OperationThread::threadSave()
{
  bool cond;
  
  while(not isStopped)
  {
    parent->bufferMutex.lock();
    cond = (parent->bufferSemaphoreCount == parent->bufferSemaphoreMax - 1 and
       not isStopped);
    
    if(cond)
    {
      parent->bufferMutex.unlock();
      ip::scoped_lock<ip::interprocess_mutex> lock(wakeMutex);      
      wakeCondition.wait(lock);
      parent->bufferMutex.lock();
    }
    else if(isStopped)
    {
      parent->bufferMutex.unlock();
      break;
    }
    
    if(parent->save())
    {
      vector<SasAtom*>& curChunk = parent->chunks.front();
      for
      (
        std::vector<SasAtom*>::iterator i = curChunk.begin();
        i < curChunk.end();
        i++
      )
        delete[] *i;
      curChunk.clear();

      parent->chunks.pop_front();
      parent->bufferSemaphore->post();
      parent->bufferSemaphoreCount++;
    }
    
    parent->bufferMutex.unlock();
  }
  
  parent->bufferMutex.lock();
  while(parent->bufferSemaphoreCount < parent->bufferSemaphoreMax)
  {
    parent->save();
    parent->chunks.pop_front();
    parent->bufferSemaphore->post();
    parent->bufferSemaphoreCount++;    
  }
  parent->bufferMutex.unlock();
}

void
SasAnalysis::OperationThread::threadOpen()
{
  bool cond;

  while(not isStopped)
  {
    parent->bufferMutex.lock();
    cond = (parent->bufferSemaphoreCount == parent->bufferSemaphoreMax - 1 and
            not isStopped);

    if(cond)
    {
      parent->bufferMutex.unlock();
      ip::scoped_lock<ip::interprocess_mutex> lock(wakeMutex);
      wakeCondition.wait(lock);
      parent->bufferMutex.lock();
    }
    else if(isStopped)
    {
      while(parent->bufferSemaphoreCount < parent->bufferSemaphoreMax)
      {
        parent->bufferSemaphore->post();
        parent->bufferSemaphoreCount++;
      }
      parent->bufferMutex.unlock();
      break;
    }

    if(not parent->open())
      isStopped = true;

    parent->bufferSemaphore->post();
    parent->bufferSemaphoreCount++;
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
    return;
  
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
  if(not changeable)
    return;

  unsigned int numFrames =  maxBytes / sizeof(SasAtom) / nAtoms;

  bufferSemaphoreMax = numFrames * nAtoms * sizeof(SasAtom) / maxChunk;
  maxFrames = numFrames / bufferSemaphoreMax;

  if(mode == MODE_SAVE)
    bufferSemaphoreCount = bufferSemaphoreMax - 1;
  else
    bufferSemaphoreCount = 0;

  chunks = boost::circular_buffer<std::vector<SasAtom*> >(bufferSemaphoreMax);
}

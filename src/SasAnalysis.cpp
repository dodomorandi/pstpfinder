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
  nAtoms = gromacs.getAtomsCount();
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
    fileIO = io::file_descriptor(filename, BOOST_IOS::trunc | BOOST_IOS::out);
    mode = MODE_SAVE;
  }
  else
  {
    fileIO = io::file_descriptor(filename, BOOST_IOS::in);
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
    saveThread->stop();
    delete saveThread;
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

    if(mode == MODE_SAVE)
      saveThread = new SaveThread(*this);
  }

  std::copy(sasAtoms, sasAtoms + nAtoms, tmpFrame);
  
  frames.push_back(tmpFrame);
  if(frames.size() == maxFrames)
    flush();
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

    if(mode == MODE_SAVE)
      saveThread = new SaveThread(*this);
  }

  bufferSemaphore->wait();
  bufferMutex.lock();
  bufferSemaphoreCount--;
  chunks.push_back(frames);
  if(curChunk == chunks.end())
    curChunk--;
  bufferMutex.unlock();
  
  saveThread->wakeUp();
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
  fileIO = io::file_descriptor(filename, BOOST_IOS::app | BOOST_IOS::out);
  
  if(not save())
  {
    fileIO.close();
    boost::filesystem::remove(filename);
    this->filename = old_filename;
    fileIO = io::file_descriptor(old_filename, BOOST_IOS::app | BOOST_IOS::out);
    return false;
  }
  
  return true;
}

bool
SasAnalysis::save()
{
  if(mode != MODE_SAVE)
  return false;

  io::filtering_ostream outFilter;
  
  outFilter.strict_sync();
  outFilter.push(io::zlib_compressor());
  outFilter.push(fileIO);
  
  // Let's fill header.
  // We need information about analysis selected options.
  // This means trajectory file, topology file and everything related
  // 
  
  archive::binary_oarchive out(outFilter);
  
  if(chunks.empty())
    return false;
  
  if(curChunk == chunks.end())
    return false;
  unsigned int size = curChunk->size();
  out << size;
  dumpChunk(*curChunk, out);
  
  return true;
}

void
SasAnalysis::dumpChunk(const vector<SasAtom*>& chunk, 
                       archive::binary_oarchive& out) const
{
  for
  (
    vector<SasAtom*>::const_iterator i = chunk.begin();
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
    {
      SasAtomSerializable atom(*j);
      out << atom;
    }
  }
}

SasAnalysis::SaveThread::SaveThread(SasAnalysis& parent)
{
  this->parent = &parent;
  isStopped = false;
  thread = boost::thread(boost::ref(*this));
}

SasAnalysis::SaveThread::~SaveThread()
{
  stop();
}

void
SasAnalysis::SaveThread::operator ()()
{
  bool cond;
  
  while(not isStopped)
  {
    parent->bufferMutex.lock();
    cond = (parent->bufferSemaphoreCount == parent->bufferSemaphoreMax and
       not isStopped);
    parent->bufferMutex.unlock();
    
    if(cond)
    {
      ip::scoped_lock<ip::interprocess_mutex> lock(wakeMutex);      
      wakeCondition.wait(lock);
    }
    else if(isStopped)
      break;
    
    parent->bufferMutex.lock();
    
    if(parent->save())
    {
      for
      (
        std::vector<SasAtom*>::iterator i = parent->curChunk->begin();
        i < parent->curChunk->end();
        i++
      )
        delete[] *i;
      parent->curChunk->clear();

      parent->curChunk++;
      parent->bufferSemaphore->post();
      parent->bufferSemaphoreCount++;
    }
    
    parent->bufferMutex.unlock();
  }
  
  parent->bufferMutex.lock();
  while(parent->bufferSemaphoreCount < parent->bufferSemaphoreMax)
  {
    parent->save();
    parent->curChunk++;
    parent->bufferSemaphore->post();
    parent->bufferSemaphoreCount++;    
  }
  parent->bufferMutex.unlock();
}

void
SasAnalysis::SaveThread::wakeUp()
{
  wakeCondition.notify_one();
}

void
SasAnalysis::SaveThread::stop()
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

  maxFrames = maxBytes / sizeof(SasAtom) / nAtoms;

  bufferSemaphoreMax = maxFrames * nAtoms * sizeof(SasAtom) / maxChunk;
  bufferSemaphoreCount = bufferSemaphoreMax;
  chunks = boost::circular_buffer<std::vector<SasAtom*> >(bufferSemaphoreMax);
  curChunk = chunks.begin();
}

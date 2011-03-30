#include "SasAnalysis.h"

using namespace Gromacs;
namespace archive = boost::archive;
namespace io = boost::iostreams;
namespace ip = boost::interprocess;

SasAnalysis::SasAnalysis(unsigned int nAtoms,
                         unsigned int maxBytes,
                         unsigned int maxChunk,
                         std::string saveFile)
{
  this->nAtoms = nAtoms;
  gromacs = 0;
  init(maxBytes, maxChunk, saveFile);
}

SasAnalysis::SasAnalysis(const Gromacs& gromacs,
                         unsigned int maxBytes,
                         unsigned int maxChunk,
                         std::string saveFile)
{
  nAtoms = gromacs.getAtomsCount();
  this->gromacs = &gromacs;
  init(maxBytes, maxChunk, saveFile);
}

void
SasAnalysis::init(unsigned int maxBytes, unsigned int maxChunk, string saveFile)
{
  fileOut = io::file_descriptor_sink("/tmp/trapof.csf",
                                     ios_base::trunc | ios_base::out);
  maxFrames = maxBytes / sizeof(SasAtom) / nAtoms;
  
  bufferSemaphoreMax = maxFrames * nAtoms * sizeof(SasAtom) / maxChunk;
  bufferSemaphoreCount = bufferSemaphoreMax;
  chunks = boost::circular_buffer<std::vector<SasAtom*> >(bufferSemaphoreMax);
  curChunk = chunks.begin();
  bufferSemaphore = new ip::interprocess_semaphore(bufferSemaphoreMax);
  saveThread = new SaveThread(*this);
}

SasAnalysis::~SasAnalysis()
{
  if(frames.size() != 0)
    flush();
  
  saveThread->stop();
  delete saveThread;
    
  delete bufferSemaphore;
}

const SasAnalysis&
SasAnalysis::operator <<(SasAtom* sasAtoms)
{
  SasAtom* tmpFrame = new SasAtom[nAtoms];
  std::copy(sasAtoms, sasAtoms + nAtoms, tmpFrame);
  
  frames.push_back(tmpFrame);
  if(frames.size() == maxFrames)
    flush();
  return *this;
}

void
SasAnalysis::flush()
{
  bufferSemaphore->wait();
  bufferMutex.lock();
  bufferSemaphoreCount--;
  chunks.push_back(frames);
  bufferMutex.unlock();
  
  saveThread->wakeUp();
}

bool
SasAnalysis::save(const std::string& filename)
{
  if (this->filename == filename)
    return true;
  
  std::string old_filename = this->filename;
  this->filename = filename;
  fileOut = io::file_descriptor_sink(filename, ios_base::trunc | ios_base::out);
  
  if(not save())
  {
    fileOut.close();
    boost::filesystem::remove(filename);
    this->filename = old_filename;
    fileOut = io::file_descriptor_sink(old_filename,
                                       std::ios_base::app | std::ios_base::out);
    return false;
  }
  
  return true;
}

bool
SasAnalysis::save()
{
  io::filtering_ostream outFilter;
  
  outFilter.strict_sync();
  outFilter.push(io::zlib_compressor());
  outFilter.push(fileOut);
  
  // Let's fill header.
  // We need information about analysis selected options.
  // This means trajectory file, topology file and everything related
  // 
  
  archive::binary_oarchive out(outFilter);
  
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
  while(not isStopped)
  {
    ip::scoped_lock<ip::interprocess_mutex> lock(parent->bufferMutex);
    parent->bufferMutex.lock();
    if(parent->bufferSemaphoreCount == parent->bufferSemaphoreMax and
       not isStopped)
      wakeCondition.wait(lock);
      
    else if(isStopped)
    {
      parent->bufferMutex.unlock();
      break;
    }
    
    parent->save();

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
    return
  
  parent->bufferMutex.lock();
  isStopped = true;
  wakeCondition.notify_one();
  parent->bufferMutex.unlock();
  
  thread.join();
}
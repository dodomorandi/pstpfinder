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
    operationThread->stop();
    delete operationThread;
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

    operationThread = new OperationThread(*this);
  }

  if(frames.size() == 0)
  {
    bufferSemaphore->wait();
    bufferMutex.lock();
    frames = *(curChunk++);
    i = frames.begin();

    bufferSemaphoreCount--;
    bufferMutex.unlock();
  }
  else if(i == frames.end())
  {
    // TODO: Control for EOF
    bufferMutex.lock();

    for(i = frames.begin(); i < frames.end(); i++)
      delete[] *i;

    if(curChunk == chunks.end())
    {
      bufferMutex.unlock();
      bufferSemaphore->wait();
      bufferMutex.lock();
      bufferSemaphoreCount--;
      if(curChunk == chunks.end())
      {
        bufferMutex.unlock();
        sasAtom = 0;
        return *this;
      }
    }

    frames = *(curChunk++);
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

    operationThread = new OperationThread(*this);
  }

  bufferSemaphore->wait();
  bufferMutex.lock();
  bufferSemaphoreCount--;
  chunks.push_back(frames);
  if(curChunk == chunks.end())
    curChunk--;
  bufferMutex.unlock();
  
  operationThread->wakeUp();
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

bool
SasAnalysis::open()
{
  if(mode != MODE_OPEN)
    return false;

  if(not fileIO.is_open())
    return false;

  io::filtering_istream inFilter;

  inFilter.strict_sync();
  inFilter.push(io::zlib_decompressor());
  inFilter.push(fileIO);

  archive::binary_iarchive in(inFilter);
  chunks.push_back(loadChunk(in));

  return true;
}

void
SasAnalysis::dumpChunk(const std::vector<SasAtom*>& chunk,
                       archive::binary_oarchive& out) const
{
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
    cond = (parent->bufferSemaphoreCount == parent->bufferSemaphoreMax and
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
SasAnalysis::OperationThread::threadOpen()
{
  bool cond;

  while(not isStopped)
  {
    parent->bufferMutex.lock();
    cond = (parent->bufferSemaphoreCount == parent->bufferSemaphoreMax and
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

    if(not parent->open())
      isStopped = true;
    else if(parent->curChunk == parent->chunks.end())
      parent->curChunk--;

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

  maxFrames = maxBytes / sizeof(SasAtom) / nAtoms;

  bufferSemaphoreMax = maxFrames * nAtoms * sizeof(SasAtom) / maxChunk;
  if(mode == MODE_SAVE)
    bufferSemaphoreCount = bufferSemaphoreMax;
  else
    bufferSemaphoreCount = 0;

  chunks = boost::circular_buffer<std::vector<SasAtom*> >(bufferSemaphoreMax);
  curChunk = chunks.begin();
}

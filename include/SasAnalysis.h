#ifndef _SASANALYSIS_H
#define _SASANALYSIS_H

#include "SasAtom.h"
#include "Gromacs.h"

#include <boost/thread.hpp>
#include <boost/circular_buffer.hpp>
#include <boost/filesystem.hpp>
#include <boost/thread/detail/thread.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/interprocess/sync/interprocess_semaphore.hpp>
#include <boost/interprocess/sync/interprocess_mutex.hpp>
#include <boost/interprocess/sync/interprocess_condition.hpp>
#include <boost/interprocess/sync/scoped_lock.hpp>

#include <vector>
#include <string>
#include <iostream>
#include <sstream>

namespace Gromacs
{
  class SasAnalysis
  {
  public:
    SasAnalysis(unsigned int nAtoms,
                std::string filename = "/tmp/sas.csf",
                bool savingMode = true);
                
    SasAnalysis(const Gromacs& gromacs,
                std::string filename = "/tmp/sas.csf",
                bool savingMode = true);
    
    ~SasAnalysis();
    const SasAnalysis& operator <<(SasAtom* sasAtoms);
    bool setMaxBytes(unsigned long bytes);
    unsigned long getMaxBytes();
    bool setMaxChunkSize(unsigned long bytes);
    unsigned long getMaxChunkSize();

  private:
    boost::circular_buffer<std::vector<SasAtom*> > chunks;
    std::vector<SasAtom*> frames;
    unsigned int nAtoms;
    std::string filename;
    boost::iostreams::file_descriptor fileIO;
    const Gromacs* gromacs;
    unsigned long maxFrames, maxBytes, maxChunk;
    boost::circular_buffer<std::vector<SasAtom* > >::iterator curChunk;
    mutable boost::interprocess::interprocess_semaphore* bufferSemaphore;
    unsigned int bufferSemaphoreCount;
    unsigned int bufferSemaphoreMax;
    mutable boost::interprocess::interprocess_mutex bufferMutex;
    bool changeable;

    enum
    {
      MODE_OPEN,
      MODE_SAVE
    } mode;

    class SaveThread
    {
    public:
      SaveThread(SasAnalysis& parent);
      ~SaveThread();
      void wakeUp();
      void stop();
      void operator ()();
    private:
      SasAnalysis* parent;
      bool isStopped;
      boost::thread thread;
      boost::interprocess::interprocess_condition wakeCondition;
      boost::interprocess::interprocess_mutex wakeMutex;
    };

    SaveThread* saveThread;

    void init(string saveFile, bool savingMode = true);
    void dumpChunk(const std::vector<SasAtom*>& chunk,
                   boost::archive::binary_oarchive& out) const;
    void flush();
    bool save(const std::string& filename);
    bool save();

    void updateChunks();
  };
};

#endif


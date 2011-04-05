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
                unsigned int maxBytes = 134217728,
                unsigned int maxChunk = 16777216,
                std::string saveFile = "/tmp/trapof.csf");
                
    SasAnalysis(const Gromacs& gromacs,
                unsigned int maxBytes = 134217728,
                unsigned int maxChunk = 16777216,
                std::string saveFile = "/tmp/trapof.csf");
    
    ~SasAnalysis();
    const SasAnalysis& operator <<(SasAtom* sasAtoms);
  private:
    boost::circular_buffer<std::vector<SasAtom*> > chunks;
    std::vector<SasAtom*> frames;
    unsigned int nAtoms;
    std::string filename;
    boost::iostreams::file_descriptor_sink fileOut;
    const Gromacs* gromacs;
    unsigned long maxFrames;
    boost::circular_buffer<std::vector<SasAtom* > >::iterator curChunk;
    mutable boost::interprocess::interprocess_semaphore* bufferSemaphore;
    unsigned int bufferSemaphoreCount;
    unsigned int bufferSemaphoreMax;
    mutable boost::interprocess::interprocess_mutex bufferMutex;

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

    void init(unsigned int maxBytes, unsigned int maxChunk, string saveFile);
    void dumpChunk(const std::vector<SasAtom*>& chunk,
                   boost::archive::binary_oarchive& out) const;
    void flush();
    bool save(const std::string& filename);
    bool save();    
  };
};

#endif


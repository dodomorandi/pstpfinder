#ifndef _SASANALYSIS_H
#define _SASANALYSIS_H

#include "SasAtom.h"
#include "Gromacs.h"

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/circular_buffer.hpp>

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
                unsigned int maxChunk = 16777216);
    SasAnalysis(const Gromacs& gromacs,
                unsigned int maxBytes = 134217728,
                unsigned int maxChunk = 16777216);
    ~SasAnalysis();
    const SasAnalysis& operator <<(SasAtom* sasAtoms);
    bool save(const std::string& filename);
    bool save() const;
  private:
    boost::circular_buffer<std::vector<SasAtom*> > chunks;
    std::vector<SasAtom*> frames;
    unsigned int nAtoms;
    std::string filename;
    const Gromacs* gromacs;
    unsigned long maxFrames;
    
    void dumpChunk(const std::vector<SasAtom*>& chunk,
                   boost::archive::binary_oarchive& out) const;
  };
};

#endif


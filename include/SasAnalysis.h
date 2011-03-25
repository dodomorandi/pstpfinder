#ifndef _SASANALYSIS_H
#define _SASANALYSIS_H

#include "SasAtom.h"

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <vector>
#include <string>
#include <iostream>
#include <sstream>

namespace Gromacs
{
  class SasAnalysis
  {
  public:
    SasAnalysis(unsigned int nAtoms);
    ~SasAnalysis();
    const SasAnalysis& operator <<(SasAtom* sasAtoms);
    bool save(const std::string& filename);
    bool save() const;
  private:
    std::vector<SasAtom*> atoms;
    unsigned int nAtoms;
    std::string filename;
  };
};

#endif


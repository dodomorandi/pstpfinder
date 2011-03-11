#ifndef _SASANALYSIS_H
#define _SASANALYSIS_H

#include "SasAtom.h"

#include <vector>

namespace Gromacs
{
  class SasAnalysis
  {
  public:
    SasAnalysis(unsigned int nAtoms);
    ~SasAnalysis();
    const SasAnalysis& operator <<(SasAtom* sasAtoms);
  private:
    std::vector<SasAtom*> atoms;
    unsigned int nAtoms;
  };
};

#endif


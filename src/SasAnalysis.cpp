#include "SasAnalysis.h"

#include <vector>

using namespace Gromacs;

SasAnalysis::SasAnalysis(unsigned int nAtoms)
{
  this->nAtoms = nAtoms;
}

SasAnalysis::~SasAnalysis()
{
  for
  (
    std::vector<SasAtom*>::iterator i = atoms.begin();
    i < atoms.end();
    i++
  )
    delete[] *i;
}

const SasAnalysis&
SasAnalysis::operator <<(SasAtom* sasAtoms)
{
  SasAtom* tmpAtoms = new SasAtom[nAtoms];
  std::copy(sasAtoms, sasAtoms + nAtoms, tmpAtoms);
  
  atoms.push_back(tmpAtoms);
  return *this;
}

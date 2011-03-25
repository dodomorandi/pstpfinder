#include "SasAnalysis.h"

using namespace Gromacs;
namespace archive = boost::archive;
namespace io = boost::iostreams;

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

bool
SasAnalysis::save(const std::string& filename)
{
  std::string old_filename = this->filename;
  this->filename = filename;
  
  if(!save())
  {
    this->filename = old_filename;
    return false;
  }
  
  return true;
}

bool
SasAnalysis::save() const
{
  using namespace std;
  
  io::filtering_ostream outFilter;
  io::file_descriptor_sink fileOut( "/tmp/trapof.csf",
                                    ios_base::trunc | ios_base::out);
  
  outFilter.strict_sync();
  outFilter.push(io::zlib_compressor());
  outFilter.push(fileOut);
  
  archive::binary_oarchive out(outFilter);
  
  for
  (
    vector<SasAtom*>::const_iterator i = atoms.begin();
    i < atoms.end();
    i++
  )
  {
    SasAtom* end = *i + nAtoms;
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
  
  return true;
}

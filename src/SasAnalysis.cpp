#include "SasAnalysis.h"

using namespace Gromacs;
namespace archive = boost::archive;
namespace io = boost::iostreams;

SasAnalysis::SasAnalysis(unsigned int nAtoms,
                         unsigned int maxBytes,
                         unsigned int maxChunk)
{
  this->nAtoms = nAtoms;
  gromacs = 0;
  maxFrames = maxBytes / sizeof(SasAtom);
}

SasAnalysis::SasAnalysis(const Gromacs& gromacs,
                         unsigned int maxBytes,
                         unsigned int maxChunk)
{
  nAtoms = gromacs.getAtomsCount();
  this->gromacs = &gromacs;
  maxFrames = maxBytes / sizeof(SasAtom);
}

SasAnalysis::~SasAnalysis()
{
  for
  (
    std::vector<SasAtom*>::iterator i = frames.begin();
    i < frames.end();
    i++
  )
    delete[] *i;
}

const SasAnalysis&
SasAnalysis::operator <<(SasAtom* sasAtoms)
{
  SasAtom* tmpFrame = new SasAtom[nAtoms];
  std::copy(sasAtoms, sasAtoms + nAtoms, tmpFrame);
  
  frames.push_back(tmpFrame);
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
  
  // Let's fill header.
  // We need information about analysis selected options.
  // This means trajectory file, topology file and everything related
  // 
  
  archive::binary_oarchive out(outFilter);
  dumpChunk(frames, out);
  
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
}
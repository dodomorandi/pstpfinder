#include "Pittpi.h"
#include <cstring>

using namespace Gromacs;
using namespace std;

Group::Group(const Residue& refResidue)
{
  reference = &refResidue;
}

Group&
Group::operator <<(const Residue& reference)
{
  residues.push_back(&reference);

  return *this;
}

const vector<const Residue*>&
Group::getResidues() const
{
  return residues;
}

Pittpi::Pittpi(const Gromacs& gromacs, float radius, unsigned long threshold)
{
  averageStructure = gromacs.getAverageStructure();

  std::vector<Group> groups = makeGroups(radius);
}

std::vector<Group>
Pittpi::makeGroups(float radius)
{
  std::vector<Group> groups;
  vector<Atom> centers;
  const vector<Residue>& residues = averageStructure.residues();

  // Calculate the center for every sidechain (excluding PRO)
  centers.reserve(residues.size());
  for
  (
    vector<Residue>::const_iterator i = residues.begin();
    i < residues.end();
    i++
  )
  {
    Atom center(0);
    const vector<PdbAtom>& atoms = i->atoms;

    if(strcmp(i->getAtomByType("H1").type, "UNK") != 0)
      continue;

    if(i->type == AA_PRO)
    {
      centers.push_back(center);
      continue;
    }
    else if(i->type == AA_GLY)
    {
      centers.push_back(i->getAtomByType("CA"));
      continue;
    }
    unsigned int count = 0;
    for(vector<PdbAtom>::const_iterator j = atoms.begin(); j < atoms.end(); j++)
    {
      if(strcmp(j->type, "N") != 0 and strcmp(j->type, "CA") != 0 and
         strcmp(j->type, "H") != 0 and strcmp(j->type, "C") != 0 and
         strcmp(j->type, "O") != 0 and strcmp(j->type, "HA") != 0)
      {
        center += *j;
        count++;
      }
    }
    center /= (float)count;

    centers.push_back(center);
  }

  vector<Atom>::const_iterator centersBegin = centers.begin();

  for
  (
    vector<Residue>::const_iterator i = residues.begin();
    i < residues.end();
    i++
  )
  {
    Group group(*i);
    const PdbAtom& hAtom = i->getAtomByType("H");
    if(strcmp(hAtom.type, "UNK") == 0)
      continue;

    for
    (
      vector<Atom>::const_iterator j = centersBegin;
      j < centers.end();
      j++
    )
    {
      if(distance(centersBegin, j) ==
         distance(vector<Residue>::const_iterator(residues.begin()), i))
        continue;

      if(hAtom.distance(*j) <= radius)
        group << residues[distance(centersBegin, j)];
    }

    groups.push_back(group);
  }

//  FIXME: NEVER TESTED!!

  return groups;
}

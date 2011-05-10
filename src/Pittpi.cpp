#include "Pittpi.h"
#include <cstring>

using namespace Gromacs;
using namespace std;

Pittpi::Pittpi(const Gromacs& gromacs, float radius, unsigned long threshold)
{
  averageStructure = gromacs.getAverageStructure();

  vector<const Residue*> groups = makeGroups(threshold);
}

vector<const Residue*>
Pittpi::makeGroups(unsigned long threshold)
{
  vector<const Residue*> groups;
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

    if(i->type == AA_PRO)
    {
      centers.push_back(center);
      continue;
    }
    else if(i->type == AA_GLY)
    {
      for
      (
        vector<PdbAtom>::const_iterator j = atoms.begin();
        j < atoms.end();
        j++
      )
        if(strcmp(j->type, "CA"))
        {
          center.x = j->x;
          center.y = j->y;
          center.z = j->z;
        }
      centers.push_back(center);
      continue;
    }
    unsigned int size = atoms.size();
    for(vector<PdbAtom>::const_iterator j = atoms.begin(); j < atoms.end(); j++)
    {
      if(strcmp(j->type, "N") != 0 and strcmp(j->type, "CA") != 0 and
         strcmp(j->type, "H") != 0 and strcmp(j->type, "C") != 0 and
         strcmp(j->type, "C") != 0 and strcmp(j->type, "O") != 0 and
         strcmp(j->type, "HA") != 0)
      {
        center += *j;
        center /= size;
      }
    }

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

      if(hAtom.distance(*j) <= threshold)
        groups.push_back(residues.data() + distance(centersBegin, j));
    }
  }

//  FIXME: NEVER TESTED!!

  return groups;
}

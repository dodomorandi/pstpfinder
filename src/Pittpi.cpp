#include "Pittpi.h"
#include "SasAtom.h"
#include "SasAnalysis.h"
#include <cstring>

using namespace Gromacs;
using namespace std;

Group::Group(const Residue& refResidue)
{
  reference = &refResidue.getAtomByType("H");
  zeros = 0;
}

Group::Group(const PdbAtom& refAtomH)
{
  reference = &refAtomH;
  zeros = 0;
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

const PdbAtom&
Group::getCentralH() const
{
  return *reference;
}

Pittpi::Pittpi(const Gromacs& gromacs,
               const std::string& sasAnalysisFileName,
               float radius,
               unsigned long threshold)
{
  m_gromacs = &gromacs;

  averageStructure = gromacs.getAverageStructure();
  vector<Group> groups = makeGroups(radius);

  fillGroups(groups, sasAnalysisFileName);
}

vector<Group>
Pittpi::makeGroups(float radius)
{
  vector<Group> groups;
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
    const PdbAtom& hAtom = i->getAtomByType("H");
    if(strcmp(hAtom.type, "UNK") == 0)
      continue;

    Group group(*i);
    for
    (
      vector<Atom>::const_iterator j = centersBegin;
      j < centers.end();
      j++
    )
    {
      if(residues[distance(centersBegin, j)].type == AA_PRO)
        continue;

      if(distance(centersBegin, j) ==
         distance(vector<Residue>::const_iterator(residues.begin()), i))
        continue;

      if(hAtom.distance(*j) <= radius)
        group << residues[distance(centersBegin, j)];
    }

    groups.push_back(group);
  }

//  FIXME: Missing sadic alghoritm
//  FIXME: NEVER TESTED!!

  return groups;
}

void
Pittpi::fillGroups(vector<Group>& groups, const string& sasAnalysisFileName)
{
  float* sas;
  float* meanSas;
  float* fIndex;
  SasAtom* sasAtoms = 0;
  SasAnalysis* sasAnalysis;

  const Gromacs& gromacs = *m_gromacs;

  const float frames = gromacs.getFramesCount();
  vector<int> protein = gromacs.getGroup("Protein");

  meanSas = new float[protein.size()]();

  /* First of all we need to calculate SAS means */
  sasAnalysis = new SasAnalysis(gromacs, sasAnalysisFileName, false);
  (*sasAnalysis) >> sasAtoms;
  while(sasAtoms != 0)
  {
    fIndex = meanSas;
    for
    (
      SasAtom* atom = sasAtoms;
      atom < sasAtoms + protein.size();
      atom++, fIndex++
    )
      *fIndex += atom->sas;

    (*sasAnalysis) >> sasAtoms;
  }

  {
    for(fIndex = meanSas; fIndex < meanSas + protein.size(); fIndex++)
      *fIndex /= frames;
  }

  delete sasAnalysis;

  /* Let's prepare groups sas vectors */
  for
  (
    vector<Group>::iterator i = groups.begin();
    i < groups.end();
    i++
  )
    i->sas.reserve(frames);

  /* Now we have to normalize values and store results per group */
  unsigned int curFrame = 0;
  sas = new float[protein.size()];
  sasAnalysis = new SasAnalysis(gromacs, sasAnalysisFileName, false);
  (*sasAnalysis) >> sasAtoms;
  while(sasAtoms != 0)
  {
    fIndex = sas;
    for
    (
      SasAtom* atom = sasAtoms;
      atom < sasAtoms + protein.size();
      atom++, fIndex++
    )
      *fIndex = atom->sas;

    for
    (
      vector<Group>::iterator i = groups.begin();
      i < groups.end();
      i++
    )
    {
      if(sas[i->getCentralH().index] == 0)
      {
        i->sas[curFrame] = 0;
        i->zeros++;
        continue;
      }

      const vector<const Residue*>& residues = i->getResidues();
      for
      (
        vector<const Residue*>::const_iterator j = residues.begin();
        j < residues.end();
        j++
      )
        for
        (
          vector<PdbAtom>::const_iterator k = (*j)->atoms.begin();
          k < (*j)->atoms.end();
          k++
        )
          if(k->type[0] == 'H')
            i->sas[curFrame] += sas[k->index - 1] / meanSas[k->index -1];

      i->sas[curFrame] /= frames;
      if(i->sas[curFrame] == 0)
        i->zeros++;
    }

    curFrame++;
    (*sasAnalysis) >> sasAtoms;
  }
  delete sasAnalysis;
}

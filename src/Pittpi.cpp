#include "Pittpi.h"
#include "SasAtom.h"
#include "SasAnalysis.h"
#include <cstring>
#include <algorithm>

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

bool
Group::sortByZeros(const Group& a, const Group& b)
{
  return (a.zeros > b.zeros);
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

  /*
   * This part is a "legacy" method. It have been implemented in perl time ago
   * and needs refactory. The main problem is math related, because we have to
   * find a good solution to take "consecutively opened pocked" above a certain
   * threshold. With every frame (and every SAS value) it could be not so easy
   * to develop a GOOD alghorithm. For now we implement only the old method used
   * until now.
   * -- Edoardo Morandi
   */

#define PS_PER_SAS 5
  unsigned int const frameStep = float(PS_PER_SAS) / gromacs.getTimeStep();
  unsigned int const frames = gromacs.getFramesCount();
  unsigned int const newSasCount = frames / frameStep +
                             (frames % frameStep == 0) ? 0 : 1;
  unsigned int const dist2frame = PS_PER_SAS * gromacs.getTimeStep();
  vector<Group> meanGroups = groups;
  for
  (
    vector<Group>::iterator i = meanGroups.begin(), j = groups.begin();
    i < meanGroups.end();
    i++, j++
  )
  {
    i->sas.reserve(newSasCount);
    for(vector<float>::iterator k = j->sas.begin(); k < j->sas.end(); k++)
    {
      float mean = 0;
      vector<float>::iterator end = k + frameStep;
      for(;k < end; k++)
        mean += *k;

      mean /= frameStep;
      i->sas.push_back(mean);
    }
  }

  sort(meanGroups.begin(), meanGroups.end(), Group::sortByZeros);

  unsigned int noZeroPass = threshold / PS_PER_SAS;
  if(noZeroPass < 20)
    noZeroPass = 0;
  else if(noZeroPass < 40)
    noZeroPass = 1;
  else if(noZeroPass < 60)
    noZeroPass = 2;
  else
    noZeroPass = 3;

  vector<Pocket> pockets;
  for
  (
    vector<Group>::iterator i = meanGroups.begin();
    i < meanGroups.end();
    i++
  )
  {
    vector<float>::iterator startPocket = i->sas.end();
    unsigned int notOpenCounter = 0;
    float* maxFrame;
    float mean = 0;
    unsigned int percOpenedNotZero = 1;

    for
    (
      vector<float>::iterator j = i->sas.begin();
      j < i->sas.end();
      j++
    )
    {
      if(startPocket == i->sas.end() and *j > 1)
      {
        startPocket = j;
        maxFrame = &(*j);
        mean = *j;
        percOpenedNotZero = 1;
      }
      else if(startPocket != i->sas.end() and *j < 1)
      {
        if(notOpenCounter < noZeroPass)
        {
          notOpenCounter++;
          mean += *j;
        }
        else
        {
          if(static_cast<unsigned int>(distance(startPocket, j) * dist2frame) >=
             threshold)
          {
            Pocket pocket;
            pocket.group = &(*i);
            pocket.startFrame = distance(i->sas.begin(), startPocket) *
                                dist2frame;
            pocket.endFrame = distance(i->sas.begin(), j) * dist2frame;
            pocket.maxAreaFrame = static_cast<int>
                                  (maxFrame - &(*i->sas.begin())) * dist2frame;
            pocket.openingFraction = static_cast<float>(percOpenedNotZero) /
                                     distance(startPocket, j);

            mean /=  distance(startPocket, j);
            vector<float>::iterator nearToMean = startPocket;
            for(vector<float>::iterator k = startPocket + 1; k < j; k++)
              if(abs(mean - *nearToMean) > abs(mean - *k))
                nearToMean = k;
            pocket.meanNearFrame = distance(i->sas.begin(), nearToMean) *
                                   dist2frame;

            pockets.push_back(pocket);
          }
          startPocket = i->sas.end();
          notOpenCounter = 0;
        }
      }
      else
      {
        if(*j > *maxFrame)
          maxFrame = &(*j);
        mean += *j;
        percOpenedNotZero++;
      }
    }
  }
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
  vector<float>::iterator curFrame;
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
      i->sas.push_back(0);
      curFrame = i->sas.end() - 1;

      if(sas[i->getCentralH().index] == 0)
      {
        i->zeros++;
        continue;
      }

      const vector<const Residue*>& residues = i->getResidues();
      unsigned int counter = 0;
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
          {
            (*curFrame) += sas[k->index - 1] / meanSas[k->index - 1];
            counter++;
          }

      (*curFrame) /= counter;
      if(*curFrame == 0)
        i->zeros++;
    }

    (*sasAnalysis) >> sasAtoms;
  }
  delete sasAnalysis;
  delete[] meanSas;
  delete[] sas;
}

/*
 *  This file is part of PSTP-finder, an user friendly tool to analyze GROMACS
 *  molecular dynamics and find transient pockets on the surface of proteins.
 *  Copyright (C) 2011 Edoardo Morandi.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "Pittpi.h"
#include "SasAtom.h"
#include "SasAnalysis.h"
#include <cstring>

#include <boost/interprocess/sync/scoped_lock.hpp>

using namespace Gromacs;
using namespace std;

Group::Group(const Residue& refResidue)
{
  referenceRes = &refResidue;
  referenceAtom = &refResidue.getAtomByType("H");
  zeros = 0;
}

Group::Group(const PdbAtom& refAtomH)
{
  referenceRes = 0;
  referenceAtom = &refAtomH;
  zeros = 0;
}

Group::Group(const PdbAtom& refAtomH, const Protein& protein)
{
  referenceRes = &protein.getResidueByAtom(refAtomH);
  referenceAtom = &refAtomH;
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
  return *referenceAtom;
}

const Residue&
Group::getCentralRes() const
{
  if(referenceRes != 0)
    return *referenceRes;
  else
  {
    static Residue nullRes;
    nullRes.type = AA_UNK;
    return nullRes;
  }
}

bool
Group::sortByZeros(const Group& a, const Group& b)
{
  return (a.zeros > b.zeros);
}

void
Group::switchReference(const Protein& structure)
{
  vector<const Residue*> oldRes = vector<const Residue*>(residues);
  residues.clear();
  for
  (
    vector<const Residue*>::const_iterator i = oldRes.begin();
    i < oldRes.end();
    i++
  )
    residues.push_back(&structure.getResidueByIndex((*i)->index));
  referenceRes = &structure.getResidueByIndex(referenceRes->index);
  referenceAtom = &referenceRes->getAtomByType("H");
}

Pittpi::Pittpi(const Gromacs& gromacs,
               const std::string& sasAnalysisFileName,
               float radius,
               unsigned long threshold) :
  gromacs(gromacs)
{
  this->sasAnalysisFileName = sasAnalysisFileName;
  this->radius = radius;
  this->threshold = threshold;
  sync = true;
  __status = 0;

  pittpiThread = boost::thread(&Pittpi::pittpiRun, boost::ref(*this));
}

Pittpi::Pittpi(const Pittpi& pittpi) :
  gromacs(pittpi.gromacs)
{
  clone(pittpi);
}

Pittpi::Pittpi(const Pittpi& pittpi, const Gromacs& gromacs) :
  gromacs(gromacs)
{
  clone(pittpi);
}

void
Pittpi::clone(const Pittpi& pittpi)
{
  sasAnalysisFileName = pittpi.sasAnalysisFileName;
  radius = pittpi.radius;
  threshold = pittpi.threshold;
  averageStructure = gromacs.getAverageStructure();
  averageStructure.forceUnlock();
  sync = pittpi.sync;
  __status = pittpi.__status;
  groups = vector<Group>(pittpi.groups);
  meanGroups = vector<Group>(pittpi.meanGroups);
  pockets = vector<Pocket>(pittpi.pockets);

  for(vector<Group>::iterator i = groups.begin(); i < groups.end(); i++)
    i->switchReference(averageStructure);

  for(vector<Group>::iterator i = meanGroups.begin(); i < meanGroups.end(); i++)
    i->switchReference(averageStructure);

  vector<Pocket>::iterator i;
  vector<Pocket>::const_iterator j;
  for
  (
    i = pockets.begin(), j = pittpi.pockets.begin();
    i < pockets.end();
    i++, j++
  )
    i->group = &meanGroups[j->group - pittpi.meanGroups.data()];
}

Pittpi::~Pittpi()
{
  join();
}

void
Pittpi::join()
{
  pittpiThread.join();
}

bool
Pittpi::isFinished()
{
  if(not sync)
    pittpiThread.join();
  return pittpiThread == boost::thread();
}

void
Pittpi::setStatus(float value)
{
  statusMutex.lock();
  __status = value;
  statusMutex.unlock();
  nextStatusCondition.notify_all();
}

float
Pittpi::getStatus() const
{
  float retval;
  statusMutex.lock();
  retval = __status;
  statusMutex.unlock();
  return retval;
}

void
Pittpi::waitNextStatus()
{
  if(isFinished())
    return;

  boost::interprocess::scoped_lock<boost::interprocess::interprocess_mutex>
                                                     lock(nextStatusMutex);
  nextStatusCondition.wait(lock);
}

void
Pittpi::pittpiRun()
{
  averageStructure = gromacs.getAverageStructure();
  groups = makeGroups(radius);

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
                             ((frames % frameStep == 0) ? 0 : 1);
  meanGroups = vector<Group>(groups);

  setStatus(0);
  for
  (
    vector<Group>::iterator i = meanGroups.begin(), j = groups.begin();
    j < groups.end();
    i++, j++
  )
  {
    i->sas.clear();
    i->zeros /= frameStep;

    i->sas.reserve(newSasCount);
    for(vector<float>::iterator k = j->sas.begin(); k < j->sas.end(); k++)
    {
      float mean = 0;

      vector<float>::iterator end = k + frameStep;
      if(end > j->sas.end())
        end = j->sas.end();

      for(; k < end; k++)
        mean += *k;

      mean /= frameStep;
      i->sas.push_back(mean);
    }

    setStatus(static_cast<float>(distance(groups.begin(), j) + 1)
              / groups.size());
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

  setStatus(0);
  for
  (
    vector<Group>::iterator i = meanGroups.begin();
    i < meanGroups.end();
    i++
  )
  {
    vector<float>::iterator startPocket = i->sas.end();
    unsigned int notOpenCounter = 0;
    float* maxFrame = 0;
    float mean = 0;

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
          if(static_cast<unsigned int>(distance(startPocket, j) *
            (float)PS_PER_SAS) >= threshold)
          {
            Pocket pocket;
            pocket.group = &(*i);
            pocket.startFrame = distance(i->sas.begin(), startPocket) *
                                frameStep + 1;
            pocket.startPs = distance(i->sas.begin(), startPocket) *
                                PS_PER_SAS;
            pocket.endFrame = distance(i->sas.begin(), j) * frameStep + 1;
            pocket.endPs = distance(i->sas.begin(), j) * PS_PER_SAS;
            pocket.width = pocket.endPs - pocket.startPs;
            pocket.maxAreaFrame = static_cast<int>
                                  (maxFrame - &(*i->sas.begin())) * frameStep
                                  + 1;
            pocket.maxAreaPs = static_cast<float>
                               (maxFrame - &(*i->sas.begin())) * PS_PER_SAS;
            pocket.openingFraction = static_cast<float>
                                     (distance(startPocket, j)) /
                                     (newSasCount - i->zeros);

            mean /=  distance(startPocket, j);
            vector<float>::iterator nearToAverage = startPocket;
            for(vector<float>::iterator k = startPocket + 1; k < j; k++)
              if(abs(mean - *nearToAverage) > abs(mean - *k))
                nearToAverage = k;
            pocket.averageNearFrame = distance(i->sas.begin(), nearToAverage) *
                          frameStep + 1;
            pocket.averageNearPs = distance(i->sas.begin(), nearToAverage) *
                                      PS_PER_SAS;

            pockets.push_back(pocket);
          }
          startPocket = i->sas.end();
          notOpenCounter = 0;
        }
      }
      else
      {
        if(maxFrame == 0 or *j > *maxFrame)
          maxFrame = &(*j);
        mean += *j;
      }
    }

    setStatus(static_cast<float>(distance(meanGroups.begin(), i) + 1)
              / meanGroups.size());
  }

  sync = false;
  setStatus(1);
}

vector<Group>
Pittpi::makeGroups(float radius)
{
  vector<Group> groups;
  vector<Atom> centers;
  const vector<Residue>& residues = averageStructure.residues();
  radius /= 10.0;

  // Calculate the center for every sidechain (excluding PRO)
  setStatus(0);
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
    setStatus(static_cast<float>(distance(residues.begin(), i) + 1) /
              residues.size());
  }

  vector<Atom>::const_iterator centersBegin = centers.begin();

  setStatus(0);
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

      if(hAtom.distance(*j) <= radius)
        group << residues[distance(centersBegin, j)];
    }

    groups.push_back(group);
    setStatus(static_cast<float>(distance(residues.begin(), i) + 1) /
              residues.size());
  }

//  FIXME: Missing sadic algorithm

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
  unsigned int counter = 0;

  const float frames = gromacs.getFramesCount();
  vector<int> protein = gromacs.getGroup("Protein");
  const int nAtoms = protein.size();

  meanSas = new float[nAtoms]();
  unsigned int* notZero = new unsigned int[nAtoms]();

  /* First of all we need to calculate SAS means */
  setStatus(0);
  sasAnalysis = new SasAnalysis(gromacs, sasAnalysisFileName, false);
  (*sasAnalysis) >> sasAtoms;
  while(sasAtoms != 0)
  {
    fIndex = meanSas;
    SasAtom* m_end = sasAtoms + nAtoms;
    unsigned int* ptr;
    SasAtom* atom;

    for(atom = sasAtoms, ptr = notZero; atom < m_end; atom++, fIndex++, ptr++)
    {
      *fIndex += atom->sas;
      //if(atom->sas >= 0.000001)
      (*ptr)++;
    }

    counter++;
    setStatus(static_cast<float>(counter) / gromacs.getFramesCount());
    (*sasAnalysis) >> sasAtoms;
  }

  {
    unsigned int* ptr;
    for
    (
      fIndex = meanSas, ptr = notZero;
      fIndex < meanSas + nAtoms;
      fIndex++, ptr++
    )
      if(*ptr != 0)
        *fIndex /= *ptr;
  }

  delete[] notZero;
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
  sas = new float[protein.size()];
  setStatus(0);
  counter = 0;
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
      float& curFrame = i->sas.back();

      if(sas[i->getCentralH().index - 1] < 0.000001)
      {
        i->zeros++;
        continue;
      }

      const vector<const Residue*>& residues = i->getResidues();
      float meanGroup = 0;
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
            if(meanSas[k->index - 1] != 0)
              curFrame += sas[k->index - 1];
            meanGroup += meanSas[k->index -1];
          }

      curFrame /= meanGroup;
      if(curFrame < 0.000001)
        i->zeros++;
    }

    counter++;
    setStatus(static_cast<float>(counter) / gromacs.getFramesCount());
    (*sasAnalysis) >> sasAtoms;
  }
  delete sasAnalysis;
  delete[] meanSas;
  delete[] sas;
}

const vector<Pocket>&
Pittpi::getPockets() const
{
  if(sync)
  {
    static vector<Pocket> static_empty_pockets;
    return static_empty_pockets;
  }

  return pockets;
}

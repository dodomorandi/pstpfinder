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

#include <utility>
#include <cassert>
#include <future>

#ifdef HAVE_PYMOD_SADIC
#include "PyIter.h"
#include "PyUtils.h"
#include <Python.h>
#endif

using namespace std;

namespace PstpFinder
{

  Group::Group(const Residue<SasPdbAtom>& refResidue) :
      referenceAtom(&refResidue.getAtomByType("H")), referenceRes(&refResidue)
  {
    zeros = 0;
  }

  Group::Group(const SasPdbAtom& refAtomH) :
      referenceAtom(&refAtomH), referenceRes()
  {
    zeros = 0;
  }

  Group::Group(const SasPdbAtom& refAtomH, const Protein<SasPdbAtom>& protein) :
      referenceAtom(&refAtomH),
      referenceRes(&protein.getResidueByAtom(refAtomH))
  {
    zeros = 0;
  }

  Group::Group(const Group& group) :
      referenceAtom(group.referenceAtom),
      referenceRes(group.referenceRes)
  {
    residues = group.residues;
    sas = group.sas;
    zeros = group.zeros;
  }

  Group::Group(Group&& group) :
      referenceAtom(group.referenceAtom),
      referenceRes(group.referenceRes)
  {
    residues = move(group.residues);
    sas = move(group.sas);
    zeros = move(group.zeros);
  }

  Group::Group(const Group& group, const Protein<SasPdbAtom>& protein) :
      referenceAtom(&protein.getAtomByIndex(group.referenceAtom->index)),
      referenceRes(&protein.getResidueByIndex(group.referenceRes->index))
  {
    for
    (
        auto residue = group.residues.cbegin();
        residue < group.residues.cend();
        residue++
    )
      residues.push_back(&protein.getResidueByIndex((*residue)->index));

    sas = group.sas;
    zeros = group.zeros;
  }

  Group&
  Group::operator =(const Group& group)
  {
    sas = group.sas;
    zeros = group.zeros;
    residues = group.residues;
    referenceAtom = group.referenceAtom;
    referenceRes = group.referenceRes;

    return *this;
  }

  Group&
  Group::operator =(Group&& group)
  {
    sas = move(group.sas);
    zeros = move(group.zeros);
    residues = move(group.residues);
    referenceAtom = group.referenceAtom;
    referenceRes = group.referenceRes;

    return *this;
  }

  Group&
  Group::operator <<(const Residue<SasPdbAtom>& reference)
  {
    residues.push_back(&reference);

    return *this;
  }

  Group&
  Group::operator <<(const Group& group)
  {
    sas = group.sas;
    zeros = group.zeros;
    residues = group.residues;

    return *this;
  }

  const vector<const Residue<SasPdbAtom>*>&
  Group::getResidues() const
  {
    return residues;
  }

  const SasPdbAtom&
  Group::getCentralH() const
  {
    return *referenceAtom;
  }

  const Residue<SasPdbAtom>&
  Group::getCentralRes() const
  {
    if(referenceRes != 0)
      return *referenceRes;
    else
    {
      static Residue<SasPdbAtom> nullRes;
      nullRes.type = AA_UNK;
      return nullRes;
    }
  }

  Pittpi::Pittpi(Gromacs& gromacs, const std::string& sessionFileName,
                 float radius, unsigned long threshold, bool runPittpi) :
      gromacs(gromacs),
      abortFlag(false)
  {
    this->sessionFileName = sessionFileName;
    this->radius = radius;
    this->threshold = threshold;
    sync = true;
    __status = 0;
    averageStructure = gromacs.getAverageStructure();

    if(runPittpi)
      pittpiThread = thread(&Pittpi::pittpiRun, ref(*this));
  }

  Pittpi::~Pittpi()
  {
    join();
  }

  void
  Pittpi::join()
  {
    if(pittpiThread.joinable())
      pittpiThread.join();
  }

  bool
  Pittpi::isFinished()
  {
    bool sync;
    {
      lock_guard<mutex> syncGuard(syncLock);
      sync = this->sync;
    }
    if(pittpiThread.joinable() and not sync)
      pittpiThread.join();
    return not pittpiThread.joinable();
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
    sessionFileName = pittpi.sessionFileName;
    radius = pittpi.radius;
    threshold = pittpi.threshold;
    averageStructure = gromacs.getAverageStructure();
    averageStructure.forceUnlock();
    sync = pittpi.sync;
    __status = pittpi.__status;
    pockets = vector<Pocket>(pittpi.pockets);

    for(auto i = pittpi.groups.cbegin(); i < pittpi.groups.cend(); i++)
    {
      groups.emplace_back(*i, averageStructure);

      for(auto j = pittpi.pockets.cbegin(); j < pittpi.pockets.cend(); j++)
      {
        if(j->group == &(*i))
        {
          pockets.emplace_back(groups.back());
          pockets.back() << *j;
        }
      }
    }
  }

  void
  Pittpi::setStatus(float value) const
  {
    statusMutex.lock();
    __status = value;
    statusMutex.unlock();
    nextStatusCondition.notify_all();
  }

  void
  Pittpi::setStatusDescription(const string& description) const
  {
    statusMutex.lock();
    __statusDescription = description;
    statusMutex.unlock();
  }

  float
  Pittpi::getStatus() const
  {
    if(not pittpiThread.joinable())
      return __status;

    float retval;
    statusMutex.lock();
    retval = __status;
    statusMutex.unlock();
    return retval;
  }

  string
  Pittpi::getStatusDescription() const
  {
    string description;
    statusMutex.lock();
    description = __statusDescription;
    statusMutex.unlock();
    return description;
  }

  void
  Pittpi::waitNextStatus()
  {
    if(isFinished())
      return;

    unique_lock<mutex> lock(nextStatusMutex);
    nextStatusCondition.wait(lock);
  }

  void
  Pittpi::pittpiRun()
  {
    makeGroups(radius);
    if(abortFlag) return;

#define PS_PER_SAS 5
    unsigned int const frameStep = float(PS_PER_SAS) / gromacs.getTimeStep();
    fillGroups(sessionFileName, frameStep);
    if(abortFlag) return;

    sort(groups.begin(), groups.end(),
         [](const Group& a, const Group& b) { return a.zeros > b.zeros; });

    unsigned int noZeroPass = threshold / PS_PER_SAS;
    if(noZeroPass < 20)
      noZeroPass = 0;
    else if(noZeroPass < 40)
      noZeroPass = 1;
    else if(noZeroPass < 60)
      noZeroPass = 2;
    else
      noZeroPass = 3;

    setStatusDescription("Searching for pockets");
    setStatus(0);
    for(vector<Group>::iterator i = groups.begin(); i < groups.end(); i++)
    {
      if(abortFlag) return;
      vector<float>::iterator startPocket = i->sas.end();
      unsigned int notOpenCounter = 0;
      float* maxFrame = 0;
      float mean = 0;

      for(vector<float>::iterator j = i->sas.begin(); j < i->sas.end(); j++)
      {
        if(abortFlag) return;
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
            if(static_cast<float>(distance(startPocket, j - noZeroPass))
               * PS_PER_SAS
               >= threshold)
            {
              Pocket pocket(*i);
              pocket.startFrame = distance(i->sas.begin(), startPocket)
                                  * frameStep
                                  + 1;
              pocket.startPs = distance(i->sas.begin(),
                                        startPocket) * PS_PER_SAS;
              pocket.endFrame = (distance(i->sas.begin(), j) - notOpenCounter
                                 - 1)
                                * frameStep
                                + 1;
              pocket.endPs = (distance(i->sas.begin(), j) - notOpenCounter - 1)
                             * PS_PER_SAS;
              pocket.width = pocket.endPs - pocket.startPs;
              pocket.maxAreaFrame = static_cast<int>(maxFrame
                                                     - &(*i->sas.begin()))
                                    * frameStep
                                    + 1;
              pocket.maxAreaPs = static_cast<float>(maxFrame
                                                    - &(*i->sas.begin()))
                                 * PS_PER_SAS;
              pocket.openingFraction = static_cast<float>(distance(startPocket,
                                                                   j)
                                                          - notOpenCounter
                                                          - 1)
                                       / (i->sas.size() - i->zeros);

              mean /= distance(startPocket, j);
              vector<float>::iterator nearToAverage = startPocket;
              for(vector<float>::iterator k = startPocket + 1; k < j; k++)
                if(abs(mean - *nearToAverage) > abs(mean - *k))
                  nearToAverage = k;
              if(abortFlag) return;
              pocket.averageNearFrame = distance(i->sas.begin(), nearToAverage)
                                        * frameStep
                                        + 1;
              pocket.averageNearPs = distance(i->sas.begin(),
                                              nearToAverage) * PS_PER_SAS;

              pockets.push_back(move(pocket));
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

      setStatus(
          static_cast<float>(distance(groups.begin(), i) + 1) / groups.size());
    }

    sort(pockets.begin(), pockets.end(),
         [](const Pocket& first, const Pocket& second)
         { return first.width > second.width;});
    for(vector<Pocket>::iterator i = pockets.begin(); i < pockets.end(); i++)
    {
      if(abortFlag) return;
      const Residue<SasPdbAtom>& centralRes = i->group->getCentralRes();
      const vector<const Residue<SasPdbAtom>*>& groupRes = i->group->getResidues();
      for(vector<const Residue<SasPdbAtom>*>::const_iterator j = groupRes.begin();
          j < groupRes.end(); j++)
      {
        if(*j == &centralRes or centralRes.index + 1 == (*j)->index
           or centralRes.index - 1 == (*j)->index)
          continue;

        for(vector<Pocket>::iterator k = i + 1; k < pockets.end(); k++)
        {
          if(&k->group->getCentralRes() != &centralRes
             and &k->group->getCentralRes() == *j)
            pockets.erase(k--);
        }
      }
    }

    {
      lock_guard<mutex> syncGuard(syncLock);
      sync = false;
      setStatusDescription("Finished");
      setStatus(1);
    }
  }

  void
  Pittpi::makeGroups(float radius)
  {
    vector<Atom> centers;
    auto& residues = averageStructure.residues();
    radius /= 10.0;

    // Calculate the center for every sidechain (excluding PRO)
    setStatusDescription("Building atom groups");
    setStatus(0);
    centers.reserve(residues.size());
    for(auto i = residues.begin(); i < residues.end(); i++)
    {
      if(abortFlag) return;
      Atom center(0);
      const vector<SasPdbAtom>& atoms = i->atoms;

      if(i->getAtomByType("H1").getTrimmedAtomType() != "UNK" or
          i->type == AA_PRO)
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
      for(auto& atom : atoms)
      {
        string atomType = atom.getTrimmedAtomType();
        if(atomType != "N" and atomType != "CA" and atomType != "H"
           and atomType != "C" and atomType != "O" and atomType != "HA")
        {
          center += atom;
          count++;
        }
      }
      center /= (float) count;

      centers.push_back(center);
      setStatus(
          static_cast<float>(distance(residues.begin(), i) + 1)
          / residues.size());
    }

    groups = makeGroupsByDistance(centers, radius);

#ifdef HAVE_PYMOD_SADIC
    Protein<SasPdbAtom> sadicStructure = runSadic(averageStructure);
    vector<SasPdbAtom> newCenters;

    setStatusDescription("Recalibrating using depth index");
    setStatus(0);
    for(vector<Group>::iterator i = groups.begin(); i < groups.end(); i++)
    {
      if(abortFlag) return;
      const vector<const Residue<SasPdbAtom>*>& groupRes = i->getResidues();
      SasPdbAtom center = i->getCentralH();
      center.x = 0;
      center.y = 0;
      center.z = 0;
      float totalDepth = 0.;

      if(groupRes.size() == 0)
      {
        newCenters.push_back(center);
        continue;
      }

      for(auto& currentGroup : groupRes)
      {
        if(abortFlag) return;
        const vector<SasPdbAtom>& atoms = currentGroup->atoms;

        if(currentGroup->type == AA_PRO)
          continue;

        for(auto& atom : atoms)
        {
          center += atom * sadicStructure.getAtomByIndex(atom.index).bFactor;
          totalDepth += sadicStructure.getAtomByIndex(atom.index).bFactor;
        }
      }

      center /= totalDepth;
      newCenters.push_back(center);
      setStatus(
          static_cast<float>(distance(groups.begin(), i) + 1) / groups.size());
    }
    groups = makeGroupsByDistance(centers, radius, newCenters);
#endif
  }

  vector<Group>
  Pittpi::makeGroupsByDistance(const vector<Atom>& centers, float radius)
  {
    vector<Group> groups;
    auto& residues = averageStructure.residues();

    for(auto i = begin(residues); i < end(residues); i++)
    {
      if(abortFlag) return vector<Group>();
      const SasPdbAtom& hAtom = i->getAtomByType("H");

      if(hAtom.getTrimmedAtomType() == "UNK")
        continue;
      /*
       * NOTE:
       * group(something); group << somethingelse;
       * is different from
       * group(somethingelse)
       * because operator << doesn't change reference residue and/or atom
       */

      Group group(*i);
      group << makeGroupByDistance(centers, hAtom, radius);
      groups.push_back(move(group));
      setStatus(
          static_cast<float>(distance(begin(residues), i) + 1)
          / residues.size());
    }

    return groups;
  }

  vector<Group>
  Pittpi::makeGroupsByDistance(const vector<Atom>& centers, float radius,
                               const vector<SasPdbAtom>& reference)
  {
    vector<Group> groups;
    auto& residues = averageStructure.residues();

    auto refIterator = begin(reference);
    for(auto resIterator = begin(residues); resIterator < end(residues);
        resIterator++, refIterator++)
    {
      if(abortFlag) return vector<Group>();
      if(resIterator->getAtomByType("H").getTrimmedAtomType() == "UNK")
      {
        refIterator--;
        continue;
      }

      Group group(*resIterator);
      group << makeGroupByDistance(centers, *refIterator, radius);
      groups.push_back(move(group));
      setStatus(
          static_cast<float>(distance(begin(residues), resIterator) + 1)
          / residues.size());
    }

    return groups;
  }

  Group
  Pittpi::makeGroupByDistance(const vector<Atom>& centers,
                              const SasPdbAtom& atom, float radius)
  {
    auto& residues = averageStructure.residues();
    Group group(atom);

    if(atom.getTrimmedAtomType() == "UNK")
      return group;

    for(vector<Atom>::const_iterator j = begin(centers); j < end(centers);
        j++)
    {
      if(abortFlag) return group;
      const Residue<SasPdbAtom>& curResidue = residues[distance(begin(centers),
                                                                j)];
      if(curResidue.type == AA_PRO)
        continue;

      if(atom.distance(*j) <= radius)
        group << curResidue;
    }

    return group;
  }

  void
  Pittpi::fillGroups(const string& sessionFileName, unsigned int timeStep)
  {
    float* sas;
    float* meanSas;
    float* fIndex;
    SasAtom* sasAtoms = 0;
    SasAnalysis<ifstream>* sasAnalysis;
    unsigned int counter = 0;

    const float frames = gromacs.getFramesCount();
    vector<int> protein = gromacs.getGroup("Protein");
    const int nAtoms = protein.size();

    meanSas = new float[nAtoms]();
    unsigned int* notZero = new unsigned int[nAtoms]();

    /* First of all we need to calculate SAS means */
    setStatusDescription("Calculating SAS means");
    setStatus(0);
    sasAnalysis = new SasAnalysis<ifstream>(gromacs, sessionFileName);
    (*sasAnalysis) >> sasAtoms;
    while(sasAtoms != 0)
    {
      if(abortFlag) return;

      fIndex = meanSas;
      SasAtom* m_end = sasAtoms + nAtoms;
      unsigned int* ptr;
      SasAtom* atom;

      for(atom = sasAtoms, ptr = notZero; atom < m_end; atom++, fIndex++, ptr++)
        *fIndex += atom->sas;

      counter++;
      setStatus(static_cast<float>(counter) / gromacs.getFramesCount());
      (*sasAnalysis) >> sasAtoms;
    }

    {
      unsigned int* ptr;
      for(fIndex = meanSas, ptr = notZero; fIndex < meanSas + nAtoms;
          fIndex++, ptr++)
        *fIndex /= counter;
    }

    delete[] notZero;
    delete sasAnalysis;

    /* Let's prepare groups sas vectors */
    for(vector<Group>::iterator i = begin(groups); i < end(groups); i++)
      i->sas.reserve(frames);

    /* Now we have to normalize values and store results per group */
    sas = new float[protein.size()];
    float* sasCounter = new float[protein.size()];
    setStatusDescription("Searching for zeros and normalizing SAS");
    setStatus(0);
    counter = 0;
    sasAnalysis = new SasAnalysis<ifstream>(gromacs, sessionFileName);
    (*sasAnalysis) >> sasAtoms;
    while(sasAtoms != 0)
    {
      if(abortFlag) return;

      /*
       * This part is a "legacy" method. It have been implemented in perl time ago
       * and needs refactoring. The main problem is math related, because we have to
       * find a good solution to take "consecutively opened pocket" above a certain
       * threshold. With every frame (and every SAS value) it could be not so easy
       * to develop a GOOD algorithm. For now we implement only the old method used
       * until now.
       *
       * 28 sep 2011: Only now I understand that I need data binning to obtain the
       * same results as the original algorithm. This must be done BEFORE
       * normalization!
       * -- Edoardo Morandi
       */

      if(counter % timeStep == 0)
      {
        for(float* i = sasCounter; i < sasCounter + protein.size(); i++)
          *i = 0.0;
      }

      fIndex = sas;
      for(SasAtom* atom = sasAtoms; atom < sasAtoms + protein.size();
          atom++, fIndex++)
        *fIndex = atom->sas;

      {
        float* i;
        float* j;
        for(i = sasCounter, j = sas; i < sasCounter + protein.size(); i++, j++)
          *i += *j;
      }

      if((counter + 1) % timeStep == 0)
      {
        for(float* i = sasCounter; i < sasCounter + protein.size(); i++)
          *i /= timeStep;

        if(abortFlag) return;

        for(vector<Group>::iterator i = groups.begin(); i < groups.end(); i++)
        {
          if(abortFlag) return;
          i->sas.push_back(0);
          float& curFrame = i->sas.back();

          if(sasCounter[i->getCentralH().index - 1] < 0.000001)
          {
            i->zeros++;
            continue;
          }

          const vector<const Residue<SasPdbAtom>*>& residues = i->getResidues();
          for(auto j = begin(residues); j < end(residues); j++)
          {
            if(abortFlag) return;
            const SasPdbAtom& atomH = (*j)->getAtomByType("H");
            if(atomH.getTrimmedAtomType() == "UNK")
              continue;

            if(meanSas[atomH.index - 1] != 0)
              curFrame += sasCounter[atomH.index - 1]
                          / meanSas[atomH.index - 1];
          }

          curFrame /= i->getResidues().size();

          if(curFrame < 0.000001)
            i->zeros++;

        }
      }

      counter++;
      setStatus(static_cast<float>(counter) / gromacs.getFramesCount());
      (*sasAnalysis) >> sasAtoms;
    }

    if(counter % timeStep != 0)
    {
      for(float* i = sasCounter; i < sasCounter + protein.size(); i++)
        *i /= (counter % timeStep);

      if(abortFlag) return;

      for(vector<Group>::iterator i = groups.begin(); i < groups.end(); i++)
      {
        if(abortFlag) return;
        i->sas.push_back(0);
        float& curFrame = i->sas.back();

        if(sasCounter[i->getCentralH().index - 1] < 0.000001)
        {
          i->zeros++;
          continue;
        }

        const vector<const Residue<SasPdbAtom>*>& residues = i->getResidues();
        for(vector<const Residue<SasPdbAtom>*>::const_iterator j = residues.begin();
            j < residues.end(); j++)
        {
          if(abortFlag) return;
          const SasPdbAtom& atomH = (*j)->getAtomByType("H");
          if(atomH.getTrimmedAtomType() == "UNK")
            continue;

          if(meanSas[atomH.index - 1] != 0)
            curFrame += sasCounter[atomH.index - 1] / meanSas[atomH.index - 1];
        }

        curFrame /= i->getResidues().size();

        if(curFrame < 0.000001)
          i->zeros++;
      }
    }

    delete[] sasCounter;
    delete sasAnalysis;
    delete[] meanSas;
    delete[] sas;
  }

#ifdef HAVE_PYMOD_SADIC
  Protein<SasPdbAtom>
  Pittpi::runSadic(const Protein<SasPdbAtom>& structure) const
  {
    Protein<SasPdbAtom> sadicProtein;
    Pdb<SasPdbAtom> pdb;
    pdb.proteins.push_back(structure);

    setStatusDescription("Running Sadic");
    setStatus(-1);
    Py_Initialize();
    PyObject* oSadicRunner = PyImport_ImportModule("sadic.runner");
    PyObject* oSadicCmdline = PyImport_ImportModule("sadic.cmdline");
    PyObject* oSadicIos = PyImport_ImportModule("sadic.ios");
    PyObject* oSadicViewer = PyImport_ImportModule("sadic.viewer");

    pdb.write("/tmp/sadic_in.pdb");

    PyObject* oSettings = py::callMethod(oSadicCmdline, "parse_command_line", "[ssss]", "--all-atoms", "-f", "pdb", "/tmp/sadic_in.pdb");
    PyObject* oOutput = py::callMethod(oSadicIos, "get_output", "O", oSettings);
    PyObject* oKeywords = Py_BuildValue("{s:O}", "settings", oSettings);
    PyObject* emptyTuple = PyTuple_New(0);
    python_iterable iterModels {PyObject_Call(PyObject_GetAttrString(oSadicRunner, "iter_models"), emptyTuple, oKeywords)};
    Py_DECREF(emptyTuple);

    static const std::function<PyObject*(python_iterable&)> getFirst {[](python_iterable& iter){return *std::begin(iter);}};
    std::future<PyObject*> firstModel = std::async(std::launch::async, getFirst, std::ref(iterModels));
    while(firstModel.wait_for(std::chrono::milliseconds(80)) != std::future_status::ready)
        setStatus(-1);

    PyObject* oViewers = PyTuple_GetItem(firstModel.get(), 1);

    Py_DECREF(py::callMethod(oOutput, "output", "O", oViewers));
    python_iterable viewers {oViewers};
    std::future<PyObject*> firstViewer = std::async(std::launch::async, getFirst, std::ref(viewers));
    while(firstViewer.wait_for(std::chrono::milliseconds(80)) != std::future_status::ready)
        setStatus(-1);
    
    PyObject* oFilename = py::callMethod(oOutput, "mangle_file_name", "O", firstViewer.get());
    sadicProtein = move(Pdb<SasPdbAtom>(py::toCString(oFilename)).proteins[0]);

    Py_DECREF(oOutput);
    Py_DECREF(oKeywords);
    Py_DECREF(oSettings);

    Py_DECREF(oSadicViewer);
    Py_DECREF(oSadicIos);
    Py_DECREF(oSadicCmdline);
    Py_DECREF(oSadicRunner);

    Py_Finalize();

    return sadicProtein;
  }
#endif /* HAVE_PYMOD_SADIC */

  void
  Pittpi::abort()
  {
    abortFlag = true;
    join();
    sync = true;
    nextStatusCondition.notify_all();
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

  Pittpi::SerializablePockets::SerializablePockets(
      const vector<Pocket>& pockets, const vector<Group>& groups)
  {
    for(auto& pocket : pockets)
    {
      assert(pocket.group >= groups.data() and
             pocket.group < groups.data() + groups.size());

      SerializablePocket serializablePocket(pocket);
      serializablePocket.groupIndex = pocket.group - groups.data();
      this->pockets.push_back(move(serializablePocket));
    }
  }

  Pittpi::SerializableGroups::SerializableGroups(const vector<Group>& groups,
                                                 const Protein<SasPdbAtom>& protein)
  {
#ifndef NDEBUG
    const vector<const SasPdbAtom*>& atoms(protein.atoms());
    const vector<Residue<SasPdbAtom>>& residues(protein.residues());
#else
    // FIXME: useless argument if not in debug mode?
    (void) protein;
#endif

    for(auto& group : groups)
    {
      SerializableGroup serializableGroup(group);

#ifndef NDEBUG
      if(serializableGroup.referenceAtom != nullptr)
      {
        bool atomFound(false);
        for(auto& atom : atoms)
        {
          if(atom == serializableGroup.referenceAtom)
          {
            atomFound = true;
            break;
          }
        }
        assert(atomFound);
      }
#endif

      assert(serializableGroup.referenceRes == nullptr or
             (serializableGroup.referenceRes >= residues.data() and
              serializableGroup.referenceRes <
              residues.data() + residues.size()));

      serializableGroup.referenceAtomIndex = serializableGroup.referenceAtom
          ->index;
      serializableGroup.referenceResIndex = serializableGroup.referenceRes
          ->index;

      for(auto& residue : serializableGroup.residues)
      {
        assert(residue >= residues.data()
               and residue < residues.data() + residues.size());

        serializableGroup.residuesIndex.push_back(residue->index);
      }

      this->groups.push_back(move(serializableGroup));
    }
  }

  void
  Pittpi::SerializableGroups::updateGroups(vector<Group>& groupsToUpdate,
                                           const Protein<SasPdbAtom>& protein) const
  {
    groupsToUpdate.clear();
    groupsToUpdate.reserve(groups.size());

    for(auto& serializableGroup : groups)
    {
      Group group(protein.getResidueByIndex(serializableGroup.referenceResIndex));
      group.sas = move(serializableGroup.sas);
      group.zeros = serializableGroup.zeros;

      for(auto& residueIndex : serializableGroup.residuesIndex)
        group << protein.getResidueByIndex(residueIndex);

      groupsToUpdate.push_back(move(group));
    }
  }

  void
  Pittpi::SerializablePockets::updatePockets(
      vector<Pocket>& pocketsToUpdate, const vector<Group>& groups) const
  {
    pocketsToUpdate.clear();
    pocketsToUpdate.reserve(pockets.size());

    for(auto& serializablePocket : pockets)
    {
      Pocket pocket(groups[serializablePocket.groupIndex]);
      pocket << serializablePocket;

      pocketsToUpdate.push_back(move(pocket));
    }
  }
}

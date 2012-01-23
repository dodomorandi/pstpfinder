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
#include <utility>

#include <boost/interprocess/sync/scoped_lock.hpp>
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>

using namespace std;
namespace py = boost::python;

namespace PstpFinder
{

  Group::Group(const Residue& refResidue) :
      referenceAtom(&refResidue.getAtomByType("H")), referenceRes(&refResidue)
  {
    zeros = 0;
  }

  Group::Group(const PdbAtom& refAtomH) :
      referenceAtom(&refAtomH), referenceRes()
  {
    zeros = 0;
  }

  Group::Group(const PdbAtom& refAtomH, const Protein& protein) :
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
  Group::operator <<(const Residue& reference)
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

  Pittpi::Pittpi(Gromacs& gromacs, const std::string& sessionFileName,
                 float radius, unsigned long threshold) :
      gromacs(gromacs),
      abortFlag(false)
  {
    this->sessionFileName = sessionFileName;
    this->radius = radius;
    this->threshold = threshold;
    sync = true;
    __status = 0;

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
      Group group(*i);

      const vector<const Residue*> residues(i->getResidues());
      for(auto j = residues.cbegin(); j < residues.cend(); j++)
        group << averageStructure.getResidueByIndex((*j)->index);
      groups.push_back(move(group));

      for(auto j = pittpi.pockets.cbegin(); j < pittpi.pockets.cend(); j++)
      {
        if(j->group == &(*i))
        {
          pockets.emplace_back(groups.back());
          pockets.back() = *j;
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
    averageStructure = gromacs.getAverageStructure();
    makeGroups(radius);
    if(abortFlag) return;

#define PS_PER_SAS 5
    unsigned int const frameStep = float(PS_PER_SAS) / gromacs.getTimeStep();
    fillGroups(sessionFileName, frameStep);
    if(abortFlag) return;

    sort(groups.begin(), groups.end(), Group::sortByZeros);

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
            if(static_cast<float>(distance(startPocket, j)) * PS_PER_SAS
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

    sort(pockets.begin(), pockets.end(), Pocket::sortByWidth);
    for(vector<Pocket>::iterator i = pockets.begin(); i < pockets.end(); i++)
    {
      if(abortFlag) return;
      const Residue& centralRes = i->group->getCentralRes();
      const vector<const Residue*>& groupRes = i->group->getResidues();
      for(vector<const Residue*>::const_iterator j = groupRes.begin();
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

    ofstream pocketLog("/tmp/pockets.log", ios::out | ios::trunc);
    ofstream pocketDetailLog("/tmp/pockets_details.log", ios::out | ios::trunc);

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
    const vector<Residue>& residues = averageStructure.residues();
    radius /= 10.0;

    // Calculate the center for every sidechain (excluding PRO)
    setStatusDescription("Building atom groups");
    setStatus(0);
    centers.reserve(residues.size());
    for(vector<Residue>::const_iterator i = residues.begin();
        i < residues.end(); i++)
    {
      if(abortFlag) return;
      Atom center(0);
      const vector<PdbAtom>& atoms = i->atoms;

      if(strcmp(i->getAtomByType("H1").type, "UNK") != 0 or i->type == AA_PRO)
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
      for(vector<PdbAtom>::const_iterator j = atoms.begin(); j < atoms.end();
          j++)
      {
        if(strcmp(j->type, "N") != 0 and strcmp(j->type, "CA") != 0
           and strcmp(j->type, "H") != 0
           and strcmp(j->type, "C") != 0
           and strcmp(j->type, "O") != 0
           and strcmp(j->type, "HA") != 0)
        {
          center += *j;
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
#if HAVE_PYMOD_SADIC == 1
    Protein sadicStructure = runSadic(averageStructure);
    vector<PdbAtom> newCenters;

    setStatusDescription("Recalibrating using depth index");
    setStatus(0);
    for(vector<Group>::iterator i = groups.begin(); i < groups.end(); i++)
    {
      if(abortFlag) return;
      const vector<const Residue*>& groupRes = i->getResidues();
      PdbAtom center = i->getCentralH();
      center.x = 0;
      center.y = 0;
      center.z = 0;
      float totalDepth = 0.;

      if(groupRes.size() == 0)
      {
        newCenters.push_back(center);
        continue;
      }

      for(vector<const Residue*>::const_iterator j = groupRes.begin();
          j < groupRes.end(); j++)
      {
        if(abortFlag) return;
        const vector<PdbAtom>& atoms = (*j)->atoms;

        if((*j)->type == AA_PRO)
          continue;

        for(vector<PdbAtom>::const_iterator k = atoms.begin(); k < atoms.end();
            k++)
        {
          center += (*k) * sadicStructure.getAtomByIndex(k->index).bFactor;
          totalDepth += sadicStructure.getAtomByIndex(k->index).bFactor;
        }
      }

      center /= totalDepth;
      newCenters.push_back(center);
      setStatus(
          static_cast<float>(distance(groups.begin(), i) + 1) / groups.size());
    }
    groups = makeGroupsByDistance(centers, radius, newCenters);
#endif
#endif
  }

  vector<Group>
  Pittpi::makeGroupsByDistance(const vector<Atom>& centers, float radius)
  {
    vector<Group> groups;
    const vector<Residue>& residues = averageStructure.residues();

    for(vector<Residue>::const_iterator i = residues.begin();
        i < residues.end(); i++)
    {
      if(abortFlag) return vector<Group>();
      const PdbAtom& hAtom = i->getAtomByType("H");

      if(strcmp(hAtom.type, "UNK") == 0)
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
          static_cast<float>(distance(residues.begin(), i) + 1)
          / residues.size());
    }

    return groups;
  }

  vector<Group>
  Pittpi::makeGroupsByDistance(const vector<Atom>& centers, float radius,
                               const vector<PdbAtom>& reference)
  {
    vector<Group> groups;
    const vector<Residue>& residues = averageStructure.residues();

    vector<Residue>::const_iterator resIterator;
    vector<PdbAtom>::const_iterator refIterator;
    for(resIterator = residues.begin(), refIterator = reference.begin();
        resIterator < residues.end(); resIterator++, refIterator++)
    {
      if(abortFlag) return vector<Group>();
      if(strcmp(resIterator->getAtomByType("H").type, "UNK") == 0)
      {
        refIterator--;
        continue;
      }
      Group group(*resIterator);
      group << makeGroupByDistance(centers, *refIterator, radius);
      groups.push_back(move(group));
      setStatus(
          static_cast<float>(distance(residues.begin(), resIterator) + 1)
          / residues.size());
    }

    return groups;
  }

  Group
  Pittpi::makeGroupByDistance(const vector<Atom>& centers, const PdbAtom& atom,
                              float radius)
  {
    const vector<Residue>& residues = averageStructure.residues();
    Group group(atom);

    if(strcmp(atom.type, "UNK") == 0)
      return group;

    for(vector<Atom>::const_iterator j = centers.begin(); j < centers.end();
        j++)
    {
      if(abortFlag) return group;
      const Residue& curResidue = residues[distance(centers.begin(), j)];
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
    SasAnalysis* sasAnalysis;
    unsigned int counter = 0;

    const float frames = gromacs.getFramesCount();
    vector<int> protein = gromacs.getGroup("Protein");
    const int nAtoms = protein.size();

    meanSas = new float[nAtoms]();
    unsigned int* notZero = new unsigned int[nAtoms]();

    /* First of all we need to calculate SAS means */
    setStatusDescription("Calculating SAS means");
    setStatus(0);
    sasAnalysis = new SasAnalysis(gromacs, sessionFileName, false);
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
    for(vector<Group>::iterator i = groups.begin(); i < groups.end(); i++)
      i->sas.reserve(frames);

    /* Now we have to normalize values and store results per group */
    sas = new float[protein.size()];
    float* sasCounter = new float[protein.size()];
    setStatusDescription("Searching for zeros and normalizing SAS");
    setStatus(0);
    counter = 0;
    sasAnalysis = new SasAnalysis(gromacs, sessionFileName, false);
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

          const vector<const Residue*>& residues = i->getResidues();
          for(vector<const Residue*>::const_iterator j = residues.begin();
              j < residues.end(); j++)
          {
            if(abortFlag) return;
            const PdbAtom& atomH = (*j)->getAtomByType("H");
            if(strcmp(atomH.type, "UNK") == 0)
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

        const vector<const Residue*>& residues = i->getResidues();
        for(vector<const Residue*>::const_iterator j = residues.begin();
            j < residues.end(); j++)
        {
          if(abortFlag) return;
          const PdbAtom& atomH = (*j)->getAtomByType("H");
          if(strcmp(atomH.type, "UNK") == 0)
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
  Protein
  Pittpi::runSadic(const Protein& structure) const
  {
    Protein sadicProtein;
    Py_Initialize();

    {
      setStatusDescription("Running Sadic");
      setStatus(-1);
      py::object sadic = py::import("sadic");
      py::object setting = py::import("sadic.setting");
      py::object viewer = py::import("sadic.viewer");
      py::object cmdline = py::import("sadic.cmdline");
      py::stl_input_iterator<py::object> iterObjEnd;

      viewer.attr("load_plugin")();
      py::list queries;

      structure.dumpPdb("/tmp/sadic_in.pdb");

      py::list argv;
      argv.append("/tmp/sadic_in.pdb");
      py::object settings = cmdline.attr("parse_command_line")(argv);
      settings.attr("entity_spec") = "/tmp/sadic_in.pdb";
      settings.attr("all_atoms") = true;
      settings.attr("file_out") = "/tmp/sadic_out.pdb";
      settings.attr("output_format") = sadic.attr("consts").attr("OUTPUT_PDB");
      settings.attr("quiet") = true;
      settings.attr("atom_name") = py::object();

      py::object out = sadic.attr("get_output")(settings);
      py::object reader = sadic.attr("get_reader")(settings);
      py::object scheme = sadic.attr("get_sampling_scheme")(settings);

      py::list models_viewers;

      py::object file = sadic.attr("iter_files")(settings).attr("next")();

      py::object models = reader.attr("get_models")(file);
      unsigned int imodel = 0;
      for(py::stl_input_iterator<py::object> model(models); model != iterObjEnd;
          model++, imodel++)
      {
        if(abortFlag) return sadicProtein;

        setStatus(-1);
        py::object query = sadic.attr("get_query")(settings, *model);
        py::object prot = sadic.attr("Protein")();
        prot.attr("add_atoms")(*model);

        py::object viewers = viewer.attr("create_viewers")(settings, query);
        models_viewers.append(viewers);

        if(imodel == 0)
        {
          for(;;)
          {
            if(abortFlag) return sadicProtein;
            setStatus(-1);
            query.attr("sample")(prot, scheme);

            if(settings.attr("radius")
               != sadic.attr("consts").attr("RADIUS_NO_INSIDE"))
              break;

            py::stl_input_iterator<py::object> i(query);
            for(; i != iterObjEnd; i++)
            {
              if(abortFlag) return sadicProtein;
              setStatus(-1);
              if(static_cast<bool>(i->attr("all_inside")))
                break;
            }
            if(i == iterObjEnd)
              break;

            float step = py::extract<float>(settings.attr("step"));
            scheme.attr("grow")(step);
          }
        }
        else
          query.attr("sample")(prot, scheme);

        out.attr("output")(viewers);
        queries.append(query);
      }

      file.attr("close")();
      setStatus(-1);

      py::object total_viewers = viewer.attr("create_total_viewers")(
          settings, models_viewers);
      for(py::stl_input_iterator<py::object> viewers(total_viewers);
          viewers != iterObjEnd; viewers++)
      {
        out.attr("output")(*viewers);
        setStatus(-1);
      }

      if(py::len(models_viewers) == 0)
        return Protein();

      for(py::stl_input_iterator<py::object> model_viewer(models_viewers);
          model_viewer != iterObjEnd; model_viewer++)
      {
        for(py::stl_input_iterator<py::object> curViewer(*model_viewer);
            curViewer != iterObjEnd; curViewer++)
        {
          if(abortFlag) return sadicProtein;
          string fileName = py::extract<string>(
              out.attr("mangle_file_name")(*curViewer));
          sadicProtein = Protein(fileName);
          break;
        }
        break; // Just the first protein of the first file... for now!
      }
    }
    Py_Finalize();

    return sadicProtein;
  }

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
}
#endif

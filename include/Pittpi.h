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

#ifndef _PITTPI_H
#define _PITTPI_H

#include "config.h"
#include "Gromacs.h"
#include <vector>
#include <boost/thread/thread.hpp>
#include <boost/interprocess/sync/interprocess_mutex.hpp>
#include <boost/interprocess/sync/interprocess_condition.hpp>

namespace Gromacs
{
  class Group
  {
  public:
    Group(const Residue& refResidue);
    Group(const PdbAtom& refAtomH);
    Group(const PdbAtom& refAtomH, const Protein& protein);
    Group& operator <<(const Residue& value);
    Group& operator =(const Group& group);
    const vector<const Residue*>& getResidues() const;
    const PdbAtom& getCentralH() const;
    const Residue& getCentralRes() const;
    static bool sortByZeros(const Group& a, const Group& b);

    vector<float> sas;
    unsigned int zeros;
  private:
    const PdbAtom* const referenceAtom;
    const Residue* const referenceRes;
    vector<const Residue*> residues;
  };

  struct Pocket
  {
    const Group& group;
    unsigned int startFrame;
    float startPs;
    unsigned int endFrame;
    float endPs;
    unsigned int width;
    float openingFraction;
    unsigned int averageNearFrame;
    float averageNearPs;
    unsigned int maxAreaFrame;
    float maxAreaPs;

    Pocket(const Group& group) : group(group) { ; }
    Pocket& operator =(const Pocket& pocket)
    {
      startFrame = pocket.startFrame;
      startPs = pocket.startPs;
      endFrame = pocket.endFrame;
      endPs = pocket.endPs;
      width = pocket.width;
      openingFraction = pocket.openingFraction;
      averageNearFrame = pocket.averageNearFrame;
      averageNearPs = pocket.averageNearPs;
      maxAreaFrame = pocket.maxAreaFrame;
      maxAreaPs = pocket.maxAreaPs;

      return *this;
    }
  };

  /**
   * @brief Related to the method developed by Matteo De Chiara and Silvia
   * Bottini
   *
   * The original code have been developed in Perl. A code refactory have been
   * done to obtain a more usable and stable program.
   *
   * @note This is only a stub. The class must be rewritten, but for now
   * the main objective remains the productivity of the software. I created
   * a class only to remember what must be done in the future. Possibly soon.
   * @author Edoardo Morandi
   */
  class Pittpi
  {
  public:
    Pittpi(Gromacs& gromacs,
           const std::string& sasAnalysisFileName,
           float radius,
           unsigned long threshold);
    ~Pittpi();

    void join();
    void setStatus(float value) const;
    float getStatus() const;
    void waitNextStatus();
    bool isFinished();
  private:
    std::vector<Group> makeGroups(float radius);
    void fillGroups(std::vector<Group>& groups,
                    const string& sasAnalysisFileName);
    std::vector<Group> makeGroupsByDistance(const std::vector<Atom>& centers,
                                            float radius);
    std::vector<Group> makeGroupsByDistance(
                                        const std::vector<Atom>& centers,
                                        float radius,
                                        const std::vector<PdbAtom>& reference);
    Group makeGroupByDistance(const std::vector<Atom>& centers,
                              const PdbAtom& atom,
                              float radius);
    void pittpiRun();
#ifdef HAVE_PYMOD_SADIC
    Protein runSadic(const Protein& structure) const;
#endif

    Gromacs* p_gromacs;
    std::string sasAnalysisFileName;
    float radius;
    unsigned long threshold;
    Protein averageStructure;
    boost::thread pittpiThread;
    mutable boost::interprocess::interprocess_mutex statusMutex;
    mutable boost::interprocess::interprocess_mutex nextStatusMutex;
    mutable boost::interprocess::interprocess_condition nextStatusCondition;
    mutable float __status;
    bool sync;
  };
}

#endif

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

#include "Gromacs.h"
#include <vector>

namespace Gromacs
{

  class Group
  {
  public:
    Group(const Residue& refResidue);
    Group(const PdbAtom& refAtomH);
    Group& operator <<(const Residue& value);
    const vector<const Residue*>& getResidues() const;
    const PdbAtom& getCentralH() const;
    static bool sortByZeros(const Group& a, const Group& b);

    vector<float> sas;
    unsigned int zeros;
  private:
    const PdbAtom* reference;
    vector<const Residue*> residues;
  };

  struct Pocket
  {
    const Group* group;
    unsigned int startFrame;
    unsigned int endFrame;
    float openingFraction;
    unsigned int meanNearFrame;
    unsigned int maxAreaFrame;
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
  private:
    std::vector<Group> makeGroups(float radius);
    void fillGroups(std::vector<Group>& groups,
                    const string& sasAnalysisFileName);

    Gromacs* m_gromacs;
    Protein averageStructure;
  };
}

#endif

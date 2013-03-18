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

#ifndef SPACECUBE_H_
#define SPACECUBE_H_

#include <array>
#include <bitset>
#include <list>

using namespace std;

namespace PstpFinder
{
  class SpaceCube
  {
    public:
      const array<float, 3> base;
      bitset<8> flags;
      list<unsigned long> involvedAtomsIndices;

      SpaceCube(float baseX, float baseY, float baseZ) noexcept :
          base({{baseX, baseY, baseZ}}) {}
      inline float x() const noexcept { return base[0]; }
      inline float y() const noexcept { return base[1]; }
      inline float z() const noexcept { return base[2]; }
  };
} /* namespace PstpFinder */
#endif /* SPACECUBE_H_ */

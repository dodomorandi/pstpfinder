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

#ifndef POINTRADIUS_H_
#define POINTRADIUS_H_

using namespace std;

namespace PstpFinder
{
  namespace MarchingCubes
  {
    struct PointRadius
    {
        PointRadius() = default;
        template<typename AtomType>
        PointRadius(const AtomType& atom, float radius) :
          x(atom.x), y(atom.y), z(atom.z), radius(radius) {}

        float x, y, z, radius;
    };
  }
}

#endif /* POINTRADIUS_H_ */

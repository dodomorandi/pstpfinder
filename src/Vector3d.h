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

#ifndef VECTOR3D_H_
#define VECTOR3D_H_

using namespace std;

#include <array>
#include <initializer_list>

namespace PstpFinder
{
  namespace MarchingCubes
  {
    struct Vector3d
    {
        Vector3d() = default;
        Vector3d(float x, float y, float z) : coords{{x, y, z}} {};
        Vector3d(initializer_list<float> init)
          { copy(begin(init), end(init), begin(coords)); }
        array<float, 3> coords;

        float& x() { return coords[0]; }
        float& y() { return coords[1]; }
        float& z() { return coords[2]; }

        const float& x() const { return coords[0]; }
        const float& y() const { return coords[1]; }
        const float& z() const { return coords[2]; }
    };

    struct PointAndNormal
    {
        array<float, 6> coords;

        float& x() { return coords[0]; }
        float& y() { return coords[1]; }
        float& z() { return coords[2]; }
        float& nx() { return coords[3]; }
        float& ny() { return coords[4]; }
        float& nz() { return coords[5]; }

        const float& x() const { return coords[0]; }
        const float& y() const { return coords[1]; }
        const float& z() const { return coords[2]; }
        const float& nx() const { return coords[3]; }
        const float& ny() const { return coords[4]; }
        const float& nz() const { return coords[5]; }
    };
  }
}

#endif /* VECTOR3D_H_ */

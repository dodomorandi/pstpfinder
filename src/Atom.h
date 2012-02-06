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

#ifndef _PSTPFINDER_ATOM_H
#define _PSTPFINDER_ATOM_H

#include <cmath>

extern "C"
{
#include <typedefs.h>
}

namespace PstpFinder
{
  struct Atom
  {
      real x, y, z;

      Atom()
      {
        ;
      }

      Atom(float value)
      {
        x = value;
        y = value;
        z = value;
      }

      Atom
      operator +(const Atom& value) const
      {
        Atom ret;
        ret.x = x + value.x;
        ret.y = y + value.y;
        ret.z = z + value.z;

        return ret;
      }

      Atom
      operator /(float value) const
      {
        Atom ret;
        ret.x = x / value;
        ret.y = y / value;
        ret.z = z / value;

        return ret;
      }

      Atom
      operator *(float value) const
      {
        Atom ret;
        ret.x = x * value;
        ret.y = y * value;
        ret.z = z * value;

        return ret;
      }

      Atom&
      operator +=(const Atom& value)
      {
        x += value.x;
        y += value.y;
        z += value.z;

        return *this;
      }

      Atom&
      operator /=(float value)
      {
        x /= value;
        y /= value;
        z /= value;

        return *this;
      }

      Atom&
      operator *=(float value)
      {
        x *= value;
        y *= value;
        z *= value;

        return *this;
      }

      float
      distance(const Atom& atom) const
      {
        return std::sqrt(std::pow(x - atom.x, 2) + std::pow(y - atom.y, 2)
            + std::pow(z - atom.z, 2));
      }
  };
}

#endif

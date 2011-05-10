#ifndef _ATOM_H
#define _ATOM_H

#include <cmath>

extern "C" {
  #include <typedefs.h>
}

namespace Gromacs
{
  struct Atom
  {
    real x, y, z;

    Atom() { return; }
    Atom(float value) { x = value; y = value; z = value; }

    Atom operator +(const Atom& value)
    {
      Atom ret;
      ret.x = x + value.x;
      ret.y = y + value.y;
      ret.z = z + value.z;

      return ret;
    }

    Atom operator /(float value)
    {
      Atom ret;
      ret.x = x / value;
      ret.y = y / value;
      ret.z = z / value;

      return ret;
    }

    Atom& operator +=(const Atom& value)
    {
      x += value.x;
      y += value.y;
      z += value.z;

      return *this;
    }

    Atom& operator /=(float value)
    {
      x /= value;
      y /= value;
      z /= value;

      return *this;
    }

    float distance(const Atom& atom) const
    {
      return std::sqrt(std::pow(x - atom.x, 2) +
                       std::pow(y - atom.y, 2) +
                       std::pow(z - atom.z, 2));
    }
  };
};

#endif

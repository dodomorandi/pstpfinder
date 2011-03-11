#ifndef _ATOM_H
#define _ATOM_H

extern "C" {
  #include <typedefs.h>
}

namespace Gromacs
{
  struct Atom
  {
    real x, y, z;
  };
};

#endif

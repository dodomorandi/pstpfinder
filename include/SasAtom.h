#ifndef _SASATOM_H
#define _SASATOM_H

#include "Atom.h"
extern "C" {
   #include <typedefs.h>
}

namespace Gromacs
{
  struct SasAtom: public Atom
  {
     real sas;
  };
};

#endif

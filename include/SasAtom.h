#ifndef _SASATOM_H
#define _SASATOM_H

#include "Atom.h"
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

extern "C" {
   #include <typedefs.h>
}

namespace Gromacs
{
  struct SasAtom: public Atom
  {
    real sas;
  };
  
  class SasAtomSerializable: SasAtom
  {
  public:
    SasAtomSerializable(SasAtom& sasAtom)
    {
      x = sasAtom.x;
      y = sasAtom.y;
      z = sasAtom.z;
      sas = sasAtom.sas;
    }
  private:
    friend class boost::serialization::access;
    
    template <typename Archive>
    void
    serialize(Archive& ar, const unsigned int version)
    {
      ar & x;
      ar & y;
      ar & z;
      ar & sas;
    }
  };
};

#endif

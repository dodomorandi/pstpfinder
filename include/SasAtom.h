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

#ifndef _SASATOM_H
#define _SASATOM_H

#include "Atom.h"
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

extern "C"
{
#include <typedefs.h>
}

namespace Gromacs
{
  struct SasAtom :
      public Atom
  {
      real sas;

    private:
      friend class boost::serialization::access;

      template<typename Archive>
        void
        serialize(Archive& ar, const unsigned int version)
        {
          ar & x;
          ar & y;
          ar & z;
          ar & sas;
        }
  };
}

#endif

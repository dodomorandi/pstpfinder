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
#include "Pdb.h"
#include "Serializer.h"

extern "C"
{
#include <gromacs/typedefs.h>
}

namespace PstpFinder
{
  struct SasPdbAtom : public PdbAtom
  {
      SasPdbAtom() : PdbAtom(), sas(0) {}
      SasPdbAtom(const SasPdbAtom& atom) = default;
      SasPdbAtom(SasPdbAtom&& atom) = default;
      SasPdbAtom& operator =(const SasPdbAtom& atom) = default;
      SasPdbAtom& operator =(SasPdbAtom&& atom) = default;

      explicit SasPdbAtom(const PdbAtom& atom) : PdbAtom(atom), sas(0) {}
      explicit SasPdbAtom(PdbAtom&& atom) : PdbAtom(move(atom)), sas(0) {}
      explicit SasPdbAtom(const ProteinAtom& atom) : PdbAtom(atom), sas(0) {}
      explicit SasPdbAtom(ProteinAtom&& atom) : PdbAtom(move(atom)), sas(0) {}

      /*
      template<typename T,
                typename enable_if<is_base_of<Atom,
                    typename remove_reference<T>::type>::value and not
                  is_same<SasPdbAtom,
                    typename remove_reference<T>::type>::value>
                ::type* = nullptr>
      explicit SasPdbAtom(T&& t) : PdbAtom(forward<T>(t)), sas(0) {}
      template<typename T,
                typename enable_if<is_base_of<Atom,
                    typename remove_reference<T>::type>::value and not
                  is_same<SasPdbAtom,
                    typename remove_reference<T>::type>::value>
                ::type* = nullptr>
      SasPdbAtom& operator =(T&& t)
      {
        ProteinAtom::operator =(forward<T>(t));
        if(is_base_of<SasPdbAtom, typename remove_reference<T>::type>::value)
        {
          const SasPdbAtom& tmp = static_cast<SasPdbAtom>(t);
          sas = tmp.sas;
        }

        return *this;
      }
      */


      explicit SasPdbAtom(int index) :
        PdbAtom(index), sas(0) {}

      float sas;
  };

  struct SasAtom : public Atom
  {
      real sas;

    private:
      template<typename, typename>
      friend class Serializer;

      template<typename Serializer>
      void
      serialize(Serializer serializer)
      {
        serializer & x;
        serializer & y;
        serializer & z;
        serializer & sas;
      }
  };
}

#endif

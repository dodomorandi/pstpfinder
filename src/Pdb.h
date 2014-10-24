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

#ifndef _PDB_H
#define _PDB_H

#include "Protein.h"
#include "utils.h"
#include <vector>
#include <tuple>
#include <type_traits>

namespace PstpFinder
{
  struct PdbAtom : public ProteinAtom
  {
      PdbAtom() : ProteinAtom(), bFactor(0), occupancy(0) {}
      PdbAtom(const PdbAtom& atom) = default;
      PdbAtom(PdbAtom&& atom) = default;
      PdbAtom& operator =(const PdbAtom& atom) = default;
      PdbAtom& operator =(PdbAtom&& atom) = default;

      explicit PdbAtom(const ProteinAtom& atom) :
          ProteinAtom(atom), bFactor(0), occupancy(0) {}
      explicit PdbAtom(ProteinAtom&& atom) :
          ProteinAtom(std::move(atom)), bFactor(0), occupancy(0) {}

      explicit PdbAtom(int index) :
          ProteinAtom(index), bFactor(0), occupancy(0) {}

      float bFactor;
      float occupancy;
  };

  template<typename AtomType = PdbAtom>
  class Pdb
  {
    public:
      Pdb() = default;
      Pdb(const std::string& fileName);
      Pdb(const Protein<AtomType>& protein) { buildPdb(protein);}
      Pdb(Protein<AtomType>&& protein) { buildPdb(std::move(protein));}

      template<typename Stream>
      Pdb(Stream&& stream, typename std::enable_if<is_stream_base_of<
          std::basic_istream, Stream>::value>::type* = nullptr);

      void write(const std::string& filename) const;

      template<typename Stream>
      typename std::enable_if<
          is_stream_base_of<std::basic_istream, Stream>::value
          or is_stream_base_of<std::basic_ostream, Stream>::value>::type
      write(Stream& stream) const;
      std::vector<Protein<AtomType>> proteins;

    private:
      template<typename Protein>
        void buildPdb(Protein&& protein);

      template<typename Stream>
      typename std::enable_if<is_stream_base_of<
          std::basic_istream, Stream>::value>::type
      readFromStream(Stream&& stream);

      inline std::tuple<float, float>
      getAuxParameters(const PdbAtom& atom) const
      {
        return std::make_tuple(atom.bFactor, atom.occupancy);
      }

      inline std::tuple<float, float>
      getAuxParameters(const Atom& atom) const
      {
        return std::tuple<float, float>({0., 1.});
      }
  };
}

#include "Pdb.cpp"

#endif

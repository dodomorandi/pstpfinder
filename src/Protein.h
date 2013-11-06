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

#ifndef _PROTEIN_H
#define _PROTEIN_H

#include "Atom.h"
#include "utils.h"
#include <string>
#include <vector>
#include <algorithm>

namespace PstpFinder
{
  enum Aminoacids
  {
    AA_UNK,
    AA_ALA,
    AA_VAL,
    AA_GLY,
    AA_LEU,
    AA_ILE,
    AA_THR,
    AA_SER,
    AA_MET,
    AA_CYS,
    AA_PRO,
    AA_TYR,
    AA_TRP,
    AA_PHE,
    AA_HIS,
    AA_ARG,
    AA_LYS,
    AA_GLU,
    AA_ASP,
    AA_GLN,
    AA_ASN
  };

  const std::string aminoacidLetter[] =
    { "X", "A", "V", "G", "L", "I", "T", "S", "M", "C", "P", "Y", "W", "F",
      "H", "R", "K", "E", "D", "Q", "N" };

  const std::string aminoacidTriplet[] =
    { "UNK", "ALA", "VAL", "GLY", "LEU", "ILE", "THR", "SER", "MET", "CYS",
      "PRO", "TYR", "TRP", "PHE", "HIS", "ARG", "LYS", "GLU", "ASP", "GLN",
      "ASN" };
  const std::string aminoacidUncommonTranslator[] =
    { "CYP", "CYS", "CYD", "CYS", "HID", "HIS", "HIE", "HIS", "HIP", "HIS",
      "LYP", "LYS", "LYD", "LYS", "LYN", "LYS", "ASH", "ASP", "GLH", "GLU" };
  const unsigned int aminoacidUncommonTranslatorSize =
      sizeof(aminoacidUncommonTranslator)
          / sizeof(*aminoacidUncommonTranslator);

  struct ProteinAtom : public Atom
  {
      char type[5];
      unsigned int index;

      ProteinAtom() = default;
      explicit ProteinAtom(const std::string& type)
      {
        if(type.size() > 4)
          type.copy(this->type, 4);
        else
          type.copy(this->type, std::string::npos);
        this->type[4] = '\0';
      }

      explicit ProteinAtom(int index) : index(index) {}

      inline bool isType(const std::string& type) const
      {
        return type == this->type;
      }

      std::string getTrimmedType() const
      {
        std::string trimmedType(std::begin(type), std::end(type));
        trimmedType.erase(remove(std::begin(trimmedType), std::end(trimmedType), ' '),
                          std::end(trimmedType));
        return trimmedType;
      }
  };

  template<typename AtomType = ProteinAtom>
  struct Residue
  {
      Aminoacids type;
      int index;
      std::vector<AtomType> atoms;
      char chain;

      Residue() = default;
      explicit Residue(Aminoacids aminoacid) : type(aminoacid) {}

      template<typename T,
        typename std::enable_if<std::is_same<typename std::remove_reference<T>::type,
          Residue<AtomType>>::value>::type* = nullptr>
      Residue(const typename std::remove_reference<T>::type& residue) :
        type(residue.type),
        index(residue.index),
        atoms(residue.atoms),
        chain(residue.chain) {}

      template<typename T,
        typename std::enable_if<std::is_same<typename std::remove_reference<T>::type,
          Residue<AtomType>>::value>::type* = nullptr>
      Residue(typename std::remove_reference<T>::type&& residue) :
        type(move(residue.type)),
        index(move(residue.index)),
        atoms(move(residue.atoms)),
        chain(move(residue.chain)) {}

      // TODO: a template template splitter
      template<template <typename> class OldResidue, typename OldAtom,
        typename std::enable_if<not std::is_same<
          typename std::remove_reference<OldResidue<OldAtom>>::type,
          Residue<AtomType>>::value>::type* = nullptr>
      Residue(OldResidue<OldAtom>&& residue) :
        type(move(residue.type)), index(move(residue.index)),
        atoms(std::begin(residue.atoms), std::end(residue.atoms)),
        chain(move(residue.chain))
      {
          static_assert(
              std::is_same<typename std::remove_reference<OldResidue<OldAtom>>::type,
              Residue<OldAtom>>::value, "Argument type is not Residue");
      }
      template<template <typename> class OldResidue, typename OldAtom,
        typename std::enable_if<not std::is_same<
          typename std::remove_reference<OldResidue<OldAtom>>::type,
          Residue<AtomType>>::value>::type* = nullptr>
      Residue(const OldResidue<OldAtom>& residue) :
        type(residue.type), index(residue.index),
        atoms(std::begin(residue.atoms), std::end(residue.atoms)),
        chain(residue.chain)
      {
          static_assert(
              std::is_same<typename std::remove_reference<OldResidue<OldAtom>>::type,
              Residue<OldAtom>>::value, "Argument type is not Residue");
      }

      template<typename T, typename std::enable_if<std::is_same<
          typename std::remove_reference<T>::type,
          Residue<AtomType>>::value>::type* = nullptr>
      Residue& operator =(const typename std::remove_reference<T>::type& residue)
      {
        type = residue.type;
        index = residue.index;
        atoms = residue.atoms;
        chain = residue.chain;

        return *this;
      }

      template<typename T, typename std::enable_if<std::is_same<
          typename std::remove_reference<T>::type,
          Residue<AtomType>>::value>::type* = nullptr>
      Residue& operator =(typename std::remove_reference<T>::type&& residue)
      {
        type = move(residue.type);
        index = move(residue.index);
        atoms = move(residue.atoms);
        chain = move(residue.chain);

        return *this;
      }

      template<template <typename> class OldResidue, typename OldAtom,
        typename std::enable_if<not std::is_same<
          typename std::remove_reference<OldResidue<OldAtom>>::type,
          Residue<AtomType>>::value>::type* = nullptr>
      Residue& operator =(OldResidue<OldAtom>&& residue);
      template<template <typename> class OldResidue, typename OldAtom,
        typename std::enable_if<not std::is_same<
          typename std::remove_reference<OldResidue<OldAtom>>::type,
          Residue<AtomType>>::value>::type* = nullptr>
      Residue& operator =(const OldResidue<OldAtom>& residue);

      const AtomType& getAtomByType(const std::string& atomType) const;
      static Aminoacids getTypeByName(std::string residueName);
  };

  template<typename AtomType = ProteinAtom>
  class Protein
  {
    public:
      typedef AtomType atom_type;
      std::string name;
      int model;

      Protein();
      template<typename OldAtom>
        Protein(const Protein<OldAtom>& protein);
      template<typename OldAtom>
        Protein& operator =(const Protein<OldAtom>& proteinr);

      const std::vector<Residue<AtomType>>& residues() const;
      std::vector<Residue<AtomType>>& residuesRW();
      std::vector<const AtomType*>& atoms() const;
      template<typename ResidueType>
        bool appendResidue(ResidueType&& residue,typename std::enable_if<
          std::is_same<typename std::remove_reference<ResidueType>::type,
          Residue<AtomType>>::value>::type* = nullptr);

      const AtomType& getAtomByIndex(unsigned int index) const;
      const Residue<AtomType>& getResidueByAtom(int atomIndex) const;
      const Residue<AtomType>& getResidueByAtom(const AtomType& atom) const;
      const Residue<AtomType>& getResidueByIndex(int index) const;
      void lock() const;
      void forceUnlock() const;

      template<typename NewAtomType>
      typename std::enable_if<not std::is_same<AtomType, NewAtomType>::value,
        std::vector<Residue<NewAtomType>>>::type
        convertResidues() const;
      template<typename NewAtomType>
      typename std::enable_if<std::is_same<AtomType, NewAtomType>::value,
        std::vector<Residue<NewAtomType>>>::type
        convertResidues() const;

    protected:
      std::vector<Residue<AtomType>> pResidues;
      mutable std::vector<const AtomType*> pAtoms;
      mutable bool locked;
  };
}

#include "Protein.cpp"

#endif

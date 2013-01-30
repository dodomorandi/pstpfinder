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

using namespace std;

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

  const string aminoacidLetter[] =
    { "X", "A", "V", "G", "L", "I", "T", "S", "M", "C", "P", "Y", "W", "F",
      "H", "R", "K", "E", "D", "Q", "N" };

  const string aminoacidTriplet[] =
    { "UNK", "ALA", "VAL", "GLY", "LEU", "ILE", "THR", "SER", "MET", "CYS",
      "PRO", "TYR", "TRP", "PHE", "HIS", "ARG", "LYS", "GLU", "ASP", "GLN",
      "ASN" };
  const string aminoacidUncommonTranslator[] =
    { "CYP", "CYS", "CYD", "CYS", "HID", "HIS", "HIE", "HIS", "HIP", "HIS",
      "LYP", "LYS", "LYD", "LYS", "LYN", "LYS", "ASH", "ASP", "GLH", "GLU" };
  const unsigned int aminoacidUncommonTranslatorSize =
      sizeof(aminoacidUncommonTranslator)
          / sizeof(*aminoacidUncommonTranslator);

  struct ProteinAtom :
      public Atom
  {
      char type[5];
      unsigned int index;

      ProteinAtom() = default;
      explicit ProteinAtom(const string& type)
      {
        if(type.size() > 4)
          type.copy(this->type, 4);
        else
          type.copy(this->type, string::npos);
        this->type[4] = '\0';
      }

      explicit ProteinAtom(int index) : index(index) {}

      inline bool isType(const string& type) const
      {
        return type == this->type;
      }

      string getTrimmedType() const
      {
        string trimmedType(begin(type), end(type));
        trimmedType.erase(remove(begin(trimmedType), end(trimmedType), ' '),
                          end(trimmedType));
        return trimmedType;
      }
  };

  template<typename AtomType = ProteinAtom>
  struct Residue
  {
      Aminoacids type;
      int index;
      vector<AtomType> atoms;
      char chain;

      Residue() = default;
      explicit Residue(Aminoacids aminoacid) : type(aminoacid) {}

      template<typename T,
        typename enable_if<is_same<typename remove_reference<T>::type,
          Residue<AtomType>>::value>::type* = nullptr>
      Residue(const typename remove_reference<T>::type& residue) :
        type(residue.type),
        index(residue.index),
        atoms(residue.atoms),
        chain(residue.chain) {}

      template<typename T,
        typename enable_if<is_same<typename remove_reference<T>::type,
          Residue<AtomType>>::value>::type* = nullptr>
      Residue(typename remove_reference<T>::type&& residue) :
        type(move(residue.type)),
        index(move(residue.index)),
        atoms(move(residue.atoms)),
        chain(move(residue.chain)) {}

      // TODO: a template template splitter
      template<template <typename> class OldResidue, typename OldAtom,
        typename enable_if<not is_same<
          typename remove_reference<OldResidue<OldAtom>>::type,
          Residue<AtomType>>::value>::type* = nullptr>
      Residue(OldResidue<OldAtom>&& residue) :
        type(move(residue.type)), index(move(residue.index)),
        atoms(begin(residue.atoms), end(residue.atoms)),
        chain(move(residue.chain))
      {
          static_assert(
              is_same<typename remove_reference<OldResidue<OldAtom>>::type,
              Residue<OldAtom>>::value, "Argument type is not Residue");
      }
      template<template <typename> class OldResidue, typename OldAtom,
        typename enable_if<not is_same<
          typename remove_reference<OldResidue<OldAtom>>::type,
          Residue<AtomType>>::value>::type* = nullptr>
      Residue(const OldResidue<OldAtom>& residue) :
        type(residue.type), index(residue.index),
        atoms(begin(residue.atoms), end(residue.atoms)),
        chain(residue.chain)
      {
          static_assert(
              is_same<typename remove_reference<OldResidue<OldAtom>>::type,
              Residue<OldAtom>>::value, "Argument type is not Residue");
      }

      template<typename T, typename enable_if<is_same<
          typename remove_reference<T>::type,
          Residue<AtomType>>::value>::type* = nullptr>
      Residue& operator =(const typename remove_reference<T>::type& residue)
      {
        type = residue.type;
        index = residue.index;
        atoms = residue.atoms;
        chain = residue.chain;

        return *this;
      }

      template<typename T, typename enable_if<is_same<
          typename remove_reference<T>::type,
          Residue<AtomType>>::value>::type* = nullptr>
      Residue& operator =(typename remove_reference<T>::type&& residue)
      {
        type = move(residue.type);
        index = move(residue.index);
        atoms = move(residue.atoms);
        chain = move(residue.chain);

        return *this;
      }

      template<template <typename> class OldResidue, typename OldAtom,
        typename enable_if<not is_same<
          typename remove_reference<OldResidue<OldAtom>>::type,
          Residue<AtomType>>::value>::type* = nullptr>
      Residue& operator =(OldResidue<OldAtom>&& residue);
      template<template <typename> class OldResidue, typename OldAtom,
        typename enable_if<not is_same<
          typename remove_reference<OldResidue<OldAtom>>::type,
          Residue<AtomType>>::value>::type* = nullptr>
      Residue& operator =(const OldResidue<OldAtom>& residue);

      const AtomType& getAtomByType(const string& atomType) const;
      static Aminoacids getTypeByName(string residueName);
      bool isComplete() const noexcept;

    private:
      bool checkAtomTypePresence(const vector<string>& atomTypes) const noexcept
        { return __checkAtomTypePresence(atomTypes);}
      bool checkAtomTypePresence(vector<string>&& atomTypes) const noexcept
        { return __checkAtomTypePresence(atomTypes);}
      template<typename VectorOfStrings>
        bool __checkAtomTypePresence(VectorOfStrings&& atomTypes) const noexcept;
  };

  template<typename AtomType = ProteinAtom>
  class Protein
  {
    public:
      typedef AtomType atom_type;
      string name;
      int model;

      Protein();
      template<typename OldAtom>
        Protein(const Protein<OldAtom>& protein);
      template<typename OldAtom>
        Protein& operator =(const Protein<OldAtom>& proteinr);

      const vector<Residue<AtomType>>& residues() const;
      vector<Residue<AtomType>>& residuesRW();
      vector<const AtomType*>& atoms() const;
      template<typename ResidueType>
        bool appendResidue(ResidueType&& residue,typename enable_if<
          is_same<typename remove_reference<ResidueType>::type,
          Residue<AtomType>>::value>::type* = nullptr);

      const AtomType& getAtomByIndex(unsigned int index) const;
      const Residue<AtomType>& getResidueByAtom(int atomIndex) const;
      const Residue<AtomType>& getResidueByAtom(const AtomType& atom) const;
      const Residue<AtomType>& getResidueByIndex(int index) const;
      void lock() const;
      void forceUnlock() const;

      template<typename NewAtomType>
      typename enable_if<not is_same<AtomType, NewAtomType>::value,
        vector<Residue<NewAtomType>>>::type
        convertResidues() const;
      template<typename NewAtomType>
      typename enable_if<is_same<AtomType, NewAtomType>::value,
        vector<Residue<NewAtomType>>>::type
        convertResidues() const;
      void removeUnknownResidues();
      void remakeIndices();

    protected:
      vector<Residue<AtomType>> pResidues;
      mutable vector<const AtomType*> pAtoms;
      mutable bool locked;
  };
}

#include "Protein.cpp"

#endif

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
#include <array>

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
    { "CYP", "CYS", "CYD", "CYS", "CYS2", "CYS", "HID", "HIS", "HIE", "HIS",
      "HIP", "HIS", "LYP", "LYS", "LYD", "LYS", "LYN", "LYS", "LYSN", "LYS",
      "ASH", "ASP", "ASPH", "ASP", "GLH", "GLU", "GLUH", "GLU", "QLN", "GLN",
      "ARGN", "ARG" };
  const unsigned int aminoacidUncommonTranslatorSize =
      sizeof(aminoacidUncommonTranslator)
          / sizeof(*aminoacidUncommonTranslator);

  struct ProteinAtom : public Atom
  {
      std::array<char, 2> name;
      std::array<char, 2> type;
      unsigned int index;

      ProteinAtom() = default;
      explicit ProteinAtom(const std::string& type) { setAtomType(type); }
      explicit ProteinAtom(int index) : index(index) {}

      std::string getAtomType() const
      {
        std::stringstream ss;
        ss << std::setiosflags(std::ios_base::right) << std::setfill(' ');
        ss << std::setw(2) << array2string(name);
        ss << std::resetiosflags(std::ios_base::right) << array2string(type);

        return ss.str();
      }

      std::string getTrimmedAtomType() const
      {
        return array2string(name) + array2string(type);
      }

      void setAtomType(const std::string& type, bool pdbFormat = true)
      {
        size_t typeSize = type.size();
        if(pdbFormat)
        {
          if(typeSize > 4)
            typeSize = 4;

          if(typeSize >= 2)
          {
            auto endIter = std::copy_if(
                std::begin(type),
                std::min(std::begin(type) + 2, std::end(type)),
                std::begin(name),
                [](char a) { return a != ' ';});
            if(endIter < std::end(name))
              *endIter = '\0';
          }
          else if(typeSize == 1)
          {
            name[0] = type[0];
            name[1] = '\0';
          }

          auto endIter = std::copy_if(std::min(begin(type) + 2, std::end(type)),
                                      std::end(type), std::begin(this->type),
                                      [](char a) { return a != ' ';});
          if(endIter < std::end(this->type))
            *endIter = '\0';
        }
        else
        {
          /*
           * We don't have indications about name-chain convention.
           * We try to evaluate atom names basing on standard atoms in proteins
           */
          switch(type[0])
          {
            case 'H':
            case 'C':
            case 'N':
            case 'O':
            case 'S':
            case 'P':
              name[0] = type[0];
              name[1] = '\0';
              {
                auto endIter = std::copy_if(
                    std::min(begin(type) + 1, std::end(type)), std::end(type),
                    std::begin(this->type),
                    [](char a) { return a != ' '; });

                if(endIter < std::end(this->type))
                  *endIter = '\0';
              }

              break;
            default:
              /* Don't know what to do. Switching to default */
              setAtomType(type, true);
              break;
          }
        }
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
        Protein& operator =(const Protein<OldAtom>& protein);

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

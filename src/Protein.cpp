/*  This file is part of PSTP-finder, an user friendly tool to analyze GROMACS
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

#include "Protein.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <locale>
#include <algorithm>
#include <cassert>

using namespace PstpFinder;
using namespace std;

template<typename AtomType>
template<template <typename> class OldResidue, typename OldAtom,
  typename enable_if<not is_same<
    typename remove_reference<OldResidue<OldAtom>>::type,
    Residue<AtomType>>::value>::type*>
Residue<AtomType>&
Residue<AtomType>::operator =(const OldResidue<OldAtom>& residue)
{
  static_assert(
      is_same<typename remove_reference<OldResidue<OldAtom>>::type,
      Residue<OldAtom>>::value, "Argument type is not Residue");

  type = residue.type;
  index = residue.index;
  atoms.reserve(residue.atoms.size());
  for(auto atom : residue.atoms)
    atoms.push_back(atom);

  return *this;
}

template<typename AtomType>
template<template <typename> class OldResidue, typename OldAtom,
  typename enable_if<not is_same<
    typename remove_reference<OldResidue<OldAtom>>::type,
    Residue<AtomType>>::value>::type*>
Residue<AtomType>&
Residue<AtomType>::operator =(OldResidue<OldAtom>&& residue)
{
  static_assert(
      is_same<typename remove_reference<OldResidue<OldAtom>>::type,
      Residue<OldAtom>>::value, "Argument type is not Residue");

  type = move(residue.type);
  index = move(residue.index);
  atoms.reserve(residue.atoms.size());
  for(auto atom : residue.atoms)
    atoms.push_back(move(atom));

  return *this;
}

template<typename AtomType>
const AtomType&
Residue<AtomType>::getAtomByType(const string& atomType) const
{
  static const AtomType unknown = static_cast<AtomType>(ProteinAtom("UNK"));
  for(const AtomType& atom : atoms)
    if(atom.getTrimmedType() == atomType)
      return atom;

  return unknown;
}

template<typename AtomType>
Aminoacids
Residue<AtomType>::getTypeByName(string residueName)
{
  transform(begin(residueName), end(residueName),
            begin(residueName), ::toupper);

  for(unsigned int i = 0; i < 21; i++)
  {
    for(unsigned int j = 0; j < 12; j += 2)
    {
      if(residueName == aminoacidUncommonTranslator[j] or
          residueName == "N" + aminoacidUncommonTranslator[j] or
          residueName == "C" + aminoacidUncommonTranslator[j])
      {
        residueName = aminoacidUncommonTranslator[j + 1];
        break;
      }
    }

    if(residueName == aminoacidTriplet[i] or
        residueName == "N" + aminoacidTriplet[i] or
        residueName == "C" + aminoacidTriplet[i])
      return static_cast<Aminoacids> (i);;
  }

  return AA_UNK;
}

template<typename AtomType>
Protein<AtomType>::Protein()
{
  locked = false;
  model = 0;
}

template<typename AtomType>
template<typename OldAtom>
Protein<AtomType>::Protein(const Protein<OldAtom>& protein) :
  name(protein.name),
  model(protein.model),
  pResidues(move(protein.convertResidues<AtomType>()))
{
  locked = false;
}

template<typename AtomType>
template<typename OldAtom>
Protein<AtomType>&
Protein<AtomType>::operator =(const Protein<OldAtom>& protein)
{
  pResidues = move(protein.convertResidues<AtomType>());
  locked = false;
  name = protein.name;
  model = protein.model;

  return *this;
}

template<typename AtomType>
const vector<Residue<AtomType>>&
Protein<AtomType>::residues() const
{
  return pResidues;
}

template<typename AtomType>
vector<Residue<AtomType>>&
Protein<AtomType>::residuesRW()
{
  locked = false;
  return pResidues;
}

template<typename AtomType>
vector<const AtomType*>&
Protein<AtomType>::atoms() const
{
  if(not locked)
  {
    locked = true;
    pAtoms.clear();
    for(auto& residue : pResidues)
      for(const AtomType& atom : residue.atoms)
        pAtoms.push_back(&atom);
  }

  return pAtoms;
}

template<typename AtomType>
template<typename ResidueType>
bool
Protein<AtomType>::appendResidue(ResidueType&& residue,
    typename enable_if<is_same<typename remove_reference<ResidueType>::type,
    Residue<AtomType>>::value>::type*)
{
  if(locked)
    return false;

  pResidues.push_back(forward<Residue<AtomType>>(residue));
  return true;
}

template<typename AtomType>
const AtomType&
Protein<AtomType>::getAtomByIndex(unsigned int index) const
{
  static AtomType unknown(-1);
  for(const AtomType* atom : pAtoms)
  {
    if(atom->index == index)
      return *atom;
  }

  return unknown;
}

template<typename AtomType>
const Residue<AtomType>&
Protein<AtomType>::getResidueByAtom(int atomIndex) const
{
  return getResidueByAtom(getAtomByIndex(atomIndex));
}

template<typename AtomType>
const Residue<AtomType>&
Protein<AtomType>::getResidueByAtom(const AtomType& atom) const
{
  static const Residue<AtomType> unknown(AA_UNK);

    for(auto residueIter = begin(pResidues); residueIter != end(pResidues);
        residueIter++)
    {
      if(residueIter->atoms[0].index > atom.index)
      {
        residueIter--;
          for(auto atomIter = begin(residueIter->atoms);
              atomIter != end(residueIter->atoms); atomIter++)
          if(atomIter->index == atom.index)
            return *residueIter;

        return unknown;
      }
    }

  return unknown;
}

template<typename AtomType>
const Residue<AtomType>&
Protein<AtomType>::getResidueByIndex(int index) const
{
  static const Residue<AtomType> unknown(AA_UNK);

  if(index >= 1 and index - 1 < static_cast<int>(pResidues.size())
     and pResidues[index - 1].index == index)
    return pResidues[index - 1];
  else if(index >= 0 and index < static_cast<int>(pResidues.size())
          and pResidues[index].index == index)
    return pResidues[index];
  else
  {
    for(auto& residue : pResidues)
    {
      if(residue.index == index)
        return residue;
    }

    return unknown;
  }
}

template<typename AtomType>
void
Protein<AtomType>::lock() const
{
  if(not locked)
  {
    (void) atoms();
    locked = true;
  }
}

template<typename AtomType>
void
Protein<AtomType>::forceUnlock() const
{
  locked = false;
}

template<typename AtomType>
template<typename NewAtomType>
typename enable_if<not is_same<AtomType, NewAtomType>::value,
  vector<Residue<NewAtomType>>>::type
Protein<AtomType>::convertResidues() const
{
  vector<Residue<NewAtomType>> newResidues;
  newResidues.reserve(pResidues.size());

  for(auto& residue : pResidues)
    newResidues.emplace_back(residue);

  return newResidues;
}

template<typename AtomType>
template<typename NewAtomType>
typename enable_if<is_same<AtomType, NewAtomType>::value,
  vector<Residue<NewAtomType>>>::type
Protein<AtomType>::convertResidues() const
{
  return vector<Residue<NewAtomType>>(pResidues);
}

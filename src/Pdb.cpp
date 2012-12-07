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

#include "Pdb.h"
#include "SasAtom.h"
#include <utility>
#include <cassert>

using namespace PstpFinder;
using namespace std;

template<typename AtomType>
Pdb<AtomType>::Pdb(const string& fileName)
{
  ifstream pdbFile(fileName);
  // FIXME: I should move pdbFile, but I can't until streams
  //        will be moveable
  readFromStream(pdbFile);
}

template<typename AtomType>
template<typename Stream>
Pdb<AtomType>::Pdb(Stream&& stream, typename enable_if<
                   is_stream_base_of<basic_istream, Stream>::value>::type*)
{
  readFromStream(forward<Stream>(stream));
}

template<typename AtomType>
Pdb<AtomType>::Pdb(Protein<AtomType>&& protein)
{
  proteins.clear();
  proteins.push_back(forward<Protein<AtomType>>(protein));
}

template<typename AtomType>
void
Pdb<AtomType>::write(const string& filename) const
{
  ofstream pdb(filename, ios_base::out | ios_base::trunc);
  write(pdb);
}

template<typename AtomType>
template<typename Stream>
  typename enable_if<is_base_of<base_stream
  (basic_istream, Stream), Stream>::value
                     or is_base_of<base_stream
                     (basic_ostream, Stream), Stream>::value>::type
Pdb<AtomType>::write(Stream& stream) const
{
  assert(stream.is_open());

  bool writeModel = true;
  if(proteins.size() <= 1)
    writeModel = false;

  for(auto& protein : proteins)
  {
    if(writeModel)
      stream << "MODEL     " << setw(4) << protein.model << endl;
    for(auto& residue : protein.residues())
    {
      for(auto& atom : residue.atoms)
      {
        string atomType(atom.type);
        atomType.erase(atomType.begin() + 1, atomType.end());
        stream << "ATOM  " << setw(5) << atom.index << " ";
        stream << setiosflags(ios::left) << setw(4) << atom.type << " ";
        stream << setw(3) << aminoacidTriplet[residue.type] << " A";
        stream << resetiosflags(ios::left) << setw(4) << residue.index;
        stream << "    " << setiosflags(ios::right);
        stream << setiosflags(ios::fixed);
        stream << setprecision(3) << setw(8) << (atom.x * 10.);
        stream << setprecision(3) << setw(8) << (atom.y * 10.);
        stream << setprecision(3) << setw(8) << (atom.z * 10.);
        stream << setprecision(2) << setw(6);

        auto aux = getAuxParameters(atom);
        stream << ((get<1>(aux) >= 100) ? 99.99 : get<1>(aux));
        stream << setprecision(2) << setw(6);
        stream << ((get<0>(aux) >= 100) ? 99.99 : get<0>(aux));
        stream << resetiosflags(ios::fixed);
        stream << "            " << setw(2) << atomType;
        stream << resetiosflags(ios::right) << endl;
      }
    }
  }
  stream.close();
}

template<typename AtomType>
template<typename Stream>
typename enable_if<is_stream_base_of<basic_istream, Stream>::value>::type
Pdb<AtomType>::readFromStream(Stream&& stream)
{
  stringstream streamLine;
  Protein<AtomType>* protein = nullptr;
  Residue<AtomType> residue;
  string buffer;
  int modelIndex = 0;
  Aminoacids residueType;

  streamLine.exceptions(ios_base::failbit);
  residue.index = 0;

  while(not stream.eof())
  {
    getline(stream, buffer);
    streamLine.clear();
    streamLine.str(buffer);
    int residueIndex;

      streamLine >> setw(6) >> buffer;
      if(buffer == "ATOM")
      {
        SasPdbAtom atom;

        // Mandatory atom section
        try
        {
          if(protein == nullptr)
          {
            protein = new Protein<AtomType>();
            protein->model = ++modelIndex;
          }

          streamLine >> setw(5) >> atom.index;
          streamLine.seekg(1, ios_base::cur); // Empty space
          streamLine >> noskipws >> setw(5) >> atom.type >> skipws;
          streamLine.seekg(1, ios_base::cur);// Alternate location indicator
          streamLine >> setw(3) >> buffer;
          residueType = Residue<AtomType>::getTypeByName(buffer);
          streamLine.seekg(2, ios_base::cur);// Empty space + Chain ID
          streamLine >> setw(4) >> residueIndex;
          streamLine.seekg(4, ios_base::cur);// Code for insertion
                                             // of residues + 3 spaces
          streamLine >> setw(8) >> atom.x;
          streamLine >> setw(8) >> atom.y;
          streamLine >> setw(8) >> atom.z;

          atom.x /= 10.;
          atom.y /= 10.;
          atom.z /= 10.;
        }
        catch(ios_base::failure& fail)
        {
          continue;
        }

        // Non-mandatory atom section
        try
        {
          streamLine >> setw(6) >> atom.occupancy;
        }
        catch(ios_base::failure& fail)
        {
          atom.occupancy = 0;
        }

        try
        {
          streamLine >> setw(6) >> atom.bFactor;
        }
        catch(ios_base::failure& fail)
        {
          atom.bFactor = 0;
        }
        // Then element and charge... but I don't mind of them

        if(residue.index == 0)
        {
          residue.type = residueType;
          residue.index = residueIndex;
        }
        else if(residueIndex != residue.index)
        {
          protein->appendResidue(move(residue));
          residue = Residue<AtomType>();
          residue.type = residueType;
          residue.index = residueIndex;
        }

        residue.atoms.push_back(move(atom));
      }
      else if(buffer.substr(0,5) == "MODEL")
      {
        if(protein)
        {
          protein->appendResidue(move(residue));
          protein->lock();
          proteins.push_back(move(*protein));
        }
        protein = new Protein<AtomType>();
        residue = Residue<AtomType>();
        residue.index = 0;

        stringstream mdlIndex;
        mdlIndex << buffer.substr(6);
        mdlIndex >> modelIndex;
        protein->model = modelIndex;
      }
      else if(protein and buffer == "ENDMDL")
      {
        protein->appendResidue(move(residue));
        protein->lock();
        proteins.push_back(move(*protein));
        protein = nullptr;
        residue = Residue<AtomType>();
        residue.index = 0;
      }
  }

  if(residue.atoms.size() > 0)
  {
    protein->appendResidue(move(residue));
  }

  if(protein and protein->residues().size() > 0)
  {
    protein->lock();
    proteins.push_back(move(*protein));
  }

  stream.close();
}

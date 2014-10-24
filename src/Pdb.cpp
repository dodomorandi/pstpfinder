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
#include "Protein.h"
#include <utility>
#include <cassert>
#include <sstream>

namespace PstpFinder
{

template<typename AtomType>
Pdb<AtomType>::Pdb(const std::string& fileName)
{
  std::ifstream pdbFile(fileName);
  // FIXME: I should move pdbFile, but I can't until streams
  //        will be moveable
  readFromStream(pdbFile);
}

template<typename AtomType>
template<typename Stream>
Pdb<AtomType>::Pdb(Stream&& stream, typename std::enable_if<
                   is_stream_base_of<std::basic_istream, Stream>::value>::type*)
{
  readFromStream(std::forward<Stream>(stream));
}

template<typename AtomType>
template<typename Protein>
void
Pdb<AtomType>::buildPdb(Protein&& protein)
{
  proteins.clear();
  proteins.push_back(std::forward<Protein>(protein));
}

template<typename AtomType>
void
Pdb<AtomType>::write(const std::string& filename) const
{
  std::ofstream pdb(filename, std::ios_base::out | std::ios_base::trunc);
  write(pdb);
}

template<typename AtomType>
template<typename Stream>
  typename std::enable_if<is_stream_base_of<std::basic_istream, Stream>::value
                     or is_stream_base_of<std::basic_ostream, Stream>::value>::type
Pdb<AtomType>::write(Stream& stream) const
{
  assert(stream.is_open());

  bool writeModel = true;
  if(proteins.size() <= 1)
    writeModel = false;

  for(auto& protein : proteins)
  {
    if(writeModel)
      stream << "MODEL     " << std::setw(4) << protein.model << std::endl;
    for(auto& residue : protein.residues())
    {
      for(auto& atom : residue.atoms)
      {
        std::string atomName(array2string(atom.name));
        std::string atomType(array2string(atom.type));
        {
          size_t pos = atomName.find('\0');
          if(pos != std::string::npos)
            atomName.erase(pos);

          pos = atomType.find('\0');
          if(pos != std::string::npos)
            atomType.erase(pos);
        }

        stream << std::resetiosflags(std::ios::adjustfield);
        stream << std::setiosflags(std::ios::right);
        stream << "ATOM  " << std::setw(5) << atom.index << " ";
        stream << std::setw(2) << atomName;
        stream << std::resetiosflags(std::ios::adjustfield) << std::setw(2);
        stream << std::setiosflags(std::ios::left);
        stream << atomType << " " << std::setw(3);
        stream << aminoacidTriplet[residue.type] << " " << residue.chain;
        stream << std::resetiosflags(std::ios::adjustfield);
        stream << std::setiosflags(std::ios::right);
        stream << std::setw(4) << residue.index << "    ";
        stream << std::setiosflags(std::ios::fixed);
        stream << std::setprecision(3) << std::setw(8) << (atom.x * 10.);
        stream << std::setprecision(3) << std::setw(8) << (atom.y * 10.);
        stream << std::setprecision(3) << std::setw(8) << (atom.z * 10.);
        stream << std::setprecision(2) << std::setw(6);

        auto aux = getAuxParameters(atom);
        stream << ((std::get<1>(aux) >= 100) ? 99.99 : std::get<1>(aux));
        stream << std::setprecision(2) << std::setw(6);
        stream << ((std::get<0>(aux) >= 100) ? 99.99 : std::get<0>(aux));
        stream << std::resetiosflags(std::ios::fixed);
        stream << std::resetiosflags(std::ios::adjustfield);
        stream << std::setiosflags(std::ios::right);
        stream << "          " << std::setw(2) << atomName;
        stream << "  " <<  std::endl; /* Atom charge */
        stream << std::resetiosflags(std::ios::adjustfield);
      }
    }
  }
  stream.close();
}

template<typename AtomType>
template<typename Stream>
typename std::enable_if<is_stream_base_of<std::basic_istream, Stream>::value>::type
Pdb<AtomType>::readFromStream(Stream&& stream)
{
  std::stringstream streamLine;
  Protein<AtomType>* protein = nullptr;
  Residue<AtomType> residue;
  std::string buffer;
  int modelIndex = 0;
  char chain;
  Aminoacids residueType;

  residue.index = 0;

  while(not stream.eof())
  {
    getline(stream, buffer);
    streamLine.clear();
    streamLine.str(buffer);
    int residueIndex;

    streamLine.exceptions(std::ios_base::goodbit);
    streamLine >> std::setw(6) >> buffer;

    if(buffer == "ATOM")
    {
      AtomType atom;
      char atomType[5];

      if(protein == nullptr)
      {
        protein = new Protein<AtomType>();
        protein->model = ++modelIndex;
      }

      // Mandatory atom section
      streamLine.exceptions(std::ios_base::failbit);
      try
      {
        streamLine >> std::setw(5) >> atom.index;
        streamLine.seekg(1, std::ios_base::cur); // Empty space
        streamLine.read(atomType, 4);
        atomType[4] = '\0';
        atom.setAtomType(std::string(atomType));

        streamLine.seekg(1, std::ios_base::cur);// Alternate location indicator
        streamLine >> std::setw(3) >> buffer;
        residueType = Residue<AtomType>::getTypeByName(buffer);
        streamLine.seekg(1, std::ios_base::cur);// Empty space
        streamLine.read(&chain, 1);
        streamLine >> std::setw(4) >> residueIndex;
        {
          char insertionCode;
          streamLine.read(&insertionCode, 1);
          if(insertionCode != ' ')
            continue;
        }
        streamLine.seekg(3, std::ios_base::cur);// Code for insertion
                                           // of residues + 3 spaces
        streamLine >> std::setw(8) >> atom.x;
        streamLine >> std::setw(8) >> atom.y;
        streamLine >> std::setw(8) >> atom.z;

        atom.x /= 10.;
        atom.y /= 10.;
        atom.z /= 10.;
      }
      catch(std::ios_base::failure& fail)
      {
        continue;
      }

      // Non-mandatory atom section
      try
      {
        streamLine >> std::setw(6) >> atom.occupancy;
      }
      catch(std::ios_base::failure& fail)
      {
        atom.occupancy = 0;
      }

      try
      {
        streamLine >> std::setw(6) >> atom.bFactor;
      }
      catch(std::ios_base::failure& fail)
      {
        atom.bFactor = 0;
      }
      // Then element and charge... but I don't mind of them

      if(residue.index == 0)
      {
        residue.type = residueType;
        residue.index = residueIndex;
        residue.chain = chain;
      }
      else if(residueIndex != residue.index)
      {
        protein->appendResidue(std::move(residue));
        residue = Residue<AtomType>();
        residue.type = residueType;
        residue.index = residueIndex;
        residue.chain = chain;
      }

      residue.atoms.push_back(std::move(atom));
    }
    else if(buffer.substr(0,5) == "MODEL")
    {
      if(protein)
      {
        protein->appendResidue(std::move(residue));
        protein->lock();
        proteins.push_back(std::move(*protein));
      }
      protein = new Protein<AtomType>();
      residue = Residue<AtomType>();
      residue.index = 0;

      std::stringstream mdlIndex;
      {
        unsigned int nChars = 6;
        if(buffer.size() < nChars)
          nChars = buffer.size();
        mdlIndex << buffer.substr(nChars);
      }
      mdlIndex >> modelIndex;
      protein->model = modelIndex;
    }
    else if(protein and buffer == "ENDMDL")
    {
      protein->appendResidue(std::move(residue));
      protein->lock();
      proteins.push_back(std::move(*protein));
      protein = nullptr;
      residue = Residue<AtomType>();
      residue.index = 0;
    }
  }

  if(residue.atoms.size() > 0)
  {
    protein->appendResidue(std::move(residue));
  }

  if(protein and protein->residues().size() > 0)
  {
    protein->lock();
    proteins.push_back(std::move(*protein));
  }

  stream.close();
}

} /* namespace PstpFinder */

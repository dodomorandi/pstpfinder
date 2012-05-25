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
#include <cstring>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <locale>

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
    { "CYP", "CYS", "CYD", "CYS", "HID", "HIS", "HIE", "HIS", "LYP", "LYS",
      "LYD", "LYS" };
  const unsigned int aminoacidUncommonTranslatorSize =
      sizeof(aminoacidUncommonTranslator)
          / sizeof(*aminoacidUncommonTranslator);

  struct PdbAtom :
      public Atom
  {
      char type[5];
      unsigned int index;
      float bFactor;
      float occupancy;
      float sas;
  };

  struct Residue
  {
      Aminoacids type;
      int index;
      vector<PdbAtom> atoms;

      const PdbAtom&
      getAtomByType(string atomType) const
      {
        for(vector<PdbAtom>::const_iterator i = atoms.begin(); i < atoms.end(); i++)
          if(i->type == atomType)
            return *i;

        static PdbAtom unk;
        strcpy(unk.type, "UNK");

        return unk;
      }

      static Aminoacids
      getTypeByName(string residueName)
      {
        for(unsigned int i = 0; i < 21; i++)
        {
          for(unsigned int j = 0; j < 12; j += 2)
          {
            if(residueName == aminoacidUncommonTranslator[j])
            {
              residueName = aminoacidUncommonTranslator[j + 1];
              break;
            }
          }

          if(residueName == aminoacidTriplet[i])
            return static_cast<Aminoacids> (i);;
        }

        return AA_UNK;
      }
  };

  class Protein
  {
    public:
      string name;

      Protein()
      {
        locked = false;
      }

      Protein(const string& fileName)
      {
        ifstream pdbFile(fileName);
        // FIXME: I should move pdbFile, but I can't until streams
        //        will be moveable
        readFromStream(pdbFile);
      }

      template<typename Stream>
      Protein(typename enable_if<is_base_of<
                base_stream(basic_istream, Stream),
                Stream>::value,
              Stream>::type&& stream)
      {
        readFromStream(forward<Stream>(stream));
      }

      const vector<Residue>&
      residues() const
      {
        return pResidues;
      }

      vector<Residue>&
      residuesRW()
      {
        locked = false;
        return pResidues;
      }

      vector<const PdbAtom*>&
      atoms() const
      {
        if(not locked)
        {
          locked = true;
          pAtoms.clear();
          for(vector<Residue>::const_iterator i = pResidues.begin(); i
              < pResidues.end(); i++)
            for(vector<PdbAtom>::const_iterator j = i->atoms.begin(); j
                < i->atoms.end(); j++)
              pAtoms.push_back(&(*j));
        }

        return pAtoms;
      }

      bool
      appendResidue(Residue& residue)
      {
        if(locked)
          return false;

        pResidues.push_back(residue);
        return true;
      }

      void
      dumpPdb(const string& filename) const
      {
        ofstream pdb(filename, ios_base::out | ios_base::trunc);
        dumpPdb(pdb);
      }

      template<typename Stream>
      typename enable_if<is_base_of<base_stream(basic_ifstream, Stream),
                                    Stream>::value
                         or is_base_of<base_stream(basic_ofstream, Stream),
                                    Stream>::value>::type
      dumpPdb(Stream& stream) const
      {
        assert(stream.is_open());

        for(auto& residue : pResidues)
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
            stream << ((atom.occupancy >= 100) ? 99.99 : atom.occupancy);
            stream << setprecision(2) << setw(6);
            stream << ((atom.bFactor >= 100) ? 99.99 : atom.bFactor);
            stream << resetiosflags(ios::fixed);
            stream << "            " << setw(2) << atomType;
            stream << resetiosflags(ios::right) << endl;
          }
        }
        stream.close();
      }

      const PdbAtom&
      getAtomByIndex(unsigned int index) const
      {
        static PdbAtom unkAtom;
        unkAtom.index = -1;

        for(vector<const PdbAtom*>::const_iterator i = atoms().begin(); i
            < atoms().end(); i++)
        {
          if((*i)->index == index)
            return **i;
        }

        return unkAtom;
      }

      const Residue&
      getResidueByAtom(int atomIndex) const
      {
        return getResidueByAtom(getAtomByIndex(atomIndex));
      }

      const Residue&
      getResidueByAtom(const PdbAtom& atom) const
      {
        static Residue unkRes;
        unkRes.type = AA_UNK;

        for(vector<Residue>::const_iterator i = pResidues.begin(); i
            < pResidues.end(); i++)
          if(i->atoms[0].index > atom.index)
          {
            i--;
            for(vector<PdbAtom>::const_iterator j = i->atoms.begin(); j
                < i->atoms.end(); j++)
              if(j->index == atom.index)
                return *i;

            return unkRes;
          }

        return unkRes;
      }

      Protein&
      operator =(const Protein& protein)
      {
        pResidues = protein.pResidues;
        locked = false;
        name = protein.name;

        return *this;
      }

      const Residue& getResidueByIndex(int index) const
      {
        if(index >= 1
           and index - 1 < static_cast<int>(pResidues.size())
           and pResidues[index - 1].index == index)
          return pResidues[index - 1];
        else if(index >= 0
                and index < static_cast<int>(pResidues.size())
                and pResidues[index].index == index)
          return pResidues[index];
        else
        {
          for
          (
            vector<Residue>::const_iterator i = pResidues.begin();
            i < pResidues.end();
            i++
          )
          {
            if(i->index == index)
              return *i;
          }

          static Residue unk;
          unk.type = AA_UNK;
          return unk;
        }
      }

      void forceUnlock()
      {
        locked = false;
      }

    private:
      vector<Residue> pResidues;
      mutable vector<const PdbAtom*> pAtoms;
      mutable bool locked;

      template<typename Stream,
               typename Type = typename remove_reference<Stream>::type>
        typename enable_if<is_base_of<base_stream(basic_ifstream, Type),
                                      Type>::value>::type
      readFromStream(Stream&& stream)
      {
        stringstream streamLine;
        Residue residue;
        string buffer;
        Aminoacids residueType;

        residue.index = 0;
        streamLine.exceptions(ios_base::failbit);

        while(not stream.eof())
        {
          getline(stream, buffer);
          streamLine.clear();
          streamLine.str(buffer);
          PdbAtom atom;
          int residueIndex;

          // Mandatory atom section
          try
          {
            streamLine >> setw(6) >> buffer;
            if(buffer != "ATOM")
              continue;

            streamLine >> setw(5) >> atom.index;
            streamLine.seekg(1, ios_base::cur); // Empty space
            streamLine >> setw(5) >> atom.type;
            streamLine.seekg(1, ios_base::cur); // Alternate location indicator
            streamLine >> setw(3) >> buffer;
            residueType = Residue::getTypeByName(buffer);
            streamLine.seekg(2, ios_base::cur); // Empty space + Chain ID
            streamLine >> setw(4) >> residueIndex;
            streamLine.seekg(4, ios_base::cur); // Code for insertion
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
            pResidues.push_back(move(residue));
            residue = Residue();
            residue.type = residueType;
            residue.index = residueIndex;
          }

          residue.atoms.push_back(move(atom));
        }

        if(residue.atoms.size() > 0)
          pResidues.push_back(move(residue));

        locked = false;
        stream.close();
      }
  };
}

#endif

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

      Protein(const string& pdbFilename)
      {
        FILE* pdb;
        pdb = fopen(pdbFilename.c_str(), "r");
        Residue res;
        res.index = 0;

        locale oldLocale;
        locale::global(locale("C"));

        locked = false;
        while(not feof(pdb))
        {
          char resName[5], altLoc, chainId, iCode, element[2], charge[2];
          char line[256];
          char* linePtr;
          int resNdx;
          Aminoacids resType;
          PdbAtom atom;
          int ret;

          linePtr = fgets(line, 256, pdb);
          if(linePtr == 0)
            break;

          ret
              = sscanf(line, "ATOM  %5d %4s%c%3s %c%4d%c   "
                "%8f%8f%8f%6f%6f%*10c%2s%2s\n", &atom.index, atom.type,
                       &altLoc, resName, &chainId, &resNdx, &iCode, &atom.x,
                       &atom.y, &atom.z, &atom.occupancy, &atom.bFactor,
                       element, charge);

          if(ret <= 0)
            continue;

          atom.x /= 10.;
          atom.y /= 10.;
          atom.z /= 10.;

          resType = Residue::getTypeByName(string(resName));
          if(res.index == 0)
          {
            res.type = resType;
            res.index = resNdx;
          }
          else if(resNdx != res.index)
          {
            pResidues.push_back(res);
            res.atoms.clear();
            res.type = resType;
            res.index = resNdx;
          }

          res.atoms.push_back(atom);
        }
        pResidues.push_back(res);

        fclose(pdb);
        locale::global(oldLocale);
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
  };
}

#endif

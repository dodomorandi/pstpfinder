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
#include <locale>

using namespace std;

namespace Gromacs
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

  const string aminoacidLetter[] = { "X", "A", "V", "G", "L", "I", "T", "S",
                                     "M", "C", "P", "Y", "W", "F", "H", "R",
                                     "K", "E", "D", "Q", "N" };

  const string aminoacidTriplet[] = { "UNK", "ALA", "VAL", "GLY", "LEU", "ILE",
                                      "THR", "SER", "MET", "CYS", "PRO", "TYR",
                                      "TRP", "PHE", "HIS", "ARG", "LYS", "GLU",
                                      "ASP", "GLN", "ASN" };
  const string aminoacidUncommonTranslator[] =
               { "CYP", "CYS", "CYD", "CYS", "HID", "HIS", "HIE", "HIS", "LYP",
                 "LYS", "LYD", "LYS" };
  const unsigned int aminoacidUncommonTranslatorSize =
    sizeof(aminoacidUncommonTranslator) / sizeof(*aminoacidUncommonTranslator);

  struct PdbAtom: public Atom
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

    const PdbAtom& getAtomByType(string atomType) const
    {
      for
      (
        vector<PdbAtom>::const_iterator i = atoms.begin();
        i < atoms.end();
        i++
      )
        if(i->type == atomType)
          return *i;

      static PdbAtom unk;
      strcpy(unk.type, "UNK");

      return unk;
    }

    static Aminoacids getTypeByName(string residueName)
    {
      for(unsigned int i = 0; i < 21; i++)
      {
        for(unsigned int j = 0; j < 12; j += 2)
        {
          if(residueName == aminoacidUncommonTranslator[j])
          {
            residueName = aminoacidUncommonTranslator[j+1];
            break;
          }
        }

        if(residueName == aminoacidTriplet[i])
          return static_cast<Aminoacids>(i);;
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

        ret = sscanf(line, "ATOM  %5d %4s%c%3s %c%4d%c   "
                          "%8f%8f%8f%6f%6f%*10c%2s%2s\n",
                          &atom.index, atom.type, &altLoc, resName, &chainId,
                          &resNdx, &iCode, &atom.x, &atom.y, &atom.z,
                          &atom.occupancy, &atom.bFactor, element, charge);

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

    const vector<Residue>& residues() const
    {
      return pResidues;
    }

    vector<Residue>& residuesRW()
    {
      locked = false;
      return pResidues;
    }

    vector<const PdbAtom*>& atoms() const
    {
      if(not locked)
      {
        locked = true;
        pAtoms.clear();
        for
        (
          vector<Residue>::const_iterator i = pResidues.begin();
          i < pResidues.end();
          i++
        )
          for
          (
            vector<PdbAtom>::const_iterator j = i->atoms.begin();
            j < i->atoms.end();
            j++
          )
            pAtoms.push_back(&(*j));
      }

      return pAtoms;
    }

    bool appendResidue(Residue& residue)
    {
      if(locked)
        return false;

      pResidues.push_back(residue);
      return true;
    }

    void dumpPdb(const string& filename) const
    {
      FILE* pdb;
      pdb = fopen(filename.c_str(), "w");
      for
      (
        vector<Residue>::const_iterator i = pResidues.begin();
        i < pResidues.end();
        i++
      )
      {
        for
        (
          vector<PdbAtom>::const_iterator j = i->atoms.begin();
          j < i->atoms.end();
          j++
        )
        {
          fprintf(pdb, "ATOM  %5d %-4s %-3s A%4d    %8.3f%8.3f%8.3f%6.2f"
                       "%6.2f            %2s\n", j->index, j->type,
                       aminoacidTriplet[i->type].c_str(), i->index, j->x * 10.,
                       j->y * 10., j->z * 10.,
                       isinf(j->occupancy)? 99.99: j->occupancy,
                       isinf(j->bFactor)? 99.99: j->bFactor,
                       j->type);
        }
      }

      fclose(pdb);
    }

    const PdbAtom& getAtomByIndex(unsigned int index) const
    {
      static PdbAtom unkAtom;
      unkAtom.index = -1;

      for
      (
        vector<const PdbAtom*>::const_iterator i = atoms().begin();
        i < atoms().end();
        i++
      )
      {
        if((*i)->index == index)
          return **i;
      }

      return unkAtom;
    }


    const Residue& getResidueByAtom(int atomIndex) const
    {
      return getResidueByAtom(getAtomByIndex(atomIndex));
    }

    const Residue& getResidueByAtom(const PdbAtom& atom) const
    {
      static Residue unkRes;
      unkRes.type = AA_UNK;

      for
      (
        vector<Residue>::const_iterator i = pResidues.begin();
        i < pResidues.end();
        i++
      )
        if(i->atoms[0].index > atom.index)
        {
          i--;
          for
          (
            vector<PdbAtom>::const_iterator j = i->atoms.begin();
            j < i->atoms.end();
            j++
          )
            if(j->index == atom.index)
              return *i;

          return unkRes;
        }

      return unkRes;
    }
  private:
    vector<Residue> pResidues;
    mutable vector<const PdbAtom*> pAtoms;
    mutable bool locked;
  };
}

#endif

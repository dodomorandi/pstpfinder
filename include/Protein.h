#ifndef _PROTEIN_H
#define _PROTEIN_H

#include "Atom.h"
#include <cstring>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

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

      locked = false;
      while(not feof(pdb))
      {
        char resName[5];
        int resNdx;
        Aminoacids resType;
        PdbAtom atom;
        int ret;

        ret = fscanf(pdb, "ATOM%7d  %4s%3s A%4d %11f %7f %7f %5f %5f%*12c\n",
                          &atom.index, atom.type, resName, &resNdx, &atom.x,
                          &atom.y, &atom.z, &atom.occupancy, &atom.bFactor);


        if(ret != 9)
          abort();

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

    vector<PdbAtom*>& atoms()
    {
      if(not locked)
      {
        locked = true;
        pAtoms.clear();
        for
        (
          vector<Residue>::iterator i = pResidues.begin();
          i < pResidues.end();
          i++
        )
          for
          (
            vector<PdbAtom>::iterator j = i->atoms.begin();
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
          fprintf(pdb, "ATOM% 7d  %-4s%-4sA%4d %11.3f %7.3f %7.3f %5.2f"
                       " %5.2f%12c\n", j->index, j->type,
                       aminoacidTriplet[i->type].c_str(), i->index, j->x, j->y,
                       j->z, j->occupancy, /*j->bFactor*/ 0., j->type[0]);
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
        vector<PdbAtom*>::const_iterator i = pAtoms.begin();
        i < pAtoms.end();
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
    vector<PdbAtom*> pAtoms;
    bool locked;
  };
}

#endif

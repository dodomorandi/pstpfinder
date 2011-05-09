#ifndef _PROTEIN_H
#define _PROTEIN_H

#include "Atom.h"
#include <string>
#include <vector>

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
  };

  struct Residue
  {
    Aminoacids type;
    int index;
    vector<PdbAtom> atoms;
  };

  class Protein
  {
  public:
    string name;

    Protein()
    {
      locked = false;
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
  private:
    vector<Residue> pResidues;
    vector<PdbAtom*> pAtoms;
    bool locked;
  };
}

#endif

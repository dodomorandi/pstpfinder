#ifndef __PSTPPROTEIN_H
#define __PSTPPROTEIN_H

#include "Atom.h"
#include "string"

using namespace std;

namespace Gromacs
{
  enum Aminoacids
  {
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

  string aminoacidLetter[] = { "A", "V", "G", "L", "I", "T", "S", "M", "C", "P",
                               "Y", "W", "F", "H", "R", "K", "E", "D", "Q",
  "N" };

  string aminoacidTriplet[] = { "ALA", "VAL", "GLY", "LEU", "ILE", "THR", "SER",
                                "MET", "CYS", "PRO", "TYR", "TRP", "PHE", "HIS",
                                "ARG", "LYS", "GLU", "ASP", "GLN", "ASN" };

  struct PdbAtom: public Atom
  {
    unsigned int index;
    float bFactor;
    float occupancy;
  };

  struct Residue
  {
    unsigned int index;
    unsigned int atoms;
    PdbAtom* atom;
  };

  struct Protein
  {
    string name;
    vector<Residue> residue;
    unsigned int atoms;
    PdbAtom* atom;
  };
}

#endif

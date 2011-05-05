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

  const string aminoacidLetter[] = { "A", "V", "G", "L", "I", "T", "S", "M",
                                     "C", "P", "Y", "W", "F", "H", "R", "K",
                                     "E", "D", "Q", "N" };

  const string aminoacidTriplet[] = { "ALA", "VAL", "GLY", "LEU", "ILE", "THR",
                                      "SER", "MET", "CYS", "PRO", "TYR", "TRP",
                                      "PHE", "HIS", "ARG", "LYS", "GLU", "ASP",
                                      "GLN", "ASN" };

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
    unsigned int index;
    vector<PdbAtom> atoms;
  };

  class Protein
  {
  public:
    string name;
    vector<Residue> residues;
    vector<PdbAtom*> atoms;

    void appendResidue(const Residue& res)
    {
      residues.push_back(res);
      Residue& last = *(residues.end() - 1);
      for
      (
        vector<PdbAtom>::iterator i = last.atoms.begin();
        i < last.atoms.end();
        i++
      )
        atoms.push_back(&(*i));
    }
  };
}

#endif

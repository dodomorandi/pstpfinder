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

#include "Gromacs.h"
#include "Session.h"
#include "SasAnalysis.h"
#include "SasAtom.h"
#include "Pdb.h"

#include <iostream>
#include <string>
#include <cctype>
#include <cmath>
#include <utility>
#include <cstdio>

#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>

using namespace std;
using namespace boost::filesystem;

struct NeighSasPdbAtom : public PstpFinder::SasPdbAtom
{
    NeighSasPdbAtom() : PstpFinder::SasPdbAtom(), neighbours(0) {}

    unsigned long neighbours;
};

struct NeighSas
{
    NeighSas() : sas(0), neighbours(0) {}
    NeighSas(const NeighSasPdbAtom& atom) :
        sas(atom.sas), neighbours(atom.neighbours) {}

    float sas;
    unsigned long neighbours;
};

string getStdAtomNameFromPdb(string atomType)
{
  if(atomType[0] == ' ')
  {
    remove(begin(atomType), end(atomType), ' ');
    atomType.resize(1);
  }
  else
  {
    if(atomType.size() > 1 and islower(atomType[1]))
      atomType.resize(2);
    else if(atomType.size() > 1)
      atomType.resize(1);
  }

  return atomType;
}

template<typename AtomType>
auto
getBoundsFromNeighbours(vector<AtomType> atoms, const size_t& neighbours)
    noexcept -> decltype(make_pair(begin(atoms), end(atoms)))
{
  auto outBegin = begin(atoms);
  for(; outBegin != end(atoms); outBegin++)
  {
    if(outBegin->neighbours == neighbours)
      break;
  }

  auto outEnd = outBegin;
  for(; outEnd != end(atoms); outEnd++)
  {
    if(outEnd->neighbours != neighbours)
      break;
  }

  return make_pair(outBegin, outEnd);
}

void usage(string programName)
{
  cout << "Usage: " << programName << " <directory>" << endl;
  cout << "\tIt will recusively read every pdb inside 'directory' and make a";
  cout << endl << "\tcorrelation between atom neighbours and SAS" << endl;
}

int main(int argc, char* argv[])
{
  if(argc != 2)
  {
    usage(argv[0]);
    return -1;
  }

  path basedir(argv[1]);
  if(not exists(basedir) or not is_directory(basedir))
  {
    cerr << "Error: specified file doesn't exists or is not a directory" << endl;
    return -1;
  }

  list<path> pdbList;
  for(recursive_directory_iterator i(basedir);
      i != recursive_directory_iterator(); i++)
  {
    string file_ext = i->path().extension().string();
    transform(begin(file_ext), end(file_ext), begin(file_ext), ::tolower);
    if(file_ext == ".pdb")
      pdbList.push_back(i->path());
  }

  vector<NeighSas> atoms;
  double threshold = 0.286;

  for(const path& pdbPath : pdbList)
  {
    PstpFinder::Gromacs gromacs(pdbPath.string(), pdbPath.string());
    {
      PstpFinder::Session<ofstream> session("/tmp/tmp.csf", gromacs, 0.14, 0);
      gromacs.calculateSas(session);
      gromacs.waitOperation();
    }

    PstpFinder::Session<ifstream> session("/tmp/tmp.csf");
    PstpFinder::Pdb<NeighSasPdbAtom> pdb(pdbPath.string());
    PstpFinder::SasAtom* sasAtoms = 0;
    bool validThreshold = false;
    PstpFinder::SasAnalysis<ifstream> sasAnalysis(gromacs, session);
    sasAnalysis >> sasAtoms;
    for(unsigned int model = 0; sasAtoms != 0; model++, sasAnalysis >> sasAtoms)
    {
      Protein<NeighSasPdbAtom>& protein = pdb.proteins[model];
      protein.removeUnknownResidues();
      if(protein.residues().size() <= 30)
      {
        cout << "Skipping model " << model + 1 << " for pdb "
             << pdbPath.filename().string() << " (only "
             << protein.residues().size() << " residues)" << endl;
        continue;
      }

      unsigned long nAtoms = protein.atoms().size();
      auto iterAtoms = begin(protein.atoms());
      PstpFinder::SasAtom* atomPtr = sasAtoms;

      // FIXME: it crashes with curAtomIndex=4089, model=0 in pdb=3adc
      for(unsigned int curAtomIndex = 0; curAtomIndex < nAtoms;
          curAtomIndex++, iterAtoms++, atomPtr++)
        (const_cast<NeighSasPdbAtom&>(**iterAtoms)).sas = atomPtr->sas;

      for(; threshold < 1; threshold += 0.001)
      {
        double threshold2 = pow(threshold, 2.);
        for(auto& atom : protein.atoms())
          const_cast<NeighSasPdbAtom*>(atom)->neighbours = 0;

        iterAtoms = begin(protein.atoms());
        for(; iterAtoms != end(protein.atoms()); iterAtoms++)
        {
          NeighSasPdbAtom& atom1 = const_cast<NeighSasPdbAtom&>(**iterAtoms);
          string atomType1 = getStdAtomNameFromPdb(atom1.type);
          if(atomType1 == "H")
            continue;

          auto iterAtoms2 = begin(protein.atoms());
          for(; iterAtoms2 != end(protein.atoms()); iterAtoms2++)
          {
            if(iterAtoms == iterAtoms2)
              continue;

            NeighSasPdbAtom& atom2 = const_cast<NeighSasPdbAtom&>(**iterAtoms2);
            string atomType2 = getStdAtomNameFromPdb(atom2.type);
            if(atomType2 == "H")
              continue;

            //if(atom1.distance2(atom2) <= 0.32 * 0.32) /* ~hydrogen bond length² */
            if(atom1.distance2(atom2) <= threshold2)
              atom1.neighbours++;
          }
        }

        if(validThreshold)
          break;

        vector<NeighSas> tmpAtoms;
        tmpAtoms.reserve(protein.atoms().size());
        for(auto& residue : protein.residues())
        {
          if(not residue.isComplete())
            continue;

          for(auto& atom : residue.atoms)
          {
            if(getStdAtomNameFromPdb(atom.type) != "H")
              tmpAtoms.push_back(atom);
          }
        }

        if(tmpAtoms.size() == 0)
          break;

        sort(begin(tmpAtoms), end(tmpAtoms),
             [&](const NeighSas& a, const NeighSas& b)
             {
                return a.neighbours < b.neighbours;
             });

        unsigned int maxNeighbours = tmpAtoms.back().neighbours;
        {
          unsigned int minNeighbours = tmpAtoms.front().neighbours;
          if(minNeighbours == 0 or maxNeighbours < 6)
          {
            cout << "(" << pdbPath.filename().string()
                 << ") Not good for threshold = " << threshold
                 << " (min neighbours = " << minNeighbours << ", max = "
                 << maxNeighbours << "). Trying next..." << endl;
            continue;
          }
        }

        auto maxNeighBounds = getBoundsFromNeighbours(tmpAtoms, maxNeighbours);
        bool repeat = false;
        for(auto iter = maxNeighBounds.first;
            not repeat and iter < maxNeighBounds.second; iter++)
        {
          if(iter->sas > 0)
          {
            repeat = true;
            break;
          }
        }
        if(not repeat)
        {
          cout << "(" << pdbPath.filename().string()
               << ") Good for threshold = " << threshold << endl;
          validThreshold = true;
          break;
        }
        else
          cout << "(" << pdbPath.filename().string()
               << ") Not good for threshold = " << threshold
               << ". Trying next..." << endl;
      }
    }

    remove("/tmp/tmp.csf");
  }

  // Recalculate neighbours using best threshold
  double threshold2 = pow(threshold, 2.);
  for(const path& pdbPath : pdbList)
  {
    PstpFinder::Gromacs gromacs(pdbPath.string(), pdbPath.string());
    {
      PstpFinder::Session<ofstream> session("/tmp/tmp.csf", gromacs, 0.14, 0);
      gromacs.calculateSas(session);
      gromacs.waitOperation();
    }

    PstpFinder::Session<ifstream> session("/tmp/tmp.csf");
    PstpFinder::Pdb<NeighSasPdbAtom> pdb(pdbPath.string());
    PstpFinder::SasAtom* sasAtoms = 0;
    PstpFinder::SasAnalysis<ifstream> sasAnalysis(gromacs, session);
    sasAnalysis >> sasAtoms;

    for(unsigned int model = 0; sasAtoms != 0; model++, sasAnalysis >> sasAtoms)
    {
      auto& protein = pdb.proteins[model];
      protein.removeUnknownResidues();
      if(protein.residues().size() <= 30)
        continue;
      unsigned long nAtoms = protein.atoms().size();
      unsigned long validAtoms = 0;

      auto iterAtoms = begin(protein.atoms());
      PstpFinder::SasAtom* atomPtr = sasAtoms;
      for(unsigned int curAtomIndex = 0; curAtomIndex < nAtoms;
          curAtomIndex++, iterAtoms++, atomPtr++)
        (const_cast<NeighSasPdbAtom&>(**iterAtoms)).sas = atomPtr->sas;

      iterAtoms = begin(protein.atoms());
      for(; iterAtoms != end(protein.atoms()); iterAtoms++)
      {
        NeighSasPdbAtom& atom1 = const_cast<NeighSasPdbAtom&>(**iterAtoms);
        string atomType1 = getStdAtomNameFromPdb(atom1.type);
        if(atomType1 == "H")
          continue;

        validAtoms++;

        auto iterAtoms2 = begin(protein.atoms());
        for(; iterAtoms2 != end(protein.atoms()); iterAtoms2++)
        {
          if(iterAtoms == iterAtoms2)
            continue;

          NeighSasPdbAtom& atom2 = const_cast<NeighSasPdbAtom&>(**iterAtoms2);
          string atomType2 = getStdAtomNameFromPdb(atom2.type);
          if(atomType2 == "H")
            continue;

          if(atom1.distance2(atom2) <= threshold2)
            atom1.neighbours++;
        }
      }

      atoms.reserve(atoms.size() + validAtoms);
      for(auto atom : protein.atoms())
      {
        if(getStdAtomNameFromPdb(atom->type) != "H")
          atoms.push_back(*atom);
      }
    }

    remove("/tmp/tmp.csf");
  }

  sort(begin(atoms), end(atoms),
       [&](const NeighSas& a, const NeighSas& b)
       {
          return a.neighbours < b.neighbours;
       });

  auto neighBounds = make_pair(atoms.front().neighbours,
                                atoms.back().neighbours + 1);
  for(unsigned int curNeighbours = neighBounds.first;
      curNeighbours < neighBounds.second; curNeighbours++)
  {
    auto atomBounds = getBoundsFromNeighbours(atoms, curNeighbours);
    long nAtoms = distance(atomBounds.first, atomBounds.second);
    if(nAtoms == 0)
      continue;
    float averageSas = 0;
    for(auto iter = atomBounds.first; iter != atomBounds.second; iter++)
      averageSas += iter->sas;
    averageSas /= nAtoms;

    cout << "Average SAS for " << curNeighbours << " neighbours (" << nAtoms
         << ") = " << averageSas << endl;
  }

  return 0;
}

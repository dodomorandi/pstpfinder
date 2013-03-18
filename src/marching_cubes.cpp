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

#include <iostream>
#include <algorithm>
#include <cmath>
#include <functional>
#include <set>

#include "marching_cubes.h"
#include "SplittedSpace.h"
#include "SpaceCube.h"
#include "SasAnalysis.h"
#include "Gromacs.h"
#include "Session.h"
#include "SasAtom.h"
#include "AtomInfo.h"

using namespace std;

int main(int argc, char* argv[])
{
  if(argc != 2)
  {
    cout << "Usage: " << argv[0] << " session.csf" << endl;
    return -1;
  }

  PstpFinder::Session<ifstream> session(argv[1]);
  PstpFinder::Gromacs gromacs(session.getTrajectoryFileName(),
                              session.getTopologyFileName(),
                              session.getRadius());

  /*
   * Gromacs class is taken as const inside sasAnalysis.
   * We need to load files BEFORE passing it to other objects
   */
  gromacs.loadTrajectoryAndTopology();
  PstpFinder::SasAnalysis<ifstream> sasAnalysis(gromacs, session);

  const matrix& box = gromacs.getBox();
  if(box[XX][YY] != 0 or box[XX][ZZ] != 0 or box[YY][XX] != 0
     or box[YY][ZZ] != 0 or box[ZZ][XX] != 0 or box[ZZ][YY] != 0)
  {
    cerr << "Sorry, only rectangular boxes supported for now." << endl;
    return -2;
  }

  PstpFinder::SplittedSpace space(box[XX][XX], box[YY][YY], box[ZZ][ZZ],
                                  PstpFinder::cubesize);

  unsigned long atomsSize = gromacs.getGroup("Protein").size();
  vector<PstpFinder::AtomInfo> atomsInfo = gromacs.getGroupInfo();
  // Add water radius to atoms radius
  for(PstpFinder::AtomInfo& info : atomsInfo)
    info.radius += 0.14;

  unsigned long frameIndex = 0;
  PstpFinder::SasAtom* sasAtom;
  sasAnalysis >> sasAtom;
  while(sasAtom != nullptr)
  {
    unsigned long atomIndex = 0;
    PstpFinder::SasAtom* atomPtr = sasAtom;
    for(; atomIndex < atomsSize; atomIndex++, atomPtr++)
    {
      if(atomsInfo[atomIndex].element == 'H')
        continue;

      real radius = atomsInfo[atomIndex].radius;
      real radius2 = radius * radius;

      list<PstpFinder::SpaceCube*> involvedCubes = space.getInvolvedCubes(
          atomPtr->x, atomPtr->y, atomPtr->z, radius);
      unordered_map<array<unsigned, 3>, float> distance2Map;
      list<PstpFinder::SpaceCube*> cubesNotInvolved;

      for(PstpFinder::SpaceCube* cube : involvedCubes)
      {
        unsigned vertexIndex = 0;
        cube->involvedAtomsIndices.push_back(atomIndex);
        array<unsigned, 3> cubeIndex = space.cubeIndexAt(cube->x(), cube->y(),
                                                         cube->z());

        bool cubeInvolved = false;
        for(const array<unsigned, 3>& vertex : PstpFinder::vertices)
        {
          array<unsigned, 3> index = cubeIndex;
          float distance2;
          {
            auto vertexIter = begin(vertex);
            auto indexIter = begin(index);
            for(; indexIter != end(index); indexIter++, vertexIter++)
              *indexIter += *vertexIter;
          }

          array<float, 3> point =
            {
              { PstpFinder::cubesize * index[0], PstpFinder::cubesize
                  * index[1],
                PstpFinder::cubesize * index[2] } };
          distance2 = pow(atomPtr->x - point[0], 2)
              + pow(atomPtr->y - point[1], 2) + pow(atomPtr->z - point[2], 2);

          if(distance2 <= radius2)
          {
            cube->flags.set(vertexIndex, true);
            cubeInvolved = true;
          }

          vertexIndex++;
        }

        if(not cubeInvolved)
          cubesNotInvolved.push_back(cube);
      }

      for(PstpFinder::SpaceCube* cubeNotInvolved : cubesNotInvolved)
      {
        cubeNotInvolved->involvedAtomsIndices.remove(atomIndex);
        involvedCubes.remove(cubeNotInvolved);
      }
    }

    set<unsigned long> involvedAtoms;
    for(PstpFinder::SpaceCube& cube : space.getCubes())
    {
      bitset<12> flags = PstpFinder::cubeEdgeFlags[cube.flags.to_ulong()];
      if(flags.any() and not flags.all())
      {
        for(unsigned int index : cube.involvedAtomsIndices)
          involvedAtoms.insert(index);
      }
    }

    cout << "cmd.select(\"frame_" << frameIndex + 1 << "\", \"ID ";
    for(unsigned long index : involvedAtoms)
      cout << index + 1 << "+";
    cout << "\", state=" << frameIndex + 1 << ")" << endl;

    space.clearFlags();
    frameIndex++;
    sasAnalysis >> sasAtom;
  }
}

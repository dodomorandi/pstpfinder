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

#include "MarchingCubes.h"
#include "SplittedSpace.h"
#include "SpaceCube.h"
#include "SasAnalysis.h"
#include "Gromacs.h"
#include "Session.h"
#include "SasAtom.h"
#include "AtomInfo.h"
#include "Viewer3D.h"

using namespace std;
using namespace PstpFinder::MarchingCubes;

constexpr static float cubesize = 0.12;

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

  vector<PstpFinder::AtomInfo> atomsInfo = gromacs.getGroupInfo();
  // Add water radius to atoms radius
  /*
  for(PstpFinder::AtomInfo& info : atomsInfo)
    info.radius += 0.14;
  */

  Viewer3D viewer;
  MarchingCubes marchingCubes(cubesize, array<float, 3>
    {{ box[XX][XX], box[YY][YY], box[ZZ][ZZ] }});

  unsigned long frameIndex = 0;
  for(const vector<PstpFinder::SasAtom>& sasAtoms : sasAnalysis)
  {
    unsigned atomIndex = 0;
    vector<PointRadius> points;
    for(const PstpFinder::SasAtom& atom : sasAtoms)
    {
      PstpFinder::AtomInfo info = atomsInfo[atomIndex++];
      if(info.element == 'H')
        continue;

      points.emplace_back(atom, info.radius);
    }
    viewer.setPointsRadii(points);
    const SplittedSpace& space = marchingCubes.run(points);

    set<unsigned long> involvedAtoms;
    for(const SpaceCube& cube : space.getCubes())
    {
      bitset<12> flags = DataSet::cubeEdgeFlags[cube.flags.to_ulong()];
      if(flags.any() and not flags.all())
      {
        for(unsigned int index : cube.involvedPointsIndices)
          involvedAtoms.insert(index);
      }
    }

    viewer.setSpace(space);

    frameIndex++;
    break;
  }

  viewer.run();
}

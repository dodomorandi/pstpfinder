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

#include "MarchingCubes.h"

#include <cmath>

namespace PstpFinder { namespace MarchingCubes {

MarchingCubes::MarchingCubes(float cubeEdge, array<float, 3> boxSize) noexcept :
  cubeEdge(cubeEdge),
  boxSize(boxSize)
{
  space = unique_ptr<SplittedSpace>
    { new SplittedSpace(boxSize[0], boxSize[1], boxSize[2], cubeEdge)};
}

const SplittedSpace&
MarchingCubes::run(const vector<PointRadius>& points)
{
  space->clear();

  unsigned long pointIndex = 0;
  for(const PointRadius& point : points)
  {
    float radius2 = pow(point.radius, 2);

    list<SpaceCube*> involvedCubes = space->getInvolvedCubes(
        point.x, point.y, point.z, point.radius);
    list<SpaceCube*> cubesNotInvolved;

    for(SpaceCube* cube : involvedCubes)
    {
      unsigned vertexIndex = 0;
      cube->involvedPointsIndices.push_back(pointIndex);
      array<unsigned, 3> cubeIndex = space->cubeIndexAt(cube->x(), cube->y(),
                                                       cube->z());

      bool cubeInvolved = false;
      for(const array<unsigned, 3>& vertex : DataSet::vertices)
      {
        array<unsigned, 3> index = cubeIndex;
        float distance2;
        {
          auto vertexIter = begin(vertex);
          auto indexIter = begin(index);
          for(; indexIter != end(index); indexIter++, vertexIter++)
            *indexIter += *vertexIter;
        }

        array<float, 3> coord = space->getPointCoordinatesAtIndex(index);
        distance2 = pow(point.x - coord[0], 2)
            + pow(point.y - coord[1], 2) + pow(point.z - coord[2], 2);
        if(distance2 <= radius2)
        {
          if(not cube->flags.test(vertexIndex) or distance2
              < cube->verticesDistance[vertexIndex].second)
          {
            cube->verticesDistance[vertexIndex].first = pointIndex;
            cube->verticesDistance[vertexIndex].second = distance2;
          }

          cube->flags.set(vertexIndex, true);
          cubeInvolved = true;
        }

        vertexIndex++;
      }

      if(not cubeInvolved)
        cubesNotInvolved.push_back(cube);
    }

    for(SpaceCube* cubeNotInvolved : cubesNotInvolved)
    {
      cubeNotInvolved->involvedPointsIndices.remove(pointIndex);
      involvedCubes.remove(cubeNotInvolved);
    }

    pointIndex++;
  }

  return *space;
}

} /* namespace MarchingCubes */ } /* namespace PstpFinder */

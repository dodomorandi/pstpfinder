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
#include <limits>

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

  {
    float minRadius = points[0].radius;
    float maxRadius = minRadius;
    for(const PointRadius& point : points)
    {
      if(point.radius < minRadius)
        minRadius = point.radius;
      if(point.radius > maxRadius)
        maxRadius = point.radius;
    }

    bool changeSpace = true;
    if(cubeEdge < minRadius * 0.4)
      cubeEdge = minRadius * 0.4;
    else if(cubeEdge > maxRadius * 0.7)
      cubeEdge = maxRadius * 0.7;
    else
      changeSpace = false;

    if(changeSpace)
        space = unique_ptr<SplittedSpace>(
            new SplittedSpace(boxSize[0], boxSize[1], boxSize[2], cubeEdge));
  }

  unsigned long pointIndex = 0;
  for(const PointRadius& point : points)
  {
    float radius2 = pow(point.radius, 2);

    list<SpaceCube*> involvedCubes = space->getInvolvedCubes(
        point.x, point.y, point.z, point.radius);
    vector<SpaceCube*> cubesNotInvolved;

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
          if(not cube->flags.test(vertexIndex) or radius2 - distance2
              < points[cube->verticesDistance[vertexIndex].first].radius -
                cube->verticesDistance[vertexIndex].second)
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

vector<PointAndNormal>
MarchingCubes::getMesh(const vector<PointRadius>& points) const
{
  if(not space)
    throw;

  vector<PointAndNormal> pointsAndNormals;
  for(const SpaceCube& cube : space->getCubes())
  {
    const bitset<12>& flags = DataSet::cubeEdgeFlags[cube.flags.to_ulong()];
    const array<int, 15>& triangleEdgesIndices =
        DataSet::triangleConnectionTable[cube.flags.to_ulong()];
    array<tuple<bool, unsigned, float>, 8> distanceCache;
    for(auto cache : distanceCache)
      get<0>(cache) = false;

    if(not flags.any() or flags.all())
      continue;

    for(unsigned flagIndex = 0; flagIndex < 12; flagIndex++)
    {
      if(not flags[flagIndex])
        continue;

      for(int edgeIndex : triangleEdgesIndices)
      {
        if(edgeIndex == -1)
          break;

        const array<unsigned, 2>& vertexIndices =
            DataSet::edgesVertices[edgeIndex];
        unsigned pointIndex;
        array<float, 2> distances2;
        array<float, 2> distances;
        array<array<float, 3>, 2> positions;
        PointRadius point;

        for(unsigned i = 0; i < 2; i++)
        {
          positions[i][0] = cube.x()
              + space->cubeEdgeSize * DataSet::vertices[vertexIndices[i]][0];
          positions[i][1] = cube.y()
              + space->cubeEdgeSize * DataSet::vertices[vertexIndices[i]][1];
          positions[i][2] = cube.z()
              + space->cubeEdgeSize * DataSet::vertices[vertexIndices[i]][2];
        }

        if(cube.flags[vertexIndices[0]])
        {
          tie(pointIndex, distances2[0]) =
              cube.verticesDistance[vertexIndices[0]];
          point = points[pointIndex];
          distances2[1] = pow(point.x - positions[1][0], 2)
                  + pow(point.y - positions[1][1], 2)
                  + pow(point.z - positions[1][2], 2);
        }
        else
        {
          tie(pointIndex, distances2[1]) =
              cube.verticesDistance[vertexIndices[1]];
          point = points[pointIndex];
          distances2[0] = pow(point.x - positions[0][0], 2)
              + pow(point.y - positions[0][1], 2)
              + pow(point.z - positions[0][2], 2);
        }

        float offset;
        distances = {{ sqrt(distances2[0]), sqrt(distances2[1]) }};
        float delta = distances[1] - distances[0];
        offset = (point.radius - distances[0])/delta;
        if(abs(offset) == numeric_limits<float>::infinity())
          offset = 0.5;

        PointAndNormal pointAndNormal;
        pointAndNormal.x() = positions[0][0]
            + cubeEdge * offset * DataSet::edgesDirections[edgeIndex][0];
        pointAndNormal.y() = positions[0][1]
            + cubeEdge * offset * DataSet::edgesDirections[edgeIndex][1];
        pointAndNormal.z() = positions[0][2]
            + cubeEdge * offset * DataSet::edgesDirections[edgeIndex][2];

        pointAndNormal.nx() = pointAndNormal.x() - point.x;
        pointAndNormal.ny() = pointAndNormal.y() - point.y;
        pointAndNormal.nz() = pointAndNormal.z() - point.z;

        float normFactor = sqrt(
            pow(pointAndNormal.nx(), 2) + pow(pointAndNormal.ny(), 2)
            + pow(pointAndNormal.nz(), 2));

        if(normFactor != 0)
        {
          pointAndNormal.nx() /= normFactor;
          pointAndNormal.ny() /= normFactor;
          pointAndNormal.nz() /= normFactor;
        }

        pointsAndNormals.push_back(pointAndNormal);
      }
    }
  }

  return pointsAndNormals;
}

} /* namespace MarchingCubes */ } /* namespace PstpFinder */

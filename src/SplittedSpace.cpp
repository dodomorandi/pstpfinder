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

#include "SplittedSpace.h"

using namespace std;
using namespace PstpFinder::MarchingCubes;

SplittedSpace::SplittedSpace(float sizeX, float sizeY, float sizeZ,
                             float cubeEdge) :
    cubeEdgeSize(cubeEdge),
    size({{sizeX, sizeY, sizeZ}}),
    nCubes(getNCubes(sizeX, sizeY, sizeZ, cubeEdgeSize))
{
  spaceCubes.reserve(nCubes[0] * nCubes[1] * nCubes[2]);
  for(unsigned z = 0; z < nCubes[2]; z++)
  {
    for(unsigned y = 0; y < nCubes[1]; y++)
    {
      for(unsigned x = 0; x < nCubes[0]; x++)
        spaceCubes.emplace_back(cubeEdgeSize * x, cubeEdgeSize * y, cubeEdgeSize * z);
    }
  }
}

SplittedSpace::SplittedSpace(float sizeX, float sizeY, float sizeZ,
                             unsigned cubesPerDimension) :
    SplittedSpace(sizeX, sizeY, sizeZ,
                  getCubeSize(sizeX, sizeY, sizeZ, cubesPerDimension))
{
}

float
SplittedSpace::getCubeSize(float sizeX, float sizeY, float sizeZ,
                         unsigned cubesPerDimension)
{
  float size = sizeX;
  if(sizeY > size)
    size = sizeY;
  if(sizeZ > size)
    size = sizeZ;

  return size / cubesPerDimension;
}

array<unsigned, 3>
SplittedSpace::getNCubes(float sizeX, float sizeY, float sizeZ, float cubeEdgeSize)
{
  array<unsigned, 3> nCubes;

  nCubes[0] = sizeX / cubeEdgeSize;
  nCubes[1] = sizeY / cubeEdgeSize;
  nCubes[2] = sizeZ / cubeEdgeSize;

  if(cubeEdgeSize * nCubes[0] < sizeX)
    nCubes[0]++;
  if(cubeEdgeSize * nCubes[1] < sizeY)
    nCubes[1]++;
  if(cubeEdgeSize * nCubes[2] < sizeZ)
    nCubes[2]++;

  return nCubes;
}

inline SpaceCube&
SplittedSpace::operator()(unsigned xIndex, unsigned yIndex, unsigned zIndex)
{
  return spaceCubes[xIndex + yIndex * nCubes[0] + zIndex * nCubes[1] * nCubes[0]];
}

const SpaceCube&
SplittedSpace::cubeAt(float x, float y, float z)
{
  return operator()(x / cubeEdgeSize, y / cubeEdgeSize,
                    z / cubeEdgeSize);
}

array<unsigned, 3>
SplittedSpace::cubeIndexAt(float x, float y, float z)
{
  return array<unsigned, 3>
    {{ static_cast<unsigned>(x / cubeEdgeSize),
       static_cast<unsigned>(y / cubeEdgeSize),
       static_cast<unsigned>(z / cubeEdgeSize) }};
}

list<SpaceCube*>
SplittedSpace::getInvolvedCubes(float x, float y, float z, float radius)
{
  list<SpaceCube*> cubesList;
  array<unsigned, 3> minIndices, maxIndices, centralIndex;

  unsigned discreteRadius = radius / cubeEdgeSize;
  if(cubeEdgeSize * discreteRadius < radius)
    discreteRadius++;

  centralIndex = cubeIndexAt(x, y, z);
  minIndices = centralIndex;
  for(unsigned& index : minIndices)
  {
    if(index >= discreteRadius)
      index -= discreteRadius;
    else
      index = 0;
  }
  maxIndices = centralIndex;
  // FIXME: use pointers or iterators for optimized code
  for(unsigned short index = 0; index < 3; index++)
  {
    unsigned tmpIndex = maxIndices[index] + discreteRadius;
    if(tmpIndex >= nCubes[index])
        tmpIndex = nCubes[index] - 1;

    maxIndices[index] = tmpIndex;
  }

  for(unsigned x = minIndices[0]; x <= maxIndices[0]; x++)
  {
    for(unsigned y = minIndices[1]; y <= maxIndices[1]; y++)
    {
      for(unsigned z = minIndices[2]; z <= maxIndices[2]; z++)
        cubesList.push_back(&operator()(x, y, z));
    }
  }

  return cubesList;
}

vector<SpaceCube>&
SplittedSpace::getCubes() noexcept
{
  return spaceCubes;
}

const vector<SpaceCube>&
SplittedSpace::getCubes() const noexcept
{
  return spaceCubes;
}

void
SplittedSpace::clear() noexcept
{
  for(SpaceCube& cube : spaceCubes)
  {
    cube.flags.reset();
    cube.involvedPointsIndices.clear();
  }
}

array<float, 3>
SplittedSpace::getPointCoordinatesAtIndex(const array<unsigned, 3>& index) const
{
  array<float, 3> point;
  point[0] = cubeEdgeSize * index[0];
  point[1] = cubeEdgeSize * index[1];
  point[2] = cubeEdgeSize * index[2];
  return point;
}

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

#ifndef SPLITTEDSPACE_H_
#define SPLITTEDSPACE_H_

#include <array>
#include <vector>
#include <list>

#include "SpaceCube.h"

using namespace std;

namespace PstpFinder
{
  namespace MarchingCubes
  {
    class SplittedSpace
    {
      public:
        SplittedSpace(float sizeX, float sizeY, float sizeZ, float cubeEdge);
        SplittedSpace(float sizeX, float sizeY, float sizeZ,
                      unsigned cubesPerDimension);
        inline SpaceCube&
          operator()(unsigned xIndex, unsigned yIndex, unsigned zIndex);
        const SpaceCube& cubeAt(float x, float y, float z);
        array<unsigned, 3>
          cubeIndexAt(float x, float y, float z);
        list<SpaceCube*>
          getInvolvedCubes(float x, float y, float z, float radius);
        vector<SpaceCube>& getCubes() noexcept;
        const vector<SpaceCube>& getCubes() const noexcept;
        void clear() noexcept;
        /* Can also be used for OOB indices. However be careful! */
        array<float, 3>
          getPointCoordinatesAtIndex(const array<unsigned, 3>& index) const;

        const float cubeEdgeSize;

      private:
        const array<float, 3> size;
        const array<unsigned, 3> nCubes;
        vector<SpaceCube> spaceCubes;

        static float
          getCubeSize(float sizeX, float sizeY, float sizeZ,
                    unsigned cubesPerDimension);
        static array<unsigned, 3>
          getNCubes(float sizeX, float sizeY, float sizeZ, float cubeSize);
    };
  } /* namespace MarchingCubes */
} /* namespace PstpFinder */
#endif /* SPLITTEDSPACE_H_ */

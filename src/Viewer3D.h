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

#ifndef VIEWER3D_H_
#define VIEWER3D_H_

#include "SplittedSpace.h"
#include "SasAtom.h"
#include "AtomInfo.h"
#include "Atom.h"
#include "MarchingCubes.h"
#include "PointRadius.h"
#include "Vector3d.h"

#include <array>
#include <vector>

#include <GL/gl.h>

namespace PstpFinder
{
  namespace MarchingCubes
  {
    class Viewer3D
    {
      public:
        Viewer3D();
        ~Viewer3D();
        void run();
        void appendFrame(const vector<PointAndNormal>& frame);

      private:
        int window;
        vector<vector<PointAndNormal>> frames;
        array<int, 2> lastMousePosition;
        array<float, 2> rotationXY;
        array<Vector3d, 2> boundaries;
        bool boundariesSet;

        void init();
        void display();
    };
  } /* namespace MarchingCubes */
} /* namespace PstpFinder */
#endif /* VIEWER3D_H_ */

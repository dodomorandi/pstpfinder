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

#include <array>
#include <list>
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
        void setPointsRadii(const vector<PointRadius>& points);
        void setSpace(const SplittedSpace& space);

      private:
        int window;
        vector<PointRadius> points;
        // A list of frame. Every frame has got lots of vertex info. Every vertex
        // info is compound of two triplet of floats. The first triplet represents
        // the coordinates, the other the normal vector. Each triplet is simply
        // x, y and z.
        list<list<array<array<GLfloat, 3>, 2>>> frames;
        array<int, 2> lastMousePosition;
        array<float, 2> rotationXZ;
        array<array<float, 3>, 2> boundaries;
        bool boundariesSet;

        void init();
        void display();
    };
  } /* namespace MarchingCubes */
} /* namespace PstpFinder */
#endif /* VIEWER3D_H_ */

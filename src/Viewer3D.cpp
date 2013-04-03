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

#include "Viewer3D.h"
#include "marching_cubes.h"

#include <iostream>
#include <tuple>

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glfw.h>

using namespace std;
using namespace PstpFinder;

Viewer3D::Viewer3D(const vector<PstpFinder::AtomInfo>& atomsInfo) :
    boundariesSet(false)
{
  this->atomsInfo = atomsInfo;
  glfwInit();
}

Viewer3D::~Viewer3D()
{
  glfwTerminate();
}

void
Viewer3D::setAtoms(const SasAtom* atoms)
{
  unsigned nAtoms = atomsInfo.size();
  this->atoms.clear();

  const SasAtom* atomPtr = atoms;
  if(not boundariesSet)
  {
    boundaries = {{ {{ atomPtr->x, atomPtr->y, atomPtr->z }},
                    {{ atomPtr->x, atomPtr->y, atomPtr->z }} }};
    boundariesSet = true;
  }
  for(unsigned atomIndex = 0; atomIndex < nAtoms; atomIndex++, atomPtr++)
  {
    this->atoms.push_back(*atomPtr);
    if(atomPtr->x < boundaries[0][0])
      boundaries[0][0] = atomPtr->x;
    if(atomPtr->x > boundaries[1][0])
      boundaries[1][0] = atomPtr->x;
    if(atomPtr->y < boundaries[0][1])
      boundaries[0][1] = atomPtr->y;
    if(atomPtr->y > boundaries[1][1])
      boundaries[1][1] = atomPtr->y;
    if(atomPtr->z < boundaries[0][2])
      boundaries[0][2] = atomPtr->z;
    if(atomPtr->z > boundaries[1][2])
      boundaries[1][2] = atomPtr->z;
  }
}

void
Viewer3D::setSpace(const SplittedSpace& space)
{
  list<array<array<GLfloat, 3>, 2>> points;
  for(const SpaceCube& cube : space.getCubes())
  {
    const bitset<12>& flags = cubeEdgeFlags[cube.flags.to_ulong()];
    const array<int, 15>& triangleEdgesIndices = triangleConnectionTable[cube.flags
        .to_ulong()];
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

        const array<unsigned, 2>& vertexIndices = edgesVertices[edgeIndex];
        unsigned atomIndex;
        array<float, 2> distances2;
        array<float, 2> distances;
        array<array<float, 3>, 2> positions;
        Atom atom;

        for(unsigned i = 0; i < 2; i++)
        {
          positions[i][0] = cube.x() + cubesize * vertices[vertexIndices[i]][0];
          positions[i][1] = cube.y() + cubesize * vertices[vertexIndices[i]][1];
          positions[i][2] = cube.z() + cubesize * vertices[vertexIndices[i]][2];
        }

        if(cube.flags[vertexIndices[0]])
        {
          tie(atomIndex, distances2[0]) =
              cube.verticesDistance[vertexIndices[0]];
          atom = atoms[atomIndex];
          distances2[1] = pow(atom.x - positions[1][0], 2)
                  + pow(atom.y - positions[1][1], 2)
                  + pow(atom.z - positions[1][2], 2);
        }
        else
        {
          tie(atomIndex, distances2[1]) =
              cube.verticesDistance[vertexIndices[1]];
          atom = atoms[atomIndex];
          distances2[0] = pow(atom.x - positions[0][0], 2)
              + pow(atom.y - positions[0][1], 2)
              + pow(atom.z - positions[0][2], 2);
        }

        float offset;
        distances = {{ sqrt(distances[0]), sqrt(distances[1]) }};
        float delta = distances[1] - distances[0];
        if(delta < 0.00001)
          offset = cubesize / 2;
        else
          offset = (atomsInfo[atomIndex].radius - distances[0])/delta;

        array<array<GLfloat, 3>, 2> pointNNormal;
        pointNNormal[0] = array<GLfloat, 3>
          {{ positions[0][0] + offset * edgesDirections[edgeIndex][0],
             positions[0][1] + offset * edgesDirections[edgeIndex][1],
             positions[0][2] + offset * edgesDirections[edgeIndex][2] }};
        pointNNormal[1] = array<GLfloat, 3>
          {{ pointNNormal[0][0] - atom.x, pointNNormal[0][1] - atom.y,
             pointNNormal[0][2] - atom.z }};

        float normFactor = sqrt(
            pow(pointNNormal[1][0], 2) + pow(pointNNormal[1][1], 2)
            + pow(pointNNormal[1][2], 2));

        if(normFactor != 0)
          pointNNormal[1] = array<GLfloat, 3>
            {{ pointNNormal[1][0] / normFactor, pointNNormal[1][1] / normFactor,
                pointNNormal[1][2] / normFactor }};

        points.push_back(pointNNormal);
      }
    }
  }

  frames.push_back(points);
}

void
Viewer3D::run()
{
  glfwOpenWindow(800,600,0,0,0,0,1,0,GLFW_WINDOW);
  init();

  while (glfwGetWindowParam(GLFW_OPENED))
  {
    display();
    glfwSwapBuffers();
  }
}

void
Viewer3D::init()
{
  rotationXZ = {{ 0, 0 }};
  glClearColor(0, 0, 0, 0);
  glShadeModel(GL_SMOOTH);
  glClearDepth(1);

  static GLfloat light_diffuse[] = {1.0, 1.0, 1.0, 1.0};
  static GLfloat light_position[] =
  { (boundaries[1][0] + boundaries[0][0]) / 2,
    (boundaries[1][1] + boundaries[0][1]) / 2 + 2.f,
    boundaries[1][2] + 2.f, 0.0};
  glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
  glLightfv(GL_LIGHT0, GL_POSITION, light_position);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);

  glEnable(GL_DEPTH_TEST);

  int width, height;
  glfwGetWindowSize(&width, &height);
  glViewport(0, 0, width, height);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(60, static_cast<float>(width) / static_cast<float>(height),
                 1, 20);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(0, 0, boundaries[1][2] + 1.f, 0, 0, 0, 0, 1, 0);
}

void
Viewer3D::display()
{
  if(frames.size() == 0)
    return;

  bool mouseClick = (glfwGetMouseButton(GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS);
  array<int, 2> mousePosition;

  if(mouseClick)
  {
    glfwGetMousePos(&mousePosition[0], &mousePosition[1]);
    rotationXZ[0] += lastMousePosition[1] - mousePosition[1];
    rotationXZ[1] += mousePosition[0] - lastMousePosition[0];
  }

  static auto frameIter = begin(frames);

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glColor3f(1, 1, 1);

  glPushMatrix();
  glRotatef(rotationXZ[0], 1, 0, 0);
  glRotatef(rotationXZ[1], 0, 0, 1);
  glTranslatef(-(boundaries[1][0] + boundaries[0][0]) / 2,
               -(boundaries[1][1] + boundaries[0][1]) / 2,
               -(boundaries[1][2] + boundaries[0][2]) / 2);

  glBegin(GL_TRIANGLES);
  for(const array<array<GLfloat, 3>, 2>& pointNNormal : *frameIter)
  {
    glNormal3fv(pointNNormal[1].data());
    glVertex3fv(pointNNormal[0].data());
  }
  glEnd();

  /*
  glBegin(GL_QUADS);
    glNormal3f(0, 0, 1);
    glVertex3f(-1, 1, 1);
    glVertex3f(1, 1, 1);
    glVertex3f(1, -1, 1);
    glVertex3f(-1, -1, 1);

    glNormal3f(0, 1, 0);
    glVertex3f(-1, 1, 1);
    glVertex3f(-1, 1, -1);
    glVertex3f(1, 1, -1);
    glVertex3f(1, 1, 1);

    glNormal3f(1, 0, 0);
    glVertex3f(1, 1, 1);
    glVertex3f(1, 1, -1);
    glVertex3f(1, -1, -1);
    glVertex3f(1, -1, 1);
  glEnd();
  */

  glPopMatrix();
  glfwGetMousePos(&lastMousePosition[0], &lastMousePosition[1]);

  frameIter++;
  if(frameIter == end(frames))
    frameIter = begin(frames);
}

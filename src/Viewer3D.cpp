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

#include <iostream>
#include <tuple>

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glfw.h>

using namespace std;
using namespace PstpFinder::MarchingCubes;

Viewer3D::Viewer3D() :
    boundariesSet(false)
{
  glfwInit();
}

Viewer3D::~Viewer3D()
{
  glfwTerminate();
}

void
Viewer3D::appendFrame(const vector<PointAndNormal>& frame)
{
  frames.push_back(frame);
}

void
Viewer3D::run()
{
  /* Calculate boundaries */
  if(frames.size() > 0)
  {
    boundaries[0] = {{ frames[0][0].x(), frames[0][0].y(), frames[0][0].z() }};
    boundaries[1] = boundaries[0];
  }
  else
    boundaries = array<Vector3d, 2>{{ {{ 0, 0, 0 }}, {{ 0, 0, 0 }} }};

  for(const vector<PointAndNormal>& frame : frames)
  {
    for(const PointAndNormal& point : frame)
    {
      if(point.x() < boundaries[0].x())
        boundaries[0].x() = point.x();
      if(point.x() > boundaries[1].x())
        boundaries[1].x() = point.x();
      if(point.y() < boundaries[0].y())
        boundaries[0].y() = point.y();
      if(point.y() > boundaries[1].y())
        boundaries[1].y() = point.y();
      if(point.z() < boundaries[0].z())
        boundaries[0].z() = point.z();
      if(point.z() > boundaries[1].z())
        boundaries[1].z() = point.z();
    }
  }

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
  rotationXY = {{ 0, 0 }};
  glClearColor(0, 0, 0, 0);
  glShadeModel(GL_SMOOTH);
  glClearDepth(1);

  static GLfloat light_diffuse[] = {1.0, 1.0, 1.0, 1.0};
  static GLfloat light_position[] =
  { (boundaries[1].x() + boundaries[0].x()) / 2,
    (boundaries[1].y()+ boundaries[0].y()) / 2 + 2.f,
    boundaries[1].z() + 2.f, 0.0};
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
                 1, 30);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(0, 0, boundaries[1].z() + 1.f, 0, 0, 0, 0, 1, 0);
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
    rotationXY[0] += mousePosition[1] - lastMousePosition[1];
    rotationXY[1] += mousePosition[0] - lastMousePosition[0];
  }

  static auto frameIter = begin(frames);

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glColor3f(1, 1, 1);

  glPushMatrix();
  glRotatef(rotationXY[0], 1, 0, 0);
  glRotatef(rotationXY[1], 0, 1, 0);
  glTranslatef(-(boundaries[1].x() + boundaries[0].x()) / 2,
               -(boundaries[1].y() + boundaries[0].y()) / 2,
               -(boundaries[1].z() + boundaries[0].z()) / 2);

  glBegin(GL_TRIANGLES);
  for(const PointAndNormal& pointNNormal : *frameIter)
  {
    glNormal3fv(pointNNormal.coords.data() + 3);
    glVertex3fv(pointNNormal.coords.data());
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

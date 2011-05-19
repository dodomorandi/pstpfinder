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

#ifndef _MAINWINDOW_H
#define _MAINWINDOW_H

#include <gtkmm.h>

#include "NewAnalysis.h"

using namespace Gtk;

class MainWindow: public Window
{
public:
  MainWindow();
  ~MainWindow();
  void createNewAnalysis();
  void destroyNewAnalysis();
  bool onDeleteEvent(GdkEventAny* event);
protected:
  HButtonBox buttonBox;
  Button buttonNew, buttonOpen;
private:
  void init();
  NewAnalysis* newAnalysis;
};

#endif

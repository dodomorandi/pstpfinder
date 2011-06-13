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

#ifndef _NEWANALYSIS_H
#define _NEWANALYSIS_H

#include "Gromacs.h"

#include <gtkmm.h>
#include <glibmm.h>

using namespace Gtk;

class NewAnalysis: public Window
{
public:
  enum parallelOperation
  {
    OPERATION_WAIT
  };

  NewAnalysis();
  NewAnalysis(Window& parent);
  void openSessionFile(const string& sessionFileName);
  Glib::Dispatcher signal_start_spin;
  Glib::Dispatcher signal_stop_spin;
  Glib::Dispatcher signal_update_limits;
private:
  Frame mainFrame;
  VBox vboxFrame1, vboxFrame2, vboxMain;
  FileChooserButton trjChooser, tprChooser;
  Label labelTrajectory, labelTopology, labelBegin, labelEnd, labelRadius,
        labelPocketThreshold, labelPs, labelAngstrom, labelSessionFile;
  HBox hboxTrajectory, hboxTopology, hboxBegin, hboxEnd, hboxFrame, hboxRadius,
       hboxPocketThreshold, hboxSession;
  Entry entrySessionFile;
  HButtonBox buttonBoxRun;
  Button  buttonRun, buttonBrowseFile;
  ProgressBar progress;
  SpinButton spinBegin, spinEnd, spinRadius, spinPocketThreshold;
  HScale hScaleBegin, hScaleEnd;
  Spinner spinnerWait;
  VSeparator vSeparator;
  int __frames;
  float __timeStep;

  void init();
  void runAnalysis();
  void chooserTrajectoryClicked();
  void threadTrajectoryClicked();
  void start_spin();
  void stop_spin();
  void checkParameters();
  void buttonBrowseFileClicked();
  void update_limits();
  Gromacs::Pittpi runPittpi(Gromacs::Gromacs& gromacs,
                            const string& analysisFileName, float radius,
                            float threshold);
};

#endif

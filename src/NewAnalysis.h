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
#include "Pittpi.h"
#include "Results.h"

#include <vector>
#include <gtkmm.h>
#include <glibmm.h>

using namespace Gtk;

namespace PstpFinder
{
  class Results;

  class NewAnalysis :
      public Window
  {
    public:
      enum parallelOperation
      {
        OPERATION_WAIT
      };

      NewAnalysis();
      NewAnalysis(Window& parent);
      void openSessionFile(const string& sessionFileName);
      void deleteResultsWindow(const Results& resultsWindow);
      Glib::Dispatcher signal_start_spin;
      Glib::Dispatcher signal_stop_spin;
      Glib::Dispatcher signal_update_limits;
    private:
      enum enumAnalysisStatus
      {
        ANALYSIS_NOT_STARTED,
        ANALYSIS_ONGOING,
        ANALYSIS_FINISHED
      };

      Frame mainFrame;
      VBox vboxFrame1, vboxFrame2, vboxMain;
      FileChooserButton trjChooser, tprChooser;
      Label labelTrajectory, labelTopology, labelBegin, labelEnd, labelRadius,
            labelPocketThreshold, labelPs, labelAngstrom, labelSessionFile;
      HBox hboxTrajectory, hboxTopology, hboxBegin, hboxEnd, hboxFrame, hboxRadius,
           hboxPocketThreshold, hboxSession;
      Entry entrySessionFile;
      HButtonBox buttonBoxRun;
      Button buttonRun, buttonBrowseFile, buttonShowResults;
      ProgressBar progress;
      Alignment progressAligner;
      SpinButton spinBegin, spinEnd, spinRadius, spinPocketThreshold;
      HScale hScaleBegin, hScaleEnd;
      Spinner spinnerWait;
      VSeparator vSeparator;
      Statusbar statusBar;
      unsigned int statusBarContext;
      shared_ptr<Pittpi> pittpiPtr;
      Gromacs* gromacs;
      enumAnalysisStatus analysisStatus;
      int __frames;
      float __timeStep;
      std::vector<Results*> resultsWindows;
      bool abortFlag;

      void init();
      void runAnalysis() throw();
      void chooserTrajectoryClicked() throw();
      void threadTrajectoryClicked() throw();
      void start_spin() throw();
      void stop_spin() throw();
      void checkParameters();
      void buttonBrowseFileClicked() throw();
      void buttonShowResultsClicked() throw();
      void update_limits() throw();
      bool close_window(GdkEventAny* event) throw();
      void runPittpi(const string& SessionFileName,
                                    float radius, float threshold);
  };
}

#endif

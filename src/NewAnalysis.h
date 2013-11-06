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

namespace PstpFinder
{
  class Results;

  class NewAnalysis : public Gtk::Window
  {
    public:
      enum parallelOperation
      {
        OPERATION_WAIT
      };

      NewAnalysis();
      NewAnalysis(Gtk::Window& parent);
      void openSessionFile(const std::string& sessionFileName);
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

      Gtk::Frame mainFrame;
      Gtk::VBox vboxFrame1, vboxFrame2, vboxMain;
      Gtk::FileChooserButton trjChooser, tprChooser;
      Gtk::Label labelTrajectory, labelTopology, labelBegin, labelEnd,
          labelRadius, labelPocketThreshold, labelPs, labelAngstrom,
          labelSessionFile;
      Gtk::HBox hboxTrajectory, hboxTopology, hboxBegin, hboxEnd, hboxFrame,
          hboxRadius, hboxPocketThreshold, hboxSession;
      Gtk::Entry entrySessionFile;
      Gtk::HButtonBox buttonBoxRun, buttonBoxBrowse;
      Gtk::Button buttonRun, buttonBrowseFile, buttonShowResults;
      Gtk::ProgressBar progress;
      Gtk::Alignment progressAligner;
      Gtk::SpinButton spinBegin, spinEnd, spinRadius, spinPocketThreshold;
      Gtk::HScale hScaleBegin, hScaleEnd;
      Gtk::Spinner spinnerWait;
      Gtk::VSeparator vSeparator;
      Gtk::Statusbar statusBar;
      unsigned int statusBarContext;
      std::shared_ptr<Pittpi> pittpiPtr;
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
      void runPittpi(const std::string& SessionFileName, float radius,
                     float threshold);
      template<typename Session>
      void calculateSas(Session& session);

      template<typename Session>
      void calculateAverageStructure(Session& session);
  };
}

#endif

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

#include "pstpfinder.h"
#include "NewAnalysis.h"
#include "Gromacs.h"
#include "Pittpi.h"
#include "Results.h"
#include "Session.h"
#include "GtkmmWrapper.h"
#include "utils.h"
#include "Pdb.h"

#include <gtkmm.h>
#include <glibmm.h>
#include <gdkmm.h>
#include <fstream>
#include <string>

namespace PstpFinder
{
  NewAnalysis::NewAnalysis()
  {
    init();
    show_all();
  }

  NewAnalysis::NewAnalysis(Window& parent)
  {
    set_transient_for(parent);
    init();
    show_all();
  }

  void
  NewAnalysis::init()
  {
    signal_start_spin.connect(sigc::mem_fun(*this, &NewAnalysis::start_spin));
    signal_stop_spin.connect(sigc::mem_fun(*this, &NewAnalysis::stop_spin));
    signal_update_limits.connect(
        sigc::mem_fun(*this, &NewAnalysis::update_limits));

    set_title("PSTP-finder");

    GtkmmWrapper<Gtk::FileFilter> trjFilter;
    trjFilter->set_name("Trajectory files");
    trjFilter->add_pattern("*.xtc");
    trjChooser.add_filter(trjFilter);
    trjChooser.signal_file_set().connect(
        sigc::mem_fun(*this, &NewAnalysis::chooserTrajectoryClicked));
    trjChooser.set_size_request(150, -1);

    labelTrajectory.set_label("Trajectory file:");

    hboxTrajectory.pack_start(labelTrajectory, Gtk::PACK_SHRINK);
    hboxTrajectory.pack_start(trjChooser);
    hboxTrajectory.set_homogeneous(false);
    hboxTrajectory.set_spacing(10);

    GtkmmWrapper<Gtk::FileFilter> tprFilter;
    tprFilter->set_name("Topology files");
    tprFilter->add_pattern("*.tpr");
    tprChooser.add_filter(tprFilter);
    tprChooser.set_size_request(150, -1);
    tprChooser.signal_file_set().connect(
        sigc::mem_fun(*this, &NewAnalysis::checkParameters));

    labelTopology.set_label("Topology file:");

    hboxTopology.pack_start(labelTopology, Gtk::PACK_SHRINK);
    hboxTopology.pack_start(tprChooser);
    hboxTopology.set_homogeneous(false);
    hboxTopology.set_spacing(10);

    labelBegin.set_label("Begin:");
#if GTKMM_MAJOR == 3
    hScaleBegin.set_adjustment(spinBegin.get_adjustment());
#else
    hScaleBegin.set_adjustment(*spinBegin.get_adjustment());
#endif
    hScaleBegin.set_draw_value(false);
    hScaleBegin.set_sensitive(false);
    spinBegin.set_sensitive(false);
    spinBegin.set_digits(3);
    hboxBegin.pack_start(labelBegin, Gtk::PACK_SHRINK);
    hboxBegin.pack_start(hScaleBegin);
    hboxBegin.pack_start(spinBegin, Gtk::PACK_SHRINK);
    hboxBegin.set_homogeneous(false);
    hboxBegin.set_spacing(10);

    labelEnd.set_label("End:");
#if GTKMM_MAJOR == 3
    hScaleEnd.set_adjustment(spinEnd.get_adjustment());
#else
    hScaleEnd.set_adjustment(*spinEnd.get_adjustment());
#endif
    hScaleEnd.set_draw_value(false);
    hScaleEnd.set_sensitive(false);
    spinEnd.set_sensitive(false);
    spinEnd.set_digits(3);
    hboxEnd.pack_start(labelEnd, Gtk::PACK_SHRINK);
    hboxEnd.pack_start(hScaleEnd);
    hboxEnd.pack_start(spinEnd, Gtk::PACK_SHRINK);
    hboxEnd.set_homogeneous(false);
    hboxEnd.set_spacing(10);

    vboxFrame1.set_spacing(10);
    vboxFrame1.pack_start(hboxTrajectory, Gtk::PACK_EXPAND_PADDING);
    vboxFrame1.pack_start(hboxTopology, Gtk::PACK_EXPAND_PADDING);
    vboxFrame1.pack_start(hboxBegin, Gtk::PACK_EXPAND_PADDING);
    vboxFrame1.pack_start(hboxEnd, Gtk::PACK_EXPAND_PADDING);

    labelRadius.set_label("Pocket radius:");
    spinRadius.set_digits(1);
    spinRadius.set_range(0.0, 20.0);
    spinRadius.set_value(7.0);
    spinRadius.set_increments(0.1, 1.0);
    labelAngstrom.set_label("Angstrom");
    hboxRadius.set_spacing(10);
    hboxRadius.pack_start(labelRadius, Gtk::PACK_SHRINK);
    hboxRadius.pack_start(spinRadius);
    hboxRadius.pack_start(labelAngstrom, Gtk::PACK_SHRINK);

    labelPocketThreshold.set_label("Pocket threshold:");
    spinPocketThreshold.set_digits(0);
    spinPocketThreshold.set_increments(1, 10);
    spinPocketThreshold.set_range(0, 20000);
    spinPocketThreshold.set_value(500);
    labelPs.set_label("ps");
    hboxPocketThreshold.set_spacing(10);
    hboxPocketThreshold.pack_start(labelPocketThreshold, Gtk::PACK_SHRINK);
    hboxPocketThreshold.pack_start(spinPocketThreshold);
    hboxPocketThreshold.pack_start(labelPs, Gtk::PACK_SHRINK);

    labelSessionFile.set_label("Session file:");
    buttonBrowseFile.set_label("Browse...");
    buttonBrowseFile.signal_clicked().connect(
        sigc::mem_fun(*this, &NewAnalysis::buttonBrowseFileClicked));
    buttonBoxBrowse.pack_start(buttonBrowseFile, Gtk::PACK_SHRINK);
    hboxSession.set_spacing(10);
    hboxSession.set_homogeneous(false);
    hboxSession.pack_start(labelSessionFile, Gtk::PACK_SHRINK);
    hboxSession.pack_start(entrySessionFile);
    hboxSession.pack_start(buttonBoxBrowse, Gtk::PACK_SHRINK);

    vboxFrame2.set_spacing(10);
    vboxFrame2.pack_start(hboxRadius);
    vboxFrame2.pack_start(hboxPocketThreshold);
    vboxFrame2.pack_start(hboxSession);

    hboxFrame.set_spacing(10);
    hboxFrame.set_border_width(10);
    hboxFrame.pack_start(vboxFrame1);
    hboxFrame.pack_start(vSeparator, Gtk::PACK_SHRINK);
    hboxFrame.pack_start(vboxFrame2);

    mainFrame.set_label("Main files");
    mainFrame.set_border_width(10);
    mainFrame.add(hboxFrame);

    buttonRun.set_label("Run!");
    buttonRun.set_sensitive(false);
    buttonRun.signal_clicked().connect(
        sigc::mem_fun(*this, &NewAnalysis::runAnalysis));
    buttonShowResults.set_label("Show results");
    buttonShowResults.set_sensitive(false);
    buttonShowResults.signal_clicked().connect(
        sigc::mem_fun(*this, &NewAnalysis::buttonShowResultsClicked));
    buttonBoxRun.set_layout(Gtk::BUTTONBOX_END);
    buttonBoxRun.set_border_width(10);
    buttonBoxRun.set_spacing(10);
    buttonBoxRun.pack_end(buttonShowResults);
    buttonBoxRun.pack_end(buttonRun, Gtk::PACK_SHRINK);

    progressAligner.add(progress);
    progressAligner.set_border_width(10);

    statusBarContext = statusBar.get_context_id("PSTP-finder status");
    statusBar.push("Welcome to PSTP-finder. Choose you options and press Run.",
                   statusBarContext);

    vboxMain.set_homogeneous(false);
    vboxMain.pack_start(mainFrame);
    vboxMain.pack_start(buttonBoxRun, Gtk::PACK_SHRINK);
    vboxMain.pack_start(progressAligner, Gtk::PACK_SHRINK);
    vboxMain.pack_start(statusBar, Gtk::PACK_SHRINK);

    add(vboxMain);
    abortFlag = false;
    analysisStatus = enumAnalysisStatus::ANALYSIS_NOT_STARTED;

    signal_delete_event().connect(
        sigc::mem_fun(*this, &NewAnalysis::close_window));
  }

  void
  NewAnalysis::start_spin() throw()
  {
    if(spinnerWait.get_parent() == 0)
      vboxMain.pack_start(spinnerWait);
    mainFrame.hide();
    buttonBoxRun.hide();
    progressAligner.hide();
    statusBar.hide();
    spinnerWait.show();
    spinnerWait.start();
  }

  void
  NewAnalysis::stop_spin() throw()
  {
    spinnerWait.stop();
    spinnerWait.hide();
    mainFrame.show();
    buttonBoxRun.show();
    progressAligner.show();
    statusBar.show();
  }

  void
  NewAnalysis::runPittpi(const std::string& SessionFileName, float radius,
                         float threshold)
  {
    progress.set_fraction(0);
    progress.set_pulse_step(0.1);
    while(Gtk::Main::events_pending())
      Gtk::Main::iteration();

    pittpiPtr = std::shared_ptr<Pittpi>
      { new Pittpi(*gromacs, SessionFileName, radius, threshold) };

    while(not pittpiPtr->isFinished())
    {
      float status = pittpiPtr->getStatus();
      float oldStatus = progress.get_fraction();
      if(status <= oldStatus)
        statusBar.push(pittpiPtr->getStatusDescription(), statusBarContext);
      if(status >= 0)
        progress.set_fraction(status);
      else
        progress.pulse();
      while(Gtk::Main::events_pending())
        Gtk::Main::iteration();
      pittpiPtr->waitNextStatus();
    }

    if(not abortFlag)
    {
      progress.set_fraction(pittpiPtr->getStatus());
      statusBar.push(pittpiPtr->getStatusDescription(), statusBarContext);
    }
  }

  void
  NewAnalysis::runAnalysis() throw()
  {
    Glib::ustring sessionFileName(entrySessionFile.get_text());
    if(sessionFileName.empty())
    {
      Gtk::MessageDialog msg("A valid session file is needed.", false,
                        Gtk::MessageType::MESSAGE_ERROR,
                        Gtk::ButtonsType::BUTTONS_OK);
      msg.run();
      return;
    }

    if(exists(sessionFileName))
    {
      Gtk::MessageDialog question("The old session file will be overwritten. "
                        "Are you sure?", false,
                        Gtk::MessageType::MESSAGE_QUESTION,
                        Gtk::ButtonsType::BUTTONS_YES_NO);
      int returnValue(question.run());
      if(returnValue == Gtk::ResponseType::RESPONSE_NO)
        return;
    }

    if(analysisStatus == enumAnalysisStatus::ANALYSIS_FINISHED)
      delete gromacs;
    analysisStatus = enumAnalysisStatus::ANALYSIS_ONGOING;
    gromacs = new Gromacs(trjChooser.get_filename(), tprChooser.get_filename());

    gromacs->setBegin(spinBegin.get_value());
    gromacs->setEnd(spinEnd.get_value());

    mainFrame.set_sensitive(false);
    buttonShowResults.set_sensitive(false);
    buttonRun.set_sensitive(false);

    Session<std::ofstream> session(sessionFileName, *gromacs,
                                   spinRadius.get_value(),
                                   spinPocketThreshold.get_value());

    calculateSas(session);
    if(abortFlag)
      return;

    calculateAverageStructure(session);
    if(abortFlag)
      return;

    runPittpi(sessionFileName, spinRadius.get_value(),
              spinPocketThreshold.get_value());
    if(abortFlag)
      return;

    pittpiPtr->save(session.getPittpiStream());

    mainFrame.set_sensitive();
    buttonRun.set_sensitive();

    buttonShowResults.set_sensitive();
    resultsWindows.push_back(new Results(*this, pittpiPtr, *gromacs));

    analysisStatus = enumAnalysisStatus::ANALYSIS_FINISHED;
  }

  void
  NewAnalysis::chooserTrajectoryClicked() throw()
  {
    if(not exists(trjChooser.get_filename()))
    {
      spinBegin.set_sensitive(false);
      spinEnd.set_sensitive(false);
      hScaleBegin.set_sensitive(false);
      hScaleEnd.set_sensitive(false);
    }
    else
      Glib::Thread::create(
          sigc::mem_fun(*this, &NewAnalysis::threadTrajectoryClicked), false);

    checkParameters();
  }

  void
  NewAnalysis::threadTrajectoryClicked() throw()
  {
    signal_start_spin();

    Gromacs tmpGromacs(trjChooser.get_filename(), "");
    __timeStep = tmpGromacs.getTimeStep();
    __frames = tmpGromacs.getFramesCount();

    signal_stop_spin();
    signal_update_limits();
  }

  void
  NewAnalysis::update_limits() throw()
  {
    spinBegin.set_sensitive();
    spinEnd.set_sensitive();
    hScaleBegin.set_sensitive();
    hScaleEnd.set_sensitive();

    spinBegin.set_range(0, (__frames - 1) * __timeStep);
    spinBegin.set_value(0);
    spinBegin.set_increments(__timeStep, (__frames - 1) / 100);

    spinEnd.set_range(0, (__frames - 1) * __timeStep);
    spinEnd.set_value((__frames - 1) * __timeStep);
    spinEnd.set_increments(__timeStep, (__frames - 1) / 100);
  }

  void
  NewAnalysis::checkParameters()
  {
    if(exists(tprChooser.get_filename()) and exists(trjChooser.get_filename()))
      buttonRun.set_sensitive(true);
    else
      buttonRun.set_sensitive(false);
  }

  void
  NewAnalysis::buttonBrowseFileClicked() throw()
  {
    GtkmmWrapper<Gtk::FileFilter> filter;
    filter->add_pattern("*.csf");
    filter->set_name("PSTP-filter compressed session file");

    Gtk::FileChooserDialog chooser("Choose a saving file for this session",
                                   Gtk::FILE_CHOOSER_ACTION_SAVE);
    chooser.add_filter(filter);
    chooser.add_button(Gtk::Stock::CANCEL, Gtk::RESPONSE_CANCEL);
    chooser.add_button(Gtk::Stock::OK, Gtk::RESPONSE_OK);
    int response = chooser.run();

    switch(response)
    {
      case Gtk::RESPONSE_OK:
      {
        std::string filename = chooser.get_filename();
        if(file_extension(filename) != ".csf")
          filename = change_extension(filename, ".csf");

        entrySessionFile.set_text(filename);
        break;
      }
    }
  }

  void
  NewAnalysis::buttonShowResultsClicked() throw()
  {
    if(analysisStatus != enumAnalysisStatus::ANALYSIS_FINISHED)
      return;

    resultsWindows.push_back(new Results(*this, pittpiPtr, *gromacs));
  }

  void
  NewAnalysis::openSessionFile(const std::string& sessionFileName)
  {
    double beginTime, endTime;

    start_spin();
    while(Gtk::Main::events_pending())
      Gtk::Main::iteration();

    Session<std::fstream> session(sessionFileName);

    trjChooser.set_filename(session.getTrajectoryFileName());
    tprChooser.set_filename(session.getTopologyFileName());

    beginTime = session.getBeginTime();
    endTime = session.getEndTime();
    spinRadius.set_value(session.getRadius());
    spinPocketThreshold.set_value(session.getPocketThreshold());
    entrySessionFile.set_text(sessionFileName);

    mainFrame.set_sensitive(false);
    buttonShowResults.set_sensitive(false);
    buttonRun.set_sensitive(false);
    start_spin();
    while(Gtk::Main::events_pending())
      Gtk::Main::iteration();

    analysisStatus = enumAnalysisStatus::ANALYSIS_ONGOING;

    gromacs = new Gromacs(session.getTrajectoryFileName(),
                          session.getTopologyFileName());
    __timeStep = gromacs->getTimeStep();
    __frames = gromacs->getFramesCount();
    gromacs->setBegin(beginTime);
    gromacs->setEnd(endTime);
    spinBegin.set_value(beginTime);
    spinEnd.set_value(endTime);

    if(not session.sasComplete())
    {
      stop_spin();

      calculateSas(session);
      if(abortFlag)
        return;

      calculateAverageStructure(session);
      if(abortFlag)
        return;
    }
    else
    {
      while(Gtk::Main::events_pending())
        Gtk::Main::iteration();
      gromacs->setAverageStructure(Pdb<>(session.getPdbStream()).proteins[0]);

      stop_spin();
      while(Gtk::Main::events_pending())
        Gtk::Main::iteration();
    }

    if(not session.pittpiComplete())
    {
      progress.set_fraction(0);
      while(Gtk::Main::events_pending())
        Gtk::Main::iteration();

      runPittpi(sessionFileName, spinRadius.get_value(),
                spinPocketThreshold.get_value());
      if(abortFlag)
        return;

      pittpiPtr->save(session.getPittpiStream());
    }
    else
    {
      progress.set_fraction(0);
      while(Gtk::Main::events_pending())
        Gtk::Main::iteration();

      pittpiPtr = std::shared_ptr<Pittpi>(
          new Pittpi(*gromacs, sessionFileName, session.getRadius(),
                     session.getPocketThreshold(), false));
      pittpiPtr->load(session.getPittpiStream());
    }

    mainFrame.set_sensitive();
    buttonRun.set_sensitive();
    update_limits();

    progress.set_fraction(0);
    progress.show();
    while(Gtk::Main::events_pending())
      Gtk::Main::iteration();

    resultsWindows.push_back(new Results(*this, pittpiPtr, *gromacs));
    
    buttonShowResults.set_sensitive();
    analysisStatus = enumAnalysisStatus::ANALYSIS_FINISHED;
  }

  bool
  NewAnalysis::close_window(GdkEventAny* event) throw()
  {
    if(not mainFrame.is_sensitive())
    {
      Gtk::MessageDialog messageAbort(
          "An analysis is running. If you quit you can resume it from "
          "'Open/resume analysis' button. Abort and quit?",
          false, Gtk::MESSAGE_QUESTION, Gtk::BUTTONS_YES_NO, true);

      int response = messageAbort.run();
      if(response == Gtk::RESPONSE_YES)
      {
        abortFlag = true;
        if(pittpiPtr)
          pittpiPtr->abort();
        if(analysisStatus == enumAnalysisStatus::ANALYSIS_ONGOING)
        {
          gromacs->abort();
          analysisStatus = enumAnalysisStatus::ANALYSIS_NOT_STARTED;
        }
        closeApplication(event);
        return false;
      }
      else
        return true;
    }

    closeApplication(event);
    return false;
  }

  void
  NewAnalysis::deleteResultsWindow(const Results& resultsWindow)
  {
    for
    (
      std::vector<Results*>::iterator i = resultsWindows.begin();
      i < resultsWindows.end();
      i++
    )
      if(*i == &resultsWindow)
      {
        delete *i;
        resultsWindows.erase(i);
        break;
      }
  }

  template<typename Session>
  void
  NewAnalysis::calculateSas(Session& session)
  {
    statusBar.push("Calculating SAS using Gromacs", statusBarContext);
    progress.set_fraction(0);
    while(Gtk::Main::events_pending())
      Gtk::Main::iteration();

    unsigned int currentFrame;
    unsigned int count = gromacs->getFramesCount();
    while(Gtk::Main::events_pending())
      Gtk::Main::iteration();

    if(abortFlag)
      return;
    gromacs->calculateSas(session);

    while((currentFrame = gromacs->getCurrentFrame()) < count)
    {
      if(abortFlag)
      {
        gromacs->abort();
        break;
      }
      progress.set_fraction(static_cast<float>(currentFrame) / count);
      while(Gtk::Main::events_pending())
        Gtk::Main::iteration();
      gromacs->waitNextFrame();
    }
    if(abortFlag)
      return;

    gromacs->waitOperation();
    while(Gtk::Main::events_pending())
      Gtk::Main::iteration();

    progress.set_fraction(1);
    while(Gtk::Main::events_pending())
      Gtk::Main::iteration();
  }

  template<typename Session>
  void
  NewAnalysis::calculateAverageStructure(Session& session)
  {
    statusBar.push("Calculating average structure", statusBarContext);
    progress.set_fraction(0);
    while(Gtk::Main::events_pending())
      Gtk::Main::iteration();

    unsigned int currentFrame;
    unsigned int count(gromacs->getFramesCount());
    gromacs->calculateAverageStructure();

    while((currentFrame = gromacs->getCurrentFrame()) < count)
    {
      if(abortFlag)
      {
        gromacs->abort();
        break;
      }
      progress.set_fraction(static_cast<float>(currentFrame) / count);
      while(Gtk::Main::events_pending())
        Gtk::Main::iteration();
      gromacs->waitNextFrame();
    }
    if(abortFlag)
      return;

    gromacs->waitOperation();
    while(Gtk::Main::events_pending())
      Gtk::Main::iteration();

    if(abortFlag)
      return;

    Pdb<> averagePdb ;
    averagePdb.proteins.push_back(gromacs->getAverageStructure());
    averagePdb.write(session.getPdbStream());

    if(abortFlag)
      return;

    progress.set_fraction(1);
    while(Gtk::Main::events_pending())
      Gtk::Main::iteration();
  }
}

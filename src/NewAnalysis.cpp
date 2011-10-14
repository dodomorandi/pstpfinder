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
#include "Session.h"

#include <gtkmm.h>
#include <glibmm.h>
#include <gdkmm.h>
#include <fstream>
#include <string>
#include <boost/filesystem.hpp>

using namespace Gtk;
namespace fs = boost::filesystem;

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

    FileFilter trjFilter;
    trjFilter.set_name("Trajectory files");
    trjFilter.add_pattern("*.xtc");
    trjFilter.add_pattern("*.trj");
    trjChooser.add_filter(trjFilter);
    trjChooser.signal_file_set().connect(
        sigc::mem_fun(*this, &NewAnalysis::chooserTrajectoryClicked));
    trjChooser.set_size_request(150, -1);

    labelTrajectory.set_label("Trajectory file:");

    hboxTrajectory.pack_start(labelTrajectory, PACK_SHRINK);
    hboxTrajectory.pack_start(trjChooser);
    hboxTrajectory.set_homogeneous(false);
    hboxTrajectory.set_spacing(10);

    FileFilter tprFilter;
    tprFilter.set_name("Topology files");
    tprFilter.add_pattern("*.tpr");
    tprChooser.add_filter(tprFilter);
    tprChooser.set_size_request(150, -1);
    tprChooser.signal_file_set().connect(
        sigc::mem_fun(*this, &NewAnalysis::checkParameters));

    labelTopology.set_label("Topology file:");

    hboxTopology.pack_start(labelTopology, PACK_SHRINK);
    hboxTopology.pack_start(tprChooser);
    hboxTopology.set_homogeneous(false);
    hboxTopology.set_spacing(10);

    labelBegin.set_label("Begin:");
    hScaleBegin.set_adjustment(*spinBegin.get_adjustment());
    hScaleBegin.set_draw_value(false);
    hScaleBegin.set_sensitive(false);
    spinBegin.set_sensitive(false);
    spinBegin.set_digits(3);
    hboxBegin.pack_start(labelBegin, PACK_SHRINK);
    hboxBegin.pack_start(hScaleBegin);
    hboxBegin.pack_start(spinBegin, PACK_SHRINK);
    hboxBegin.set_homogeneous(false);
    hboxBegin.set_spacing(10);

    labelEnd.set_label("End:");
    hScaleEnd.set_adjustment(*spinEnd.get_adjustment());
    hScaleEnd.set_draw_value(false);
    hScaleEnd.set_sensitive(false);
    spinEnd.set_sensitive(false);
    spinEnd.set_digits(3);
    hboxEnd.pack_start(labelEnd, PACK_SHRINK);
    hboxEnd.pack_start(hScaleEnd);
    hboxEnd.pack_start(spinEnd, PACK_SHRINK);
    hboxEnd.set_homogeneous(false);
    hboxEnd.set_spacing(10);

    vboxFrame1.set_spacing(10);
    vboxFrame1.pack_start(hboxTrajectory, PACK_EXPAND_PADDING);
    vboxFrame1.pack_start(hboxTopology, PACK_EXPAND_PADDING);
    vboxFrame1.pack_start(hboxBegin, PACK_EXPAND_PADDING);
    vboxFrame1.pack_start(hboxEnd, PACK_EXPAND_PADDING);

    labelRadius.set_label("Pocket radius:");
    spinRadius.set_digits(1);
    spinRadius.set_range(0.0, 20.0);
    spinRadius.set_value(7.0);
    spinRadius.set_increments(0.1, 1.0);
    labelAngstrom.set_label("Angstrom");
    hboxRadius.set_spacing(10);
    hboxRadius.pack_start(labelRadius, PACK_SHRINK);
    hboxRadius.pack_start(spinRadius);
    hboxRadius.pack_start(labelAngstrom, PACK_SHRINK);

    labelPocketThreshold.set_label("Pocket threshold:");
    spinPocketThreshold.set_digits(0);
    spinPocketThreshold.set_increments(1, 10);
    spinPocketThreshold.set_range(0, 20000);
    spinPocketThreshold.set_value(500);
    labelPs.set_label("ps");
    hboxPocketThreshold.set_spacing(10);
    hboxPocketThreshold.pack_start(labelPocketThreshold, PACK_SHRINK);
    hboxPocketThreshold.pack_start(spinPocketThreshold);
    hboxPocketThreshold.pack_start(labelPs, PACK_SHRINK);

    labelSessionFile.set_label("Session file:");
    buttonBrowseFile.set_label("Browse...");
    buttonBrowseFile.signal_clicked().connect(
        sigc::mem_fun(*this, &NewAnalysis::buttonBrowseFileClicked));
    hboxSession.set_spacing(10);
    hboxSession.set_homogeneous(false);
    hboxSession.pack_start(labelSessionFile, PACK_SHRINK);
    hboxSession.pack_start(entrySessionFile);
    hboxSession.pack_start(buttonBrowseFile);

    vboxFrame2.set_spacing(10);
    vboxFrame2.pack_start(hboxRadius);
    vboxFrame2.pack_start(hboxPocketThreshold);
    vboxFrame2.pack_start(hboxSession);

    hboxFrame.set_spacing(10);
    hboxFrame.set_border_width(10);
    hboxFrame.pack_start(vboxFrame1);
    hboxFrame.pack_start(vSeparator, PACK_SHRINK);
    hboxFrame.pack_start(vboxFrame2);

    mainFrame.set_label("Main files");
    mainFrame.add(hboxFrame);

    buttonRun.set_label("Run!");
    buttonRun.set_sensitive(false);
    buttonRun.signal_clicked().connect(
        sigc::mem_fun(*this, &NewAnalysis::runAnalysis));
    buttonBoxRun.set_layout(BUTTONBOX_END);
    buttonBoxRun.pack_end(buttonRun);

    vboxMain.set_homogeneous(false);
    vboxMain.set_spacing(10);
    vboxMain.set_border_width(10);
    vboxMain.pack_start(mainFrame);
    vboxMain.pack_start(buttonBoxRun, PACK_SHRINK);

    add(vboxMain);
    pittpi = 0;

    signal_delete_event().connect(sigc::ptr_fun(&closeApplication));
  }

  void
  NewAnalysis::start_spin()
  {
    if(spinnerWait.get_parent() == 0)
      vboxMain.pack_start(spinnerWait);
    mainFrame.hide();
    buttonBoxRun.hide();
    spinnerWait.show();
    spinnerWait.start();
  }

  void
  NewAnalysis::stop_spin()
  {
    spinnerWait.stop();
    spinnerWait.hide();
    mainFrame.show();
    buttonBoxRun.show();
  }

  void
  NewAnalysis::runPittpi(Gromacs& gromacs, const string& SessionFileName,
                         float radius, float threshold)
  {
    progress.set_fraction(0);
    progress.set_pulse_step(0.1);
    while(Main::events_pending())
      Main::iteration();

    pittpi = new Pittpi(gromacs, SessionFileName, radius, threshold);

    while(not pittpi->isFinished())
    {
      float status = pittpi->getStatus();
      if(status >= 0)
        progress.set_fraction(status);
      else
        progress.pulse();
      while(Main::events_pending())
        Main::iteration();
      pittpi->waitNextStatus();
    }

    delete pittpi;
    pittpi = 0;
  }

  void
  NewAnalysis::runAnalysis()
  {
    bool writeSession = true;
    std::locale oldLocale;
    std::locale::global(std::locale("C"));

    Gromacs gromacs(trjChooser.get_filename(), tprChooser.get_filename());

    gromacs.setBegin(spinBegin.get_value());
    gromacs.setEnd(spinEnd.get_value());

    if(not progress.is_ancestor(vboxMain))
    {
      vboxMain.pack_end(progress, PACK_EXPAND_WIDGET, 10);
      int x, y, width, height;
      get_position(x, y);
      get_size(width, height);
      Glib::RefPtr<Gdk::Screen> screen = get_screen();
      if(screen->get_height() < y + height)
      {
        y = screen->get_height() - height;
        move(x, y);
      }
    }

    set_sensitive(false);

    if(fs::exists(fs::path("/tmp/sas.psf")))
      fs::remove(fs::path("/tmp/sas.psf"));

    progress.set_fraction(0);
    progress.show();
    while(Main::events_pending())
      Main::iteration();

    if(entrySessionFile.get_text().empty())
      writeSession = false;

    std::ofstream sessionFile;
    if(writeSession)
    {
      sessionFile.open(entrySessionFile.get_text().c_str(),
                       std::ios::trunc | std::ios::out | std::ios::binary);

      sessionFile << trjChooser.get_filename() << endl;
      sessionFile << tprChooser.get_filename() << endl;
      sessionFile << spinBegin.get_value() << endl;
      sessionFile << spinEnd.get_value() << endl;
      sessionFile << spinRadius.get_value() << endl;
      sessionFile << spinPocketThreshold.get_value() << endl;
    }

    while(Main::events_pending())
      Main::iteration();

    unsigned int currentFrame;
    unsigned int count = gromacs.getFramesCount();
    while(Main::events_pending())
      Main::iteration();

    gromacs.calculateSas();

    while((currentFrame = gromacs.getCurrentFrame()) < count)
    {
      progress.set_fraction(static_cast<float>(currentFrame) / count);
      while(Main::events_pending())
        Main::iteration();
      gromacs.waitNextFrame();
    }

    gromacs.waitOperation();
    while(Main::events_pending())
      Main::iteration();

    if(writeSession)
    {
      sessionFile << fs::file_size(fs::path("/tmp/sas.psf")) << endl;
      std::ifstream sasFile("/tmp/sas.psf", std::ios::in | std::ios::binary);
      sessionFile << sasFile.rdbuf();
      sessionFile << endl;
    }

    progress.set_fraction(1);
    while(Main::events_pending())
      Main::iteration();

    progress.set_fraction(0);
    while(Main::events_pending())
      Main::iteration();
    gromacs.calculateAverageStructure();

    while((currentFrame = gromacs.getCurrentFrame()) < count)
    {
      progress.set_fraction(static_cast<float>(currentFrame) / count);
      while(Main::events_pending())
        Main::iteration();
      gromacs.waitNextFrame();
    }

    gromacs.waitOperation();
    while(Main::events_pending())
      Main::iteration();

    gromacs.getAverageStructure().dumpPdb("/tmp/aver.pdb");

    if(writeSession)
    {
      sessionFile << fs::file_size(fs::path("/tmp/aver.pdb")) << endl;
      std::ifstream pdbFile("/tmp/aver.pdb");
      sessionFile << pdbFile.rdbuf();
      sessionFile << endl;

      sessionFile.flush();
      sessionFile.close();
    }
    fs::remove(fs::path("/tmp/aver.pdb"));

    progress.set_fraction(1);
    while(Main::events_pending())
      Main::iteration();

    runPittpi(gromacs, "/tmp/sas.psf", spinRadius.get_value(),
              spinPocketThreshold.get_value());
    fs::remove(fs::path("/tmp/sas.psf"));

    progress.hide();
    set_sensitive(true);
    std::locale::global(oldLocale);

    MessageDialog msg("Two log files have been written in /tmp/pockets.log and "
                      "/tmp/pockets_details.log");
    msg.run();
  }

  void
  NewAnalysis::chooserTrajectoryClicked()
  {
    if(not fs::exists(fs::path(static_cast<string>(trjChooser.get_filename()))))
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
  NewAnalysis::threadTrajectoryClicked()
  {
    signal_start_spin();

    Gromacs tmpGromacs(trjChooser.get_filename(), "");
    __frames = tmpGromacs.getFramesCount();
    __timeStep = tmpGromacs.getTimeStep();

    signal_stop_spin();
    signal_update_limits();
  }

  void
  NewAnalysis::update_limits()
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
    if(fs::exists(fs::path(static_cast<string>(tprChooser.get_filename())))
       and fs::exists(fs::path(static_cast<string>(trjChooser.get_filename()))))
      buttonRun.set_sensitive(true);
    else
      buttonRun.set_sensitive(false);
  }

  void
  NewAnalysis::buttonBrowseFileClicked()
  {
    FileFilter filter;
    filter.add_pattern("*.csf");
    filter.set_name("PSTP-filter compressed session file");

    FileChooserDialog chooser("Choose a saving file for this session",
                              FILE_CHOOSER_ACTION_SAVE);
    chooser.add_filter(filter);
    chooser.add_button(Stock::CANCEL, RESPONSE_CANCEL);
    chooser.add_button(Stock::OK, RESPONSE_OK);
    int response = chooser.run();

    switch(response)
    {
      case RESPONSE_OK:
      {
        string filename = chooser.get_filename();
        if(fs::extension(fs::path(filename)) != ".csf")
          filename = fs::change_extension(fs::path(filename), ".csf").string();

        entrySessionFile.set_text(filename);
        break;
      }
    }
  }

  void
  NewAnalysis::openSessionFile(const string& sessionFileName)
  {
    std::locale oldLocale;
    std::locale::global(std::locale("C"));
    double beginTime, endTime;

    start_spin();
    while(Main::events_pending())
      Main::iteration();

    Session sessionFile(sessionFileName);
    std::string tmpString;

    trjChooser.set_filename(sessionFile.getTrajectoryFileName());
    tprChooser.set_filename(sessionFile.getTopologyFileName());

    beginTime = sessionFile.getBeginTime();
    endTime = sessionFile.getEndTime();
    spinRadius.set_value(sessionFile.getRadius());
    spinPocketThreshold.set_value(sessionFile.getPocketThreshold());
    entrySessionFile.set_text(sessionFileName);

    set_sensitive(false);
    start_spin();
    while(Main::events_pending())
      Main::iteration();

    MetaStream& pdbStream = sessionFile.getPdbStream();
    char* chunk = new char[1024 * 1024 * 128];
    unsigned long nChunks = sessionFile.getPdbSize() / (1024 * 1024 * 128);
    unsigned long remainChunk = sessionFile.getPdbSize() % (1024 * 1024 * 128);
    std::ofstream streamPdb("/tmp/aver.pdb",
                            std::ios::trunc | std::ios::out | std::ios::binary);
    for(unsigned long i = 0; i < nChunks; i++)
    {
      while(Main::events_pending())
        Main::iteration();
      pdbStream.read(chunk, 1024 * 1024 * 128);
      while(Main::events_pending())
        Main::iteration();
      streamPdb.write(chunk, 1024 * 1024 * 128);
    }

    if(remainChunk != 0)
    {
      while(Main::events_pending())
        Main::iteration();
      pdbStream.read(chunk, remainChunk);
      while(Main::events_pending())
        Main::iteration();
      streamPdb.write(chunk, remainChunk);
    }
    streamPdb.flush();
    streamPdb.close();

    delete[] chunk;

    Gromacs gromacs(sessionFile.getTrajectoryFileName(),
                    sessionFile.getTopologyFileName());
    __timeStep = gromacs.getTimeStep();
    __frames = gromacs.getFramesCount();
    gromacs.setBegin(beginTime);
    gromacs.setEnd(endTime);
    gromacs.setAverageStructure(Protein("/tmp/aver.pdb"));

    stop_spin();
    while(Main::events_pending())
      Main::iteration();

    if(not progress.is_ancestor(vboxMain))
    {
      vboxMain.pack_end(progress, PACK_EXPAND_WIDGET, 10);
      int x, y, width, height;
      get_position(x, y);
      get_size(width, height);
      Glib::RefPtr<Gdk::Screen> screen = get_screen();
      if(screen->get_height() < y + height)
      {
        y = screen->get_height() - height;
        move(x, y);
      }
    }

    progress.set_fraction(0);
    progress.show();
    while(Main::events_pending())
      Main::iteration();

    runPittpi(gromacs, sessionFileName, spinRadius.get_value(),
              spinPocketThreshold.get_value());

    progress.hide();
    std::locale::global(oldLocale);

    set_sensitive();
    update_limits();
    spinBegin.set_value(beginTime);
    spinEnd.set_value(endTime);
    buttonRun.set_sensitive();

    while(Main::events_pending())
      Main::iteration();

    MessageDialog msg("Two log files have been written in /tmp/pockets.log and "
                      "/tmp/pockets_details.log");
    msg.run();
  }

}

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

#include <gtkmm.h>

#include "pstpfinder.h"
#include "NewAnalysis.h"
#include "Gromacs.h"
#include "Pittpi.h"

#include <boost/filesystem.hpp>

using namespace Gtk;
namespace fs = boost::filesystem;

NewAnalysis::NewAnalysis(Window& parent)
{
  set_transient_for(parent);
  init();
  show_all();
}

NewAnalysis::NewAnalysis()
{
  init();
  show_all();
}

void
NewAnalysis::init()
{
  signal_start_spin.connect(sigc::mem_fun(*this, &NewAnalysis::start_spin));
  signal_stop_spin.connect(sigc::mem_fun(*this, &NewAnalysis::stop_spin));

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
  hboxRadius.set_spacing(10);
  hboxRadius.pack_start(labelRadius, PACK_SHRINK);
  hboxRadius.pack_start(spinRadius);

  labelPocketThreshold.set_label("Pocket threshold:");
  spinPocketThreshold.set_digits(0);
  spinPocketThreshold.set_increments(1, 10);
  spinPocketThreshold.set_range(0, 20000);
  labelPs.set_label("ps");
  hboxPocketThreshold.set_spacing(10);
  hboxPocketThreshold.pack_start(labelPocketThreshold, PACK_SHRINK);
  hboxPocketThreshold.pack_start(spinPocketThreshold);
  hboxPocketThreshold.pack_start(labelPs, PACK_SHRINK);

  vboxFrame2.set_spacing(10);
  vboxFrame2.pack_start(hboxRadius);
  vboxFrame2.pack_start(hboxPocketThreshold);

  hboxFrame.set_spacing(10);
  hboxFrame.set_border_width(10);
  hboxFrame.pack_start(vboxFrame1);
  hboxFrame.pack_start(vSeparator, PACK_SHRINK);
  hboxFrame.pack_start(vboxFrame2);

  mainFrame.set_label("Main files");
  mainFrame.add(hboxFrame);
  
  buttonRun.set_label("Run!");
  buttonRun.signal_clicked().
    connect(sigc::mem_fun(*this, &NewAnalysis::runAnalysis));
  buttonBoxRun.set_layout(BUTTONBOX_END);
  buttonBoxRun.pack_end(buttonRun);

  vboxMain.set_homogeneous(false);
  vboxMain.set_spacing(10);
  vboxMain.set_border_width(10);
  vboxMain.pack_start(mainFrame);
  vboxMain.pack_start(buttonBoxRun, PACK_SHRINK);

  add(vboxMain);

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

  spinBegin.set_sensitive();
  spinEnd.set_sensitive();
  hScaleBegin.set_sensitive();
  hScaleEnd.set_sensitive();

  spinBegin.set_range(0, (tmpGromacsFrames - 1) * tmpGromacs->getTimeStep());
  spinBegin.set_value(0);
  spinBegin.set_increments(tmpGromacs->getTimeStep(),
                           (tmpGromacsFrames - 1) / 100);

  spinEnd.set_range(0, (tmpGromacsFrames - 1) * tmpGromacs->getTimeStep());
  spinEnd.set_value((tmpGromacsFrames - 1) * tmpGromacs->getTimeStep());
  spinEnd.set_increments(tmpGromacs->getTimeStep(),
                         (tmpGromacsFrames - 1) / 100);

  delete tmpGromacs;
}

void
NewAnalysis::runAnalysis()
{
  Gromacs::Gromacs gromacs( trjChooser.get_filename(),
                            tprChooser.get_filename());

  gromacs.setBegin(spinBegin.get_value());
  gromacs.setEnd(spinEnd.get_value());

  if(not progress.is_ancestor(vboxMain))
    vboxMain.pack_end(progress, PACK_EXPAND_WIDGET, 10);

  set_sensitive(false);

  progress.set_fraction(0);
  progress.show();
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
  progress.set_fraction(1);
  while(Main::events_pending())
    Main::iteration();

  Gromacs::Pittpi pittpi(gromacs, "/tmp/sas.csf", spinRadius.get_value(),
                         spinPocketThreshold.get_value());

  progress.hide();
  set_sensitive(true);
}

void
NewAnalysis::chooserTrajectoryClicked()
{
  if(not fs::exists(fs::path(trjChooser.get_filename())))
  {
    spinBegin.set_sensitive(false);
    spinEnd.set_sensitive(false);
    hScaleBegin.set_sensitive(false);
    hScaleEnd.set_sensitive(false);
  }
  else
    Glib::Thread::create(
      sigc::mem_fun(*this, &NewAnalysis::threadTrajectoryClicked), false);
}

void
NewAnalysis::threadTrajectoryClicked()
{
  signal_start_spin();

  tmpGromacs = new Gromacs::Gromacs(trjChooser.get_filename(), "");
  tmpGromacsFrames = tmpGromacs->getFramesCount();

  signal_stop_spin();
}

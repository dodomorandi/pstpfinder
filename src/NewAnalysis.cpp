#include <gtkmm.h>

#include "trapof.h"
#include "NewAnalysis.h"
#include "Gromacs.h"

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
  if(Glib::thread_supported())
    Glib::thread_init();

  FileFilter trjFilter;
  trjFilter.set_name("Trajectory files");
  trjFilter.add_pattern("*.xtc");
  trjFilter.add_pattern("*.trj");
  trjChooser.add_filter(trjFilter);
  trjChooser.signal_file_set().connect(
    sigc::mem_fun(*this, &NewAnalysis::chooserTrajectoryClicked));
  
  labelTrajectory.set_label("Trajectory file:");
  
  hboxTrajectory.pack_start(labelTrajectory, PACK_SHRINK);
  hboxTrajectory.pack_start(trjChooser);
  hboxTrajectory.set_homogeneous(false);
  hboxTrajectory.set_spacing(10);
  
  FileFilter tprFilter;
  tprFilter.set_name("Topology files");
  tprFilter.add_pattern("*.tpr");
  tprChooser.add_filter(tprFilter);
  
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

  vboxFrame.pack_start(hboxTrajectory, PACK_EXPAND_PADDING);
  vboxFrame.pack_start(hboxTopology, PACK_EXPAND_PADDING);
  vboxFrame.pack_start(hboxBegin, PACK_EXPAND_PADDING);
  vboxFrame.pack_start(hboxEnd, PACK_EXPAND_PADDING);
  vboxFrame.set_spacing(10);
  vboxFrame.set_border_width(10);
  
  mainFrame.set_label("Main files");
  mainFrame.add(vboxFrame);
  
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
  set_size_request(300);

  signal_delete_event().connect(sigc::ptr_fun(&closeApplication));
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
  progress.set_fraction(1);
  while(Main::events_pending())
    Main::iteration();

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
  if(spinnerWait.get_parent() == 0)
    vboxMain.pack_start(spinnerWait);
  mainFrame.hide();
  buttonBoxRun.hide();
  spinnerWait.show();
  spinnerWait.start();

  Gromacs::Gromacs gromacs(trjChooser.get_filename(),"");
  int frames = gromacs.getFramesCount();

  spinnerWait.stop();
  spinnerWait.hide();
  mainFrame.show();
  buttonBoxRun.show();

  spinBegin.set_sensitive();
  spinEnd.set_sensitive();
  hScaleBegin.set_sensitive();
  hScaleEnd.set_sensitive();

  spinBegin.set_range(0, (frames - 1) * gromacs.getTimeStep());
  spinBegin.set_value(0);
  spinBegin.set_increments(gromacs.getTimeStep(), (frames - 1) / 100);

  spinEnd.set_range(0, (frames - 1) * gromacs.getTimeStep());
  spinEnd.set_value((frames - 1) * gromacs.getTimeStep());
  spinEnd.set_increments(gromacs.getTimeStep(), (frames - 1) / 100);
}

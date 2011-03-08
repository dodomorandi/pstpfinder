#include <gtkmm.h>

#include "trapof.h"
#include "NewAnalysis.h"
#include "Gromacs.h"

using namespace Gtk;

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
  FileFilter trjFilter;
  trjFilter.set_name("Trajectory files");
  trjFilter.add_pattern("*.xtc");
  trjFilter.add_pattern("*.trj");
  trjChooser.add_filter(trjFilter);
  
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

  vboxFrame.pack_start(hboxTrajectory, PACK_EXPAND_PADDING);
  vboxFrame.pack_start(hboxTopology, PACK_EXPAND_PADDING);
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
  
  gromacs.calculateSas();
}
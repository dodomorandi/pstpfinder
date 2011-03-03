#include <gtkmm.h>

#include "NewAnalysis.h"

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

  vbox1.pack_start(hboxTrajectory, PACK_EXPAND_PADDING);
  vbox1.pack_start(hboxTopology, PACK_EXPAND_PADDING);
  vbox1.set_spacing(10);
  vbox1.set_border_width(10);
  
  mainFrame.set_label("Main files");
  mainFrame.set_border_width(10);
  mainFrame.add(vbox1);
  
  add(mainFrame);
  set_size_request(300);
}

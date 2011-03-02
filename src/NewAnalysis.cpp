#include <gtkmm.h>

#include "NewAnalysis.h"

using namespace Gtk;

NewAnalysis::NewAnalysis()
{
  FileFilter fileFilter;
  fileFilter.set_name("Trajectory files");
  fileFilter.add_pattern("*.xtc");
  fileFilter.add_pattern("*.trj");
  trjChooser.add_filter(fileFilter);
  
  vbox1.pack_start(trjChooser);
  
  mainFrame.set_label("Main files");
  mainFrame.set_border_width(10);
  mainFrame.add(vbox1);
  
  add(mainFrame);
  show_all();
}
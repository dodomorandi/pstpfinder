#ifndef _NEWANALYSIS_H
#define _NEWANALYSIS_H

#include <gtkmm.h>

using namespace Gtk;

class NewAnalysis: public Window
{
public:
  NewAnalysis();
  NewAnalysis(Window& parent);
private:
  Frame mainFrame;
  VBox vboxFrame, vboxMain;
  FileChooserButton trjChooser, tprChooser;
  Label labelTrajectory, labelTopology;
  HBox hboxTrajectory, hboxTopology;
  HButtonBox buttonBoxRun;
  Button  buttonRun;
  
  void init();
  void runAnalysis();
};

#endif
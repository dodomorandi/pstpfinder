#ifndef _NEWANALYSIS_H
#define _NEWANALYSIS_H

#include <gtkmm.h>

using namespace Gtk;

class NewAnalysis: public Window
{
public:
  enum parallelOperation
  {
    OPERATION_WAIT
  };

  NewAnalysis();
  NewAnalysis(Window& parent);
private:
  Frame mainFrame;
  VBox vboxFrame, vboxMain;
  FileChooserButton trjChooser, tprChooser;
  Label labelTrajectory, labelTopology, labelBegin, labelEnd;
  HBox hboxTrajectory, hboxTopology, hboxBegin, hboxEnd;
  HButtonBox buttonBoxRun;
  Button  buttonRun;
  ProgressBar progress;
  SpinButton spinBegin, spinEnd;
  HScale hScaleBegin, hScaleEnd;
  Spinner spinnerWait;

  void init();
  void runAnalysis();
  void chooserTrajectoryClicked();
  void threadTrajectoryClicked();
};

#endif

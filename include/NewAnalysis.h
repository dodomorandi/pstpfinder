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
  VBox vboxFrame1, vboxFrame2, vboxMain;
  FileChooserButton trjChooser, tprChooser;
  Label labelTrajectory, labelTopology, labelBegin, labelEnd, labelRadius,
        labelPocketThreshold, labelPs;
  HBox hboxTrajectory, hboxTopology, hboxBegin, hboxEnd, hboxFrame, hboxRadius,
       hboxPocketThreshold;
  HButtonBox buttonBoxRun;
  Button  buttonRun;
  ProgressBar progress;
  SpinButton spinBegin, spinEnd, spinRadius, spinPocketThreshold;
  HScale hScaleBegin, hScaleEnd;
  Spinner spinnerWait;
  VSeparator vSeparator;

  void init();
  void runAnalysis();
  void chooserTrajectoryClicked();
  void threadTrajectoryClicked();
};

#endif

#ifndef _NEWANALYSIS_H
#define _NEWANALYSIS_H

#include "Gromacs.h"

#include <gtkmm.h>
#include <glibmm.h>

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
  Glib::Dispatcher signal_start_spin;
  Glib::Dispatcher signal_stop_spin;
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
  Gromacs::Gromacs* tmpGromacs;
  int tmpGromacsFrames;

  void init();
  void runAnalysis();
  void chooserTrajectoryClicked();
  void threadTrajectoryClicked();
  void start_spin();
  void stop_spin();
};

#endif

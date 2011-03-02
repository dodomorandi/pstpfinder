#ifndef _MAINWINDOW_H
#define _MAINWINDOW_H

#include <gtkmm.h>

#include "NewAnalysis.h"

using namespace Gtk;

class MainWindow: public Window
{
public:
  MainWindow();
  ~MainWindow();
  void createNewAnalysis();
  void destroyNewAnalysis();
protected:
  HButtonBox buttonBox;
  Button buttonNew, buttonOpen;
private:
  void init();
  NewAnalysis* newAnalysis;
};

#endif
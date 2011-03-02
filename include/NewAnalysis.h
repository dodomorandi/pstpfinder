#ifndef _NEWANALYSIS_H
#define _NEWANALYSIS_H

#include <gtkmm.h>

using namespace Gtk;

class NewAnalysis: public Window
{
public:
  NewAnalysis();
private:
  Frame mainFrame;
  VBox vbox1;
  FileChooserButton trjChooser;
};

#endif
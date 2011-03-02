#include <gtkmm.h>
#include <vector>

#include "MainWindow.h"

using namespace Gtk;

MainWindow::MainWindow()
{
  this->init();
  this->show_all();
}

void
MainWindow::init()
{
  this->buttonNew.set_label("New analysis");
  this->buttonNew.signal_clicked().
    connect(sigc::mem_fun(*this, &MainWindow::createNewAnalysis));
  this->buttonOpen.set_label("Open analysis...");
  
  this->buttonBox.add((Widget&)this->buttonNew);
  this->buttonBox.add(this->buttonOpen);
  this->buttonBox.set_layout(BUTTONBOX_SPREAD);
  this->buttonBox.set_spacing(10);
  this->buttonBox.set_border_width(10);
  this->buttonBox.set_child_min_height(40);
  
  this->add(this->buttonBox);
}

void
MainWindow::createNewAnalysis()
{
  MessageDialog msg(*this,"prova");
  msg.run();
}
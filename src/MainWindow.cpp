#include <gtkmm.h>
#include <vector>

#include "MainWindow.h"
#include "NewAnalysis.h"

using namespace Gtk;

MainWindow::MainWindow()
{
  init();
  show_all();
}

MainWindow::~MainWindow()
{
  destroyNewAnalysis();
}

void
MainWindow::init()
{
  buttonNew.set_label("New analysis");
  buttonNew.signal_clicked().
    connect(sigc::mem_fun(*this, &MainWindow::createNewAnalysis));
  buttonOpen.set_label("Open analysis...");
  
  buttonBox.add((Widget&)buttonNew);
  buttonBox.add(buttonOpen);
  buttonBox.set_layout(BUTTONBOX_SPREAD);
  buttonBox.set_spacing(10);
  buttonBox.set_border_width(10);
  buttonBox.set_child_min_height(40);
  
  add(buttonBox);
}

void
MainWindow::createNewAnalysis()
{
  if(newAnalysis == 0)
  {
    newAnalysis = new NewAnalysis();
    newAnalysis->signal_unmap().
      connect(sigc::mem_fun(*this, &MainWindow::destroyNewAnalysis));
  }
}

void
MainWindow::destroyNewAnalysis()
{
  if(newAnalysis != 0)
  {
    delete(newAnalysis);
    newAnalysis = 0;
  }
}
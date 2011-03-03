#include <gtkmm.h>

#include "main.h"
#include "MainWindow.h"

using namespace Gtk;

static Main* kit;

int main(int argc, char* argv[])
{
  kit = new Main(argc, argv);

  MainWindow win;
  kit->run();
  
  delete kit;
  
  return 0;
}

bool close_app(GdkEventAny* event)
{
  kit->quit();
  return false;
}

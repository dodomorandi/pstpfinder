#include <gtkmm.h>
#include "MainWindow.h"

using namespace Gtk;

int main(int argc, char* argv[])
{
  Main kit(argc, argv);

  MainWindow win;
  Main::run(win);
}

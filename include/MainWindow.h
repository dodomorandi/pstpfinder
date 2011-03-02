#include <gtkmm.h>

using namespace Gtk;

class MainWindow: public Window
{
public:
  MainWindow();
  void createNewAnalysis();
protected:
  HButtonBox buttonBox;
  Button buttonNew, buttonOpen;
private:
  void init();
};
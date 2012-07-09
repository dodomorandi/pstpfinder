/*
 *  This file is part of PSTP-finder, an user friendly tool to analyze GROMACS
 *  molecular dynamics and find transient pockets on the surface of proteins.
 *  Copyright (C) 2011 Edoardo Morandi.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "pstpfinder.h"
#include "MainWindow.h"
#include "NewAnalysis.h"
#include "GtkmmWrapper.h"
#include "utils.h"

#include <gtkmm.h>
#include <vector>

using namespace Gtk;

namespace PstpFinder
{

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
    set_title("PSTP-finder");
    try
    {
      Glib::RefPtr<IconTheme> defaultIconTheme(IconTheme::get_default());
      auto searchPath(defaultIconTheme->get_search_path());
#if GTKMM_MAJOR == 3
      searchPath.push_back("./");
      searchPath.push_back("../");
#endif
      defaultIconTheme->set_search_path(searchPath);
      Glib::RefPtr<Gdk::Pixbuf> icon(
          defaultIconTheme->load_icon("pstpfinder", 64));
      set_default_icon(icon);
    }
    catch(IconThemeError& error)
    {
      cerr << "Warning: can't set PSTP-finder icon." << endl;
    }
    buttonNew.set_label("New analysis");
    buttonNew.signal_clicked().connect(
        sigc::mem_fun(*this, &MainWindow::createNewAnalysis));
    buttonOpen.set_label("Open analysis...");
    buttonOpen.signal_clicked().connect(
        sigc::mem_fun(*this, &MainWindow::buttonOpenClicked));

    buttonBox.add((Widget&) buttonNew);
    buttonBox.add(buttonOpen);
    buttonBox.set_layout(BUTTONBOX_SPREAD);
    buttonBox.set_spacing(10);
    buttonBox.set_border_width(10);
    //buttonBox.set_child_min_height(40);

    add(buttonBox);

    newAnalysis = 0;
    signal_delete_event().connect(sigc::ptr_fun(&closeApplication));
  }

  void
  MainWindow::createNewAnalysis()
  {
    if(newAnalysis == 0)
    {
      newAnalysis = new NewAnalysis(*this);
      newAnalysis->set_position(WIN_POS_CENTER_ON_PARENT);
      newAnalysis->signal_unmap().connect(
          sigc::mem_fun(*this, &MainWindow::destroyNewAnalysis));
      hide();
    }
  }

  void
  MainWindow::destroyNewAnalysis()
  {
    if(newAnalysis != 0)
    {
      newAnalysis->hide();
      delete (newAnalysis);
      newAnalysis = 0;
      if(not quitting)
        show();
    }
  }

  void
  MainWindow::buttonOpenClicked()
  {
    int response;
    string filename;
    {
      GtkmmWrapper<FileFilter> filter;
      filter->add_pattern("*.csf");
      filter->set_name("PSTP-filter compressed session file");

      FileChooserDialog chooser("Choose a saving file for this session",
                                FILE_CHOOSER_ACTION_OPEN);
      chooser.add_filter(filter);
      chooser.add_button(Stock::CANCEL, RESPONSE_CANCEL);
      chooser.add_button(Stock::OPEN, RESPONSE_OK);
      response = chooser.run();
      filename = chooser.get_filename();
    }

    switch(response)
    {
      case RESPONSE_OK:
        if(exists(filename) and newAnalysis == 0)
        {
          hide();
          newAnalysis = new NewAnalysis(*this);
          newAnalysis->set_position(WIN_POS_CENTER_ON_PARENT);
          newAnalysis->signal_unmap().connect(
              sigc::mem_fun(*this, &MainWindow::destroyNewAnalysis));
          newAnalysis->openSessionFile(filename);
        }
        break;
    }
  }

}

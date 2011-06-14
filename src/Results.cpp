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

#include "Results.h"
#include <gdkmm.h>
#include <cairomm/cairomm.h>

using namespace Gromacs;

Results::Results(NewAnalysis& parent, const Pittpi& pittpi)
  : pittpi(pittpi), parent(parent)
{
  init();
}

void
Results::init()
{
  signal_delete_event()
      .connect(sigc::mem_fun(*this, &Results::removeFromParent));

  drawResultsGraph.set_size_request(500, 200);
  drawResultsGraph.signal_expose_event()
    .connect(sigc::mem_fun(*this, &Results::drawResultsGraphExposeEvent));

  notebook.append_page(drawResultsGraph, "Results");

  add(notebook);
  show_all();
}

bool
Results::removeFromParent(GdkEventAny* event)
{
  parent.deleteResultsWindow(*this);
  return true;
}

bool
Results::drawResultsGraphExposeEvent(GdkEventExpose* event)
{
  Glib::RefPtr<Gdk::Window> window = drawResultsGraph.get_window();
  Cairo::RefPtr<Cairo::Context> context = window->create_cairo_context();

  window->clear();

  Gdk::Rectangle area_paint(0, 0, window->get_height(), window->get_width());
  window->begin_paint_rect(area_paint);
  // Color background
  context->set_source_rgb(1.0, 1.0, 1.0);
  context->paint();

  // Draw graph
  // TODO

  window->end_paint();

  return true;
}

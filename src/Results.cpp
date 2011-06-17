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
#include <cstdlib>
#include <cmath>
#include <vector>
#include <gdkmm.h>
#include <cairomm/cairomm.h>

using namespace std;
using namespace Gromacs;

Results::Results(NewAnalysis& parent, const Pittpi& pittpi,
                 const Gromacs& gromacs)
  : pittpi(pittpi, gromacs), gromacs(gromacs), parent(parent)
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

  fillResidues();
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

  Gdk::Rectangle area_paint(0, 0, window->get_width(), window->get_height());
  window->begin_paint_rect(area_paint);
  // Color background
  context->set_source_rgb(1.0, 1.0, 1.0);
  context->paint();

  // Draw graph
  // Axis
  context->set_source_rgb(0.0, 0.0, 0.0);
  context->set_line_width(1);
  context->move_to(20, area_paint.get_height() - 20);
  context->line_to(area_paint.get_width() - 20, area_paint.get_height() - 20);
  context->move_to(20, area_paint.get_height() - 20);
  context->line_to(20, 10);
  context->stroke();

  // Pockets
  int graphOffsetStart = 20 + 1;
  float columnModuleX = (float)(area_paint.get_width() - graphOffsetStart * 2)
                        / (residues.size() * 3 + 1);
  float columnModuleY = (float)(area_paint.get_height() - graphOffsetStart * 2)
                        / maxPocketLength;
  for
  (
    vector<PocketResidue>::const_iterator i = residues.begin();
    i < residues.end();
    i++
  )
  {
    int columnOffsetX = graphOffsetStart + columnModuleX *
        (3 * distance(static_cast<vector<PocketResidue>::const_iterator>
        (residues.begin()), i) + 1);
    int columnOffsetY = area_paint.get_height() - graphOffsetStart;

    for
    (
      vector<const Pocket*>::const_iterator j = i->pockets.begin();
      j < i->pockets.end();
      j++
    )
    {
      int columnHeight = (float)columnModuleY * (*j)->width;
      const Gdk::Color& color = colors[distance(i->pockets.begin(), j)];

      context->set_source_rgb(color.get_red_p(),
                              color.get_green_p(),
                              color.get_blue_p());
      context->rectangle(columnOffsetX, columnOffsetY - columnHeight,
                         columnModuleX * 2, columnHeight);
      context->fill();
      columnOffsetY -= columnHeight;
    }
  }

  window->end_paint();

  return true;
}

void
Results::fillResidues()
{
  const vector<Pocket>& pockets = pittpi.getPockets();
  maxPocketLength = 0;
  unsigned int maxPocketsPerResidue = 0;

  for
  (
    vector<Pocket>::const_iterator i = pockets.begin();
    i < pockets.end();
    i++
  )
  {
    bool exists = false;
    const Residue& currentRes = i->group->getCentralRes();

    for
    (
      vector<PocketResidue>::iterator j = residues.begin();
      j < residues.end();
      j++
    )
    {
      if(j->residue.index == currentRes.index)
      {
        j->pockets.push_back(&(*i));

        if(j->pockets.size() > maxPocketsPerResidue)
          maxPocketsPerResidue = j->pockets.size();

        exists = true;
        break;
      }
    }

    if(not exists)
    {
      PocketResidue pocket(currentRes);
      pocket.pockets.push_back(&(*i));

      residues.push_back(pocket);
      if(maxPocketsPerResidue == 0)
        maxPocketsPerResidue = 1;
    }

    if(i->width > maxPocketLength)
      maxPocketLength = i->width;
  }

  for(unsigned int i = 0; i < maxPocketsPerResidue; i++)
  {
    Gdk::Color color;
    color.set_rgb_p((float)rand() / RAND_MAX * 0.8,
                    (float)rand() / RAND_MAX * 0.8,
                    (float)rand() / RAND_MAX * 0.8);

    colors.push_back(color);
  }
}

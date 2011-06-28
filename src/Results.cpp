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
#include <algorithm>
#include <string>
#include <gdkmm.h>
#include <cairomm/cairomm.h>

using namespace std;
using namespace Gromacs;

Results::Results(NewAnalysis& parent, const Pittpi& pittpi,
                 const Gromacs& gromacs)
  : gromacs(gromacs), pittpi(pittpi, this->gromacs), parent(parent)
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
  int graphLineWidth = 1;
  int graphBorder = 10;
  int graphOffsetStart = graphBorder + graphLineWidth;
  int graphFooterHeight = 0.08 * area_paint.get_height(); // 8%
  context->select_font_face("Sans", Cairo::FONT_SLANT_NORMAL,
                            Cairo::FONT_WEIGHT_NORMAL);
  // Axis
  context->set_source_rgb(0.0, 0.0, 0.0);
  context->set_line_width(graphLineWidth);
  context->move_to(graphBorder,
                   area_paint.get_height() - graphBorder - graphFooterHeight);
  context->line_to(area_paint.get_width() - graphBorder,
                   area_paint.get_height() - graphBorder - graphFooterHeight);
  context->move_to(graphBorder,
                   area_paint.get_height() - graphBorder - graphFooterHeight);
  context->line_to(graphBorder, graphBorder);
  context->stroke();

  // Pockets
  float columnModuleX = (float)(area_paint.get_width() - graphOffsetStart * 2)
                        / (residues.size() * 3 + 1);
  float columnModuleY = (float)(area_paint.get_height() - graphOffsetStart * 2
                        - graphFooterHeight) / maxPocketLength;
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
    int columnOffsetY = area_paint.get_height() - graphOffsetStart
                        - graphFooterHeight;

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

    stringstream index;
    Cairo::TextExtents extents;
    index << i->residue->index;
    string strIndex = index.str();
    context->set_source_rgb(0, 0, 0);
    context->set_font_size(graphFooterHeight * 0.6);
    context->get_text_extents(strIndex, extents);
    if(extents.width > columnModuleX * 2)
    {
      Cairo::Matrix matrix;
      context->get_font_matrix(matrix);
      matrix.scale((float)columnModuleX * 2 / extents.width, 1);
      context->set_font_matrix(matrix);
      extents.width = columnModuleX * 2;
    }
    context->move_to(columnOffsetX + columnModuleX - extents.width / 2,
                     area_paint.get_height() - graphOffsetStart
                     - graphFooterHeight * 0.2 );
    context->show_text(strIndex);
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
      if(j->residue->index == currentRes.index)
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
  }

  for
  (
    vector<PocketResidue>::const_iterator i = residues.begin();
    i < residues.end();
    i++
  )
  {
    unsigned int residueLength = 0;
    for
    (
      vector<const Pocket*>::const_iterator j = i->pockets.begin();
      j < i->pockets.end();
      j++
    )
      residueLength += (*j)->width;

    if(maxPocketLength < residueLength)
      maxPocketLength = residueLength;
  }


  vector<Gdk::Color> unscrambledColors;
  for(unsigned int i = 0; i < maxPocketsPerResidue; i++)
    unscrambledColors.push_back(rainbow(1.0 / (maxPocketsPerResidue - 1) * i));

  vector<unsigned int> indexes;
  indexes.reserve(maxPocketsPerResidue);
  for(unsigned int i = 0; i < maxPocketsPerResidue; i++)
  {
    unsigned int index;
    bool gotIt;
    do
    {
      gotIt = false;
      index = (double)rand() / RAND_MAX * maxPocketsPerResidue;
      for
      (
        vector<unsigned int>::const_iterator j = indexes.begin();
        j < indexes.end();
        j++
      )
      {
        if(index == *j)
        {
          gotIt = true;
          break;
        }
      }
    } while(gotIt);

    indexes.push_back(index);
    colors.push_back(unscrambledColors[index]);
  }

  sort(residues.begin(), residues.end(), PocketResidue::sortByResidueIndex);
}

Gdk::Color
Results::rainbow(double value)
{
  Gdk::Color color;
  double red, green, blue;

  if(value <= 0.25)
  {
    red = 1.0;
    green = value * 4;
    blue = 0;
  }
  else if(value <= 0.5)
  {
    red = 3.0 - value * 4;
    green = 1.0;
    blue = 0;
  }
  else if(value <= 0.75)
  {
    red = 0;
    green = 1.0;
    blue = (value - 0.5) * 4;
  }
  else
  {
    red = 0;
    green = 4.0 - value * 4;
    blue = 1.0;
  }

  color.set_rgb_p(red, green, blue);
  return color;
}

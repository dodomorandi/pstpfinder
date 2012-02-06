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
#include <sstream>
#include <iomanip>
#include <gdkmm.h>
#include <cairomm/cairomm.h>
#include <pangomm.h>

using namespace std;
using namespace PstpFinder;

Results::Results(NewAnalysis& parent, const shared_ptr<Pittpi>& pittpi,
                 const Gromacs& gromacs) :
    gromacs(gromacs),
    pittpi(pittpi),
    parent(parent),
    graphLineWidth(1),
    graphBorder(10),
    graphOffsetStart(graphBorder + graphLineWidth)
{
  labelYMultiplier = 1.;
  labelXMultiplier = 0.1;
  graphModifier = enumModifier::NOTHING;
  graphLeftBorder = 0;

  init();
}

void
Results::init() throw()
{
  fillResidues();
  signal_delete_event()
      .connect(sigc::mem_fun(*this, &Results::removeFromParent));

  drawResultsGraph.set_size_request(500, 200);
#if GTKMM_MAJOR == 3
  drawResultsGraph.signal_draw().connect(
      sigc::mem_fun(*this, &Results::drawResultsGraphExposeEvent));
#else
  drawResultsGraph.signal_expose_event()
    .connect(sigc::mem_fun(*this, &Results::drawResultsGraphExposeEvent));
#endif

  drawResultsGraph.signal_scroll_event().connect(
      sigc::mem_fun(*this, &Results::drawResultsGraphScrollEvent));
  drawResultsGraph.signal_motion_notify_event().connect(
      sigc::mem_fun(*this, &Results::drawResultsGraphMotionEvent));
  drawResultsGraph.add_events(
      Gdk::EventMask::SCROLL_MASK | Gdk::EventMask::POINTER_MOTION_MASK);

  drawResultsStatusBar.push("Move the pointer over graph bars to get more information");

  drawResultsVBox.pack_start(drawResultsGraph);
  drawResultsVBox.pack_start(drawResultsStatusBar, false, true);
  notebook.append_page(drawResultsVBox, "Results");

  stringstream streamData, streamDetails;
  streamData << setfill(' ') << setw(11) << left << "zeros";
  streamData << setfill(' ') << setw(11) << left << "center";
  streamData << setfill(' ') << setw(11) << left << "start";
  streamData << setfill(' ') << setw(11) << left << "end";
  streamData << setfill(' ') << setw(11) << left << "duration";
  streamData << "group members" << endl;

  for(auto i = residues.cbegin(); i < residues.cend(); i++)
  {
    vector<const Pocket*>::const_iterator bestPocket = i->pockets.end();
    bool writtenHeader = false;

    for
    (
      vector<const Pocket*>::const_iterator j = i->pockets.begin();
      j < i->pockets.end();
      j++
    )
    {
      if(not writtenHeader)
      {
        streamDetails << "Pocket centered on "
                        << aminoacidTriplet[i->residue->type]
                        << i->residue->index << ":";

        const vector<const Residue*>& residuesRef = (*j)->group->getResidues();
        for
        (
          vector<const Residue*>::const_iterator k = residuesRef.begin();
          k < residuesRef.end();
          k++
        )
          streamDetails << " " << (*k)->index;
        streamDetails << endl << endl;

        streamDetails << setfill(' ') << setw(11) << left << "start";
        streamDetails << setfill(' ') << setw(11) << left << "duration";
        streamDetails << setfill(' ') << setw(12) << left << "ps max area";
        streamDetails << setfill(' ') << setw(11) << left << "ps average";
        streamDetails << "percentage" << endl;

        writtenHeader = true;
      }

      if(bestPocket == i->pockets.end() or (*bestPocket)->width < (*j)->width)
        bestPocket = j;

      stringstream perc;
      perc << static_cast<int>((*j)->openingFraction * 100) << "%";
      streamDetails << setfill('0') << setw(7) << right << (int)(*j)->startPs
                      << "    ";
      streamDetails << setfill('0') << setw(7) << right << (int)(*j)->width
                      << "    ";
      streamDetails << setfill('0') << setw(7) << right << (int)(*j)->maxAreaPs
                      << "     ";
      streamDetails << setfill('0') << setw(7) << right
                      << (int)(*j)->averageNearPs << "    ";
      streamDetails << perc.str();

      streamDetails << endl;
    }

    if(bestPocket != i->pockets.end())
    {
      const Residue& centralRes = (*bestPocket)->group->getCentralRes();
      stringstream aaRef;
      aaRef << aminoacidTriplet[centralRes.type] << centralRes.index;

      streamData << setfill('0') << setw(7) << right << (*bestPocket)->group->zeros
                                << "    ";
      streamData << setfill(' ') << setw(11) << left << aaRef.str();
      streamData << setfill('0') << setw(7) << right << (int)(*bestPocket)->startPs
                << "    ";
      streamData << setfill('0') << setw(7) << right << (int)(*bestPocket)->endPs
                << "    ";
      streamData << setfill('0') << setw(7) << right << (int)(*bestPocket)->width
                << "   ";

      const vector<const Residue*>& pocketResidues =
        (*bestPocket)->group->getResidues();
      for
      (
        vector<const Residue*>::const_iterator k = pocketResidues.begin();
        k < pocketResidues.end();
        k++
      )
        streamData << " " << (*k)->index;

      streamData << endl;
    }

    if(writtenHeader)
      streamDetails << endl;
  }

  Pango::FontDescription fontMono;
  fontMono.set_family("mono");
  fontMono.set_size(8 * Pango::SCALE);

  textViewData.get_buffer()->set_text(streamData.str());
  textViewData.set_editable(false);
#if GTKMM_MAJOR == 3
  textViewData.override_font(fontMono);
#else
  textViewData.modify_font(fontMono);
#endif
  scrollData.add(textViewData);
  notebook.append_page(scrollData, "Data");

  textViewDetails.get_buffer()->set_text(streamDetails.str());
  textViewDetails.set_editable(false);
#if GTKMM_MAJOR == 3
  textViewDetails.override_font(fontMono);
#else
  textViewDetails.modify_font(fontMono);
#endif
  scrollDetails.add(textViewDetails);
  notebook.append_page(scrollDetails, "Details");

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
#if GTKMM_MAJOR == 3
Results::drawResultsGraphExposeEvent(
    const Cairo::RefPtr<Cairo::Context>& event) throw ()
#else
Results::drawResultsGraphExposeEvent(GdkEventExpose* event) throw()
#endif
{
  Glib::RefPtr<Gdk::Window> window = drawResultsGraph.get_window();
  Cairo::RefPtr<Cairo::Context> context = window->create_cairo_context();

#if GTKMM_MAJOR == 2
  window->clear();
#endif

  Gdk::Rectangle area_paint(0, 0, window->get_width(), window->get_height());
  window->begin_paint_rect(area_paint);
  // Color background
  context->set_source_rgb(1.0, 1.0, 1.0);
  context->paint();

  // Draw graph
  graphFooterHeight = area_paint.get_height() * labelXMultiplier;
  if(area_paint.get_height() - graphFooterHeight / 0.8 - graphBorder
     - graphOffsetStart
     < 0)
  {
    graphFooterHeight = (area_paint.get_height() - graphOffsetStart
                         - graphBorder) * 0.8;
    labelXMultiplier = static_cast<float>(graphFooterHeight)
                       / area_paint.get_height();
  }
  else if(graphFooterHeight < 8)
  {
    graphFooterHeight = 8;
    labelXMultiplier = static_cast<float>(graphFooterHeight)
                       / area_paint.get_height();
  }

  float graphLabelYSize;
  context->set_identity_matrix();
  context->save();
  {
    Cairo::TextExtents extents;
    context->rotate(-3.1415 / 2);
    // We use an approximate value to evaluate graphHeaderHeight
    context->get_text_extents("pocket opening time(ps)", extents);
    graphLabelYSize = 10. / extents.width
                      * (area_paint.get_height() - graphFooterHeight
                         - graphBorder
                         - graphOffsetStart);

    context->set_font_size(graphLabelYSize * labelYMultiplier);
    context->get_text_extents("000", extents);
    if(graphLeftBorder != 0
       and (extents.height > graphLeftBorder * 0.4 or extents.height < 5))
    {
      float fontSize;
      if(extents.height < 5.)
        fontSize = 5.;
      else
        fontSize = graphLeftBorder * 0.4 * graphLabelYSize
                       * labelYMultiplier
                       / extents.height;
      labelYMultiplier = fontSize / graphLabelYSize;
      context->set_font_size(fontSize);
      context->get_text_extents("000", extents);
    }

    float numberWidth = extents.width;
    context->get_text_extents("00", extents);
    numberWidth -= extents.width;
    graphHeaderHeight = numberWidth * (ceil(log10(maxPocketLength + 1)) / 2);

    // Now we can calculate the true values
    context->set_font_size(graphLabelYSize);
    context->get_text_extents("pocket opening time(ps)", extents);
    graphLabelYSize = graphLabelYSize / extents.width
                      * (area_paint.get_height() - graphFooterHeight
                         - graphBorder
                         - graphOffsetStart
                         - graphHeaderHeight) * 0.9;
    graphLeftBorder = extents.height * 2.5;
    context->set_font_size(graphLabelYSize * labelYMultiplier);
    context->get_text_extents("0", extents);
    graphHeaderHeight = extents.width * (ceil(log10(maxPocketLength)) / 2 + 1);
  }
  context->restore();

  context->select_font_face("Sans", Cairo::FONT_SLANT_NORMAL,
                            Cairo::FONT_WEIGHT_NORMAL);
  // Axis
  context->set_source_rgb(0.0, 0.0, 0.0);
  context->set_line_width(graphLineWidth);
  context->move_to(graphBorder + graphLeftBorder,
                   area_paint.get_height() - graphBorder - graphFooterHeight);
  context->line_to(area_paint.get_width() - graphBorder,
                   area_paint.get_height() - graphBorder - graphFooterHeight);
  context->move_to(graphBorder + graphLeftBorder ,
                   area_paint.get_height() - graphBorder - graphFooterHeight);
  context->line_to(graphBorder + graphLeftBorder,
                   graphBorder + graphHeaderHeight);
  context->stroke();

  // Pockets
  float columnModuleX = (float)(area_paint.get_width() - graphOffsetStart * 2 -
                        graphLeftBorder) / (residues.size() * 3 + 1);
  float columnModuleY = (float)(area_paint.get_height() - graphOffsetStart * 2
                                - graphHeaderHeight - graphFooterHeight)
                        / maxPocketLength;
  for(auto i = residues.cbegin(); i < residues.cend(); i++)
  {
    int columnOffsetX = graphOffsetStart + graphLeftBorder + columnModuleX *
        (3 * distance(residues.cbegin(), i) + 1);
    int columnOffsetY = area_paint.get_height() - graphOffsetStart
                        - graphFooterHeight;

    for(auto j = i->pockets.cbegin(); j < i->pockets.cend(); j++)
    {
      int columnHeight = (float)columnModuleY * (*j)->width;
      const Gdk::Color& color = colors[distance(i->pockets.cbegin(), j)];

      context->set_source_rgb(color.get_red_p(),
                              color.get_green_p(),
                              color.get_blue_p());
      context->rectangle(columnOffsetX, columnOffsetY - columnHeight,
                         columnModuleX * 2, columnHeight);
      context->fill();
      columnOffsetY -= columnHeight;
    }

    // X Axis text
    stringstream index;
    Cairo::TextExtents extents;
    index << i->residue->index;
    string strIndex = index.str();
    if(graphModifier == enumModifier::LABEL_X)
      context->set_source_rgb(1, 0, 0);
    else
      context->set_source_rgb(0, 0, 0);
    context->set_font_size(graphFooterHeight);
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
                     area_paint.get_height() - graphOffsetStart);
    context->show_text(strIndex);
  }

  // Y Axis text
  context->save();
  {
    Cairo::TextExtents extents;
    context->set_identity_matrix();
    context->rotate(-3.1415 / 2);
    context->set_font_size(graphLabelYSize);
    if(graphModifier == enumModifier::LABEL_Y)
      context->set_source_rgb(1, 0, 0);
    else
      context->set_source_rgb(0, 0, 0);
    context->get_text_extents("pocket opening time(ps)", extents);

    unsigned int spaceBeforeLabelY = (area_paint.get_height()
                                      - graphHeaderHeight
                                      - graphFooterHeight
                                      - graphOffsetStart
                                      - graphBorder
                                      - extents.width)
                                     / 2;

    context->move_to(
        -graphHeaderHeight - graphBorder - extents.width - spaceBeforeLabelY,
        graphBorder + extents.height);
    context->show_text("pocket opening time(ps)");

    context->set_font_size(graphLabelYSize * labelYMultiplier);
    context->get_text_extents("000", extents);
    float numberWidth = extents.width;
    context->get_text_extents("00", extents);
    numberWidth -= extents.width;
    context->move_to(
        -area_paint.get_height() + graphOffsetStart + graphFooterHeight
        - numberWidth / 2,
        graphBorder + graphLeftBorder * 0.8);
    context->show_text("0");


    unsigned int graphMaxVerticalStep = (float)(area_paint.get_height()
      - graphOffsetStart * 2 - graphHeaderHeight - graphFooterHeight)
      / (numberWidth * (log10(maxPocketLength + 1) + 2));
    for(unsigned int i = 1; i <= graphMaxVerticalStep; i++)
    {
      stringstream pocketSizeStream;
      string pocketSize;
      pocketSizeStream << (int)(maxPocketLength / graphMaxVerticalStep * i);
      pocketSize = pocketSizeStream.str();

      context->get_text_extents(pocketSize, extents);
      context->move_to(- (area_paint.get_height() - graphOffsetStart
                          - graphFooterHeight)
                       + (float)(area_paint.get_height() - graphOffsetStart * 2
                                 - graphHeaderHeight - graphFooterHeight)
                       / graphMaxVerticalStep * i - extents.width / 2,
                       graphBorder + graphLeftBorder * 0.8);
      context->show_text(pocketSize);
    }

  }
  context->restore();

  window->end_paint();

  return true;
}

bool
Results::drawResultsGraphScrollEvent(GdkEventScroll* event) throw ()
{
  float* multiplier;

  if(graphModifier == enumModifier::LABEL_Y)
    multiplier = &labelYMultiplier;
  else if(graphModifier == enumModifier::LABEL_X)
    multiplier = &labelXMultiplier;
  else
    return true;

  if(event->direction == GdkScrollDirection::GDK_SCROLL_UP)
    *multiplier *= 1.2;
  else if(event->direction == GdkScrollDirection::GDK_SCROLL_DOWN)
    *multiplier /= 1.2;

  drawResultsGraph.queue_draw();
  return true;
}

bool
Results::drawResultsGraphMotionEvent(GdkEventMotion* event) throw()
{
  enumModifier oldModifier = graphModifier;

  if(event->x >= graphBorder and event->x < graphBorder + graphLeftBorder
     and event->y >= graphHeaderHeight + graphBorder
     and event->y
         < drawResultsGraph.get_height() - graphOffsetStart - graphFooterHeight)
    graphModifier = enumModifier::LABEL_Y;
  else if(event->y
          >= drawResultsGraph.get_height() - graphBorder - graphFooterHeight
          and event->y < drawResultsGraph.get_height() - graphBorder
          and event->x >= graphBorder + graphLeftBorder
          and event->x < drawResultsGraph.get_width() - graphBorder)
    graphModifier = enumModifier::LABEL_X;
  else
    graphModifier = enumModifier::NOTHING;

  if(graphModifier != oldModifier)
    drawResultsGraph.queue_draw();
  return true;
}

void
Results::fillResidues()
{
  const vector<Pocket>& pockets = pittpi->getPockets();
  maxPocketLength = 0;
  unsigned int maxPocketsPerResidue = 0;

  for(auto& pocket : pockets)
  {
    bool exists = false;
    const Residue& currentRes = pocket.group->getCentralRes();

    for(auto& residue : residues)
    {
      if(residue.residue->index == currentRes.index)
      {
        residue.pockets.push_back(&pocket);

        if(residue.pockets.size() > maxPocketsPerResidue)
          maxPocketsPerResidue = residue.pockets.size();

        exists = true;
        break;
      }
    }

    if(not exists)
    {
      PocketResidue pocketResidue(currentRes);
      pocketResidue.pockets.push_back(&pocket);

      residues.push_back(std::move(pocketResidue));
      if(maxPocketsPerResidue == 0)
        maxPocketsPerResidue = 1;
    }
  }

  for(auto& residue : residues)
  {
    unsigned int residueLength = 0;
    for(auto& pocket : residue.pockets)
      residueLength += pocket->width;

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
      index = (float)rand() / RAND_MAX * maxPocketsPerResidue;
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

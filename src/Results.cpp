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
#include "ColorsChooser.h"
#include <cstdlib>
#include <cmath>
#include <vector>
#include <utility>
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
    statusBarMessages
      { "Move the pointer over graph bars or axis labels to get more information",
        "Scroll mouse wheel to change font size",
        "Click to fix selection on this pocket",
        "Click to release selection",
        "Click to change selection to this pocket",
        "Move the pointer over axis labels or click on the graph to release selection"},
    gromacs(gromacs),
    pittpi(pittpi),
    parent(parent),
    graphLineWidth(1),
    graphBorder(10),
    graphOffsetStart(graphBorder + graphLineWidth)
{
  labelYMultiplier = 0.5;
  labelXMultiplier = 0.5;
  graphModifier = enumModifier::NOTHING;
  graphLeftBorder = 0;
  fixedSelection = false;

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
  drawResultsGraph.signal_draw().
      connect(sigc::mem_fun(*this, &Results::drawResultsGraphDrawEvent));
#else
  drawResultsGraph.signal_expose_event().
      connect(sigc::mem_fun(*this, &Results::drawResultsGraphExposeEvent));
#endif

  drawResultsGraph.signal_scroll_event().connect(
      sigc::mem_fun(*this, &Results::drawResultsGraphScrollEvent));
  drawResultsGraph.signal_motion_notify_event().connect(
      sigc::mem_fun(*this, &Results::drawResultsGraphMotionEvent));
  drawResultsGraph.signal_button_press_event().connect(
      sigc::mem_fun(*this, &Results::drawResultsGraphButtonPressEvent));
  drawResultsGraph.add_events(
      Gdk::EventMask::SCROLL_MASK | Gdk::EventMask::POINTER_MOTION_MASK
      | Gdk::EventMask::BUTTON_PRESS_MASK);

  labelPocketCenter.set_text("Pocket centered on");
  labelPocketStart.set_text("Pocket starts at ps");
  labelPocketEnd.set_text("Pocket ends at ps");
  labelPocketWidth.set_text("Pocket length in ps");
  entryPocketCenter.set_width_chars(8);
  entryPocketCenter.set_editable(false);
  entryPocketStart.set_width_chars(8);
  entryPocketStart.set_editable(false);
  entryPocketEnd.set_width_chars(8);
  entryPocketEnd.set_editable(false);
  entryPocketWidth.set_width_chars(8);
  entryPocketWidth.set_editable(false);

  hboxPocketCenter.set_spacing(5);
  hboxPocketCenter.pack_start(labelPocketCenter);
  hboxPocketCenter.pack_start(entryPocketCenter, false, false, 0);
  hboxPocketStart.set_spacing(5);
  hboxPocketStart.pack_start(labelPocketStart);
  hboxPocketStart.pack_start(entryPocketStart, false, false, 0);
  hboxPocketEnd.set_spacing(5);
  hboxPocketEnd.pack_start(labelPocketEnd);
  hboxPocketEnd.pack_start(entryPocketEnd, false, false, 0);
  hboxPocketWidth.set_spacing(5);
  hboxPocketWidth.pack_start(labelPocketWidth);
  hboxPocketWidth.pack_start(entryPocketWidth, false, false, 0);

  labelPocketResidues.set_text("Involved residues:");
  labelPocketResidues.set_alignment(0, 0.5);
  textPocketResidues.set_wrap_mode(Gtk::WrapMode::WRAP_WORD);
  textPocketResidues.set_editable(false);
  textPocketResidues.set_border_width(1);
#if GTKMM_MAJOR == 3
  textPocketResidues.set_hexpand(false);
#endif

  vboxPocketInformation.pack_start(hboxPocketCenter, false, false, 0);
  vboxPocketInformation.pack_start(hboxPocketStart, false, false, 0);
  vboxPocketInformation.pack_start(hboxPocketEnd, false, false, 0);
  vboxPocketInformation.pack_start(hboxPocketWidth, false, false, 0);
  vboxPocketInformation.pack_start(labelPocketResidues, false, false, 0);
  vboxPocketInformation.pack_start(textPocketResidues, false, false, 0);
  vboxPocketInformation.set_border_width(5);

  panedMain.pack1(drawResultsGraph, true, false);
  panedMain.pack2(vboxPocketInformation, true, true);
  drawResultsStatusBar.push(statusBarMessages[0]);
  drawResultsVBox.pack_start(panedMain);
  drawResultsVBox.pack_start(drawResultsStatusBar, false, true);
  notebook.append_page(drawResultsVBox, "Results");

  stringstream streamData, streamDetails;
  streamData << setfill(' ') << setw(11) << left << "zeros";
  streamData << setfill(' ') << setw(11) << left << "center";
  streamData << setfill(' ') << setw(11) << left << "start";
  streamData << setfill(' ') << setw(11) << left << "end";
  streamData << setfill(' ') << setw(11) << left << "duration";
  streamData << "group members" << endl;

  selectedPocket = nullptr;
  hoveringOnPocket = nullptr;
  for(auto i = begin(residues); i != end(residues); i++)
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

        auto& residuesRef = (*j)->group->getResidues();
        for(auto& residueRef : residuesRef)
          streamDetails << " " << residueRef->index;
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
      const Residue<SasPdbAtom>& centralRes = (*bestPocket)->group
          ->getCentralRes();
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

      for(auto& pocketResidue : (*bestPocket)->group->getResidues())
        streamData << " " << pocketResidue->index;

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

  buildMenu();

  vboxMain.pack_start(*uiManager->get_widget("/menuBar"), false, false);
  vboxMain.pack_start(notebook);
  add(vboxMain);
  set_default_size(800, 320);
  set_title("PSTP-finder results");

  show_all();
}

bool
Results::removeFromParent(GdkEventAny* event)
{
  (void) event;
  parent.deleteResultsWindow(*this);
  return true;
}

bool
#if GTKMM_MAJOR == 3
Results::drawResultsGraphDrawEvent(
    const Cairo::RefPtr<Cairo::Context>& context) throw ()
#else
Results::drawResultsGraphExposeEvent(GdkEventExpose*) throw()
#endif
{
#if GTKMM_MAJOR == 2
  Glib::RefPtr<Gdk::Window> window(drawResultsGraph.get_window());
  window->clear();
  Cairo::RefPtr<Cairo::Context> context(window->create_cairo_context());
#endif

  int height(drawResultsGraph.get_height());
  int width(drawResultsGraph.get_width());

#if GTKMM_MAJOR == 2
  Gdk::Rectangle area_paint(0, 0, width, height);
  window->begin_paint_rect(area_paint);
#endif

  // Color background
  context->set_source_rgb(1.0, 1.0, 1.0);
  context->paint();

  context->select_font_face("Sans", Cairo::FONT_SLANT_NORMAL,
                            Cairo::FONT_WEIGHT_NORMAL);

  float graphLabelXSize, graphLabelYSize;
  context->save();
  {
    Cairo::TextExtents extents;
    context->get_text_extents("Residue number", extents);
    float labelXProportion = extents.width / extents.height;
    graphLabelXSize = width / 3. / labelXProportion;

    if(labelXMultiplier != 0)
    {
      context->set_font_size(graphLabelXSize * labelXMultiplier);
      context->get_text_extents("000", extents);
      graphLabelXSize = extents.height / labelXMultiplier * 2;
    }
    else
    {
      context->set_font_size(graphLabelXSize);
      context->get_text_extents("Residue number", extents);
      graphLabelXSize = extents.height * 2;
    }

    context->set_font_size(graphLabelXSize * labelXMultiplier);
    context->get_text_extents("000", extents);
    graphBottomBorder = extents.height + 3;

    context->set_font_size(graphLabelXSize * (1. - labelXMultiplier));
    context->get_text_extents("Residue number", extents);
    graphBottomBorder += extents.height + 5;
  }
  context->restore();

  context->save();
  const float graphHeight = height - 2. * graphBorder - graphOffsetStart
                            - graphBottomBorder;
  {
    Cairo::TextExtents extents;
    context->rotate(-3.1415 / 2);
    context->get_text_extents("pocket disclosure time(ps)", extents);
    float labelYProportion = extents.width / extents.height;

    if(labelYMultiplier != 0)
    {
      graphLabelYSize = graphHeight / labelYProportion;
      context->set_font_size(graphLabelYSize * labelYMultiplier);
      context->get_text_extents("0", extents);

      graphLabelYSize = extents.height / labelYMultiplier * 2;
    }

    context->set_font_size(graphLabelYSize * (1. - labelYMultiplier));
    context->get_text_extents("pocket disclosure time(ps)", extents);
    if(labelYMultiplier == 0 or extents.width > graphHeight)
      graphLabelYSize = graphHeight / labelYProportion;

    context->set_font_size(graphLabelYSize * (1. - labelYMultiplier));
    context->get_text_extents("pocket disclosure time(ps)", extents);

    context->set_font_size(graphLabelYSize * labelYMultiplier);
    context->get_text_extents("000", extents);
    graphLeftBorder = extents.height + 3;

    context->set_font_size(graphLabelYSize * (1. - labelYMultiplier));
    context->get_text_extents("pocket disclosure time(ps)", extents);
    graphLeftBorder += extents.height + 5;
  }
  context->restore();

  // Axis
  context->set_source_rgb(0.0, 0.0, 0.0);
  context->set_line_width(graphLineWidth);
  context->move_to(graphBorder + graphLeftBorder,
                   height - graphBorder - graphBottomBorder);
  context->line_to(width - graphBorder,
                   height - graphBorder - graphBottomBorder);
  context->move_to(graphBorder + graphLeftBorder ,
                   height - graphBorder - graphBottomBorder);
  context->line_to(graphBorder + graphLeftBorder,
                   graphBorder);
  context->stroke();

  // Pockets
  float columnModuleX = (float)(width - graphOffsetStart * 2 -
                        graphLeftBorder) / (residues.size() * 3 + 1);
  float columnModuleY = (float)(height - graphOffsetStart * 2
                                - graphBottomBorder)
                        / maxPocketLength;

  bool foundSelection(false);
  array<float, 4> selectionCoordinates;

  float numbersHeight;
  {
    Cairo::TextExtents extents;
    context->set_font_size(graphLabelXSize * labelXMultiplier);
    context->get_text_extents("0123456789", extents);
    numbersHeight = extents.height;
  }

  for(auto i = begin(residues); i != end(residues); i++)
  {
    int columnOffsetX = graphOffsetStart + graphLeftBorder + columnModuleX *
        (3 * distance(begin(residues), i) + 1);
    int columnOffsetY = height - graphOffsetStart
                        - graphBottomBorder;

    for(auto j = begin(i->pockets); j != end(i->pockets); j++)
    {
      int columnHeight = (float)columnModuleY * (*j)->width;
      const Color& color = colors[distance(begin(i->pockets), j)];

      if((not fixedSelection and *j == hoveringOnPocket)
         or (fixedSelection and *j == selectedPocket))
      {
        selectionCoordinates = array<float, 4>
          {{ float(columnOffsetX), float(columnOffsetY - columnHeight),
             columnModuleX * 2, float(columnHeight) }};
        foundSelection = true;
      }

#if GTKMM_MAJOR == 2
      context->set_source_rgb(color.get_red_p(),
                              color.get_green_p(),
                              color.get_blue_p());
#else
      context->set_source_rgb(color.get_red(),
                              color.get_green(),
                              color.get_blue());
#endif
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
    context->set_font_size(graphLabelXSize * labelXMultiplier);
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
                     height - graphBorder - graphBottomBorder +
                     numbersHeight + 3);
    context->show_text(strIndex);
  }

  if(foundSelection)
  {
    context->rectangle(selectionCoordinates[0], selectionCoordinates[1],
                       selectionCoordinates[2], selectionCoordinates[3]);
    context->set_source_rgb(0, 0, 0);
    context->stroke();
  }

  // X Axis label text
  context->save();
  {
    Cairo::TextExtents extents;
    context->set_font_size(graphLabelXSize * (1. - labelXMultiplier));
    if(graphModifier == enumModifier::LABEL_X)
      context->set_source_rgb(1, 0, 0);
    else
      context->set_source_rgb(0, 0, 0);

    context->get_text_extents("Residue number", extents);

    int spaceBeforeLabelX = (width - graphBorder - graphOffsetStart
                             - extents.width)
                            / 2;
    context->move_to(graphOffsetStart + spaceBeforeLabelX,
                     height - graphBorder);
    context->show_text("Residue number");
  }
  context->restore();

  // Y Axis text
  context->save();
  {
    Cairo::TextExtents extents;
    context->rotate(-3.1415 / 2);
    context->set_font_size(graphLabelYSize * (1. - labelYMultiplier));
    if(graphModifier == enumModifier::LABEL_Y)
      context->set_source_rgb(1, 0, 0);
    else
      context->set_source_rgb(0, 0, 0);
    context->get_text_extents("pocket disclosure time(ps)", extents);

    int spaceBeforeLabelY = (height - graphBottomBorder
                             - graphOffsetStart - 2 * graphBorder - extents.width)
                            / 2;
    if(spaceBeforeLabelY < 0)
        spaceBeforeLabelY = 0;

    context->move_to(
        - graphBorder - extents.width - spaceBeforeLabelY,
        graphBorder + extents.height);
    context->show_text("pocket disclosure time(ps)");

    context->set_font_size(graphLabelYSize * labelYMultiplier);
    context->get_text_extents("000", extents);
    float numberWidth = extents.width;
    context->get_text_extents("00", extents);
    numberWidth -= extents.width;
    context->move_to(
        -height + graphOffsetStart + graphBottomBorder
        - numberWidth / 2,
        graphBorder + graphLeftBorder - 3);
    context->show_text("0");

    unsigned int lastPocket = (graphHeight -
        numberWidth * log10(maxPocketLength + 1) / 2.) *
        maxPocketLength / graphHeight;
    unsigned graphHeightWithoutLastPocket = graphHeight -
      numberWidth * log10(lastPocket + 1) / 2.;

    unsigned int graphMaxVerticalStep = graphHeightWithoutLastPocket 
      / (numberWidth * (log10(lastPocket + 1) + 2));
    for(unsigned int i = 1; i <= graphMaxVerticalStep; i++)
    {
      stringstream pocketSizeStream;
      string pocketSize;
      pocketSizeStream << static_cast<int>(lastPocket /
          graphMaxVerticalStep * i);
      pocketSize = pocketSizeStream.str();

      context->get_text_extents(pocketSize, extents);
      context->move_to(- (height - graphOffsetStart
                          - graphBottomBorder)
                       + static_cast<float>(graphHeightWithoutLastPocket)
                       / graphMaxVerticalStep * i -
                       extents.width / 2,
                       graphBorder + graphLeftBorder - 3);
      context->show_text(pocketSize);
    }

  }
  context->restore();

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
  {
    if(*multiplier == 0)
      *multiplier = 0.085;
    else
      *multiplier *= 1.2;
  }
  else if(event->direction == GdkScrollDirection::GDK_SCROLL_DOWN)
    *multiplier /= 1.2;

  if(*multiplier > 1.)
    *multiplier = 1.;
  else if(*multiplier < 0.08)
    *multiplier = 0;
  drawResultsGraph.queue_draw();
  return true;
}

bool
Results::drawResultsGraphMotionEvent(GdkEventMotion* event) throw()
{
  enumModifier oldModifier = graphModifier;
  double cursorX = event->x;
  double cursorY = event->y;

  int width(drawResultsGraph.get_width());
  int height(drawResultsGraph.get_height());
  bool gotcha(false);
  bool newSelection(false);

  if(event->x >= graphBorder and event->x < graphBorder + graphLeftBorder
     and event->y >= graphBorder
     and event->y
         < height - graphOffsetStart - graphBottomBorder)
    graphModifier = enumModifier::LABEL_Y;
  else if(event->y
          >= height - graphBorder - graphBottomBorder
          and event->y < height - graphBorder
          and event->x >= graphBorder + graphLeftBorder
          and event->x < width - graphBorder)
    graphModifier = enumModifier::LABEL_X;
  else
  {
    float columnModuleX = (float) (width - graphOffsetStart * 2
                                   - graphLeftBorder)
                          / (residues.size() * 3 + 1);
    float columnModuleY = (float) (height - graphOffsetStart * 2
                                   - graphBottomBorder)
                          / maxPocketLength;
    for(auto i = begin(residues); i != end(residues); i++)
    {
      if(gotcha) break;
      int columnOffsetX = graphOffsetStart + graphLeftBorder
                          + columnModuleX
                            * (3 * distance(begin(residues), i) + 1);
      int columnOffsetY = height - graphOffsetStart
                          - graphBottomBorder;
      for(auto j = begin(i->pockets); j != end(i->pockets); j++)
      {
        int columnHeight = (float)columnModuleY * (*j)->width;
        if(not gotcha and cursorX >= columnOffsetX
           and cursorX < columnOffsetX + columnModuleX * 2
           and cursorY >= columnOffsetY - columnHeight
           and cursorY < columnOffsetY)
        {
          if(hoveringOnPocket != *j)
          {
            newSelection = true;
            hoveringOnPocket = *j;
          }

          graphModifier = enumModifier::POCKET_BAR;
          gotcha = true;
          break;
        }

        columnOffsetY -= columnHeight;
      }
    }

    if(not gotcha)
      graphModifier = enumModifier::NOTHING;
  }

  if(not gotcha)
  {
    if(hoveringOnPocket != nullptr)
    {
      newSelection = true;
      hoveringOnPocket = nullptr;
    }
  }

  if(graphModifier != oldModifier or newSelection)
  {
    updateInformation();
    drawResultsGraph.queue_draw();
  }
  return true;
}

bool
Results::drawResultsGraphButtonPressEvent(GdkEventButton* event) throw ()
{
  if(event->button == 1)
  {
    if(fixedSelection)
    {
      if(selectedPocket == hoveringOnPocket or hoveringOnPocket == nullptr)
        fixedSelection = false;

      selectedPocket = hoveringOnPocket;
    }
    else
    {
      if(hoveringOnPocket != nullptr)
      {
        selectedPocket = hoveringOnPocket;
        fixedSelection = true;
      }
    }

    updateInformation();
    drawResultsGraph.queue_draw();
  }

  return true;
}

void
Results::fillResidues()
{
  const vector<Pocket>& pockets = pittpi->getPockets();
  maxPocketLength = 1;
  unsigned int maxPocketsPerResidue = 0;

  for(auto& pocket : pockets)
  {
    bool exists = false;
    const Residue<SasPdbAtom>& currentRes = pocket.group->getCentralRes();

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


  vector<Color> unscrambledColors;
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

  residues.sort(PocketResidue::sortByResidueIndex);
}

Results::Color
Results::rainbow(double value)
{
  Color color;
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

#if GTKMM_MAJOR == 2
  color.set_rgb_p(red, green, blue);
#else
  color.set_rgba(red, green, blue);
#endif
  return color;
}

void
Results::updateInformation()
{
  const Pocket* pocket;
  if(fixedSelection)
    pocket = selectedPocket;
  else
    pocket = hoveringOnPocket;

  if(pocket)
  {
    stringstream ss;
    const Residue<SasPdbAtom> centralRes(pocket->group->getCentralRes());
    ss << centralRes.index << aminoacidTriplet[centralRes.type];
    entryPocketCenter.set_text(ss.str());
    ss.str("");
    ss << pocket->startPs;
    entryPocketStart.set_text(ss.str());
    ss.str("");
    ss << pocket->endPs;
    entryPocketEnd.set_text(ss.str());
    ss.str("");
    ss << pocket->endPs - pocket->startPs;
    entryPocketWidth.set_text(ss.str());
    ss.str("");
    auto residues(pocket->group->getResidues());
    for(auto residue : residues)
      ss << residue->index << aminoacidTriplet[residue->type] << " ";
    textPocketResidues.get_buffer()->set_text(ss.str());
  }
  else
  {
    entryPocketCenter.set_text("");
    entryPocketStart.set_text("");
    entryPocketEnd.set_text("");
    entryPocketWidth.set_text("");
    textPocketResidues.get_buffer()->set_text("");
  }

  switch(graphModifier)
  {
    case enumModifier::LABEL_X:
    case enumModifier::LABEL_Y:
      drawResultsStatusBar.push(statusBarMessages[1]);
      break;
    case enumModifier::POCKET_BAR:
      if(fixedSelection)
      {
        if(selectedPocket == hoveringOnPocket)
          drawResultsStatusBar.push(statusBarMessages[3]);
        else
          drawResultsStatusBar.push(statusBarMessages[4]);
      }
      else
        drawResultsStatusBar.push(statusBarMessages[2]);
      break;
    case enumModifier::NOTHING:
    default:
      if(fixedSelection)
        drawResultsStatusBar.push(statusBarMessages[5]);
      else
        drawResultsStatusBar.push(statusBarMessages[0]);
      break;
  }
}

void
Results::buildMenu()
{
  static const Glib::ustring menuString =
      "<ui>"
      "  <menubar name='menuBar'>"
      "    <menu action='graphMenu'>"
      "      <menuitem action='changeColors'/>"
      "    </menu>"
      "  </menubar>"
      "</ui>";

  actionGroup = Gtk::ActionGroup::create();
  actionGroup->add(Gtk::Action::create("graphMenu", "_Graph"));
  actionGroup->add(
      Gtk::Action::create("changeColors", "_Change colors",
                     "Change default colors for every bar"),
      sigc::mem_fun(*this, &Results::runColorsChooserDialog));

  uiManager = Gtk::UIManager::create();
  uiManager->insert_action_group(actionGroup);
  uiManager->add_ui_from_string(menuString);
}

void
Results::runColorsChooserDialog()
{
  ColorsChooser colorsChooser;
  colorsChooser.set_colors(colors);
  if(colorsChooser.run() != Gtk::ResponseType::RESPONSE_OK)
    return;

  colors = colorsChooser.get_colors();
  drawResultsGraph.queue_draw();
}

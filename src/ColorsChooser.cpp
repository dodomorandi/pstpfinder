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

#include "ColorsChooser.h"

using namespace PstpFinder;
using namespace std;

ColorsChooser::ColorsChooser()
{
  set_title("Choose your colors");
  set_modal();

#if GTKMM_MAJOR == 2
  Gtk::Box* mainBox(get_vbox());
#else
  Gtk::Box* mainBox(get_content_area());
#endif

  hScalePocketChooser.signal_value_changed().connect(
      sigc::mem_fun(*this,
                    &ColorsChooser::signalHScalePocketChooserValueChanged));
  signal_response().connect(
      sigc::mem_fun(*this, &ColorsChooser::signalResponse));

  mainBox->pack_start(hScalePocketChooser, false, false);
  mainBox->pack_start(colorSelection);

  add_button("Cancel", Gtk::ResponseType::RESPONSE_CANCEL);
  add_button("Ok", Gtk::ResponseType::RESPONSE_OK);
  set_default_response(Gtk::ResponseType::RESPONSE_OK);
}

void
ColorsChooser::set_colors(const vector<Results::Color>& pocketColors)
{
  colors = pocketColors;
  oldColors = colors;
}

vector<Results::Color>
ColorsChooser::get_colors() const
{
  return colors;
}

int
ColorsChooser::run()
{
  oldColorIndex = 0;
  hScalePocketChooser.set_range(1, colors.size());
  hScalePocketChooser.set_increments(1, 1);
  hScalePocketChooser.set_digits(0);
  hScalePocketChooser.set_value(0);
#if GTKMM_MAJOR == 2
  colorSelection.set_current_color(colors[0]);
  colorSelection.set_previous_color(oldColors[0]);
#else
  colorSelection.set_current_rgba(colors[0]);
  colorSelection.set_previous_rgba(oldColors[0]);
#endif

  show_all();
  return Dialog::run();
}

void
ColorsChooser::signalHScalePocketChooserValueChanged()
{
  size_t value(hScalePocketChooser.get_value() - 1);

  if(value == oldColorIndex)
    return;

#if GTKMM_MAJOR == 2
  colors[oldColorIndex] = colorSelection.get_current_color();
  colorSelection.set_current_color(colors[value]);
  colorSelection.set_previous_color(oldColors[value]);
#else
  colors[oldColorIndex] = colorSelection.get_current_rgba();
  colorSelection.set_current_rgba(colors[value]);
  colorSelection.set_previous_rgba(oldColors[value]);
#endif
  oldColorIndex = value;
}

void
ColorsChooser::signalResponse(int responseID)
{
  if(responseID == Gtk::ResponseType::RESPONSE_OK)
#if GTKMM_MAJOR == 2
    colors[oldColorIndex] = colorSelection.get_current_color();
#else
    colors[oldColorIndex] = colorSelection.get_current_rgba();
#endif
}

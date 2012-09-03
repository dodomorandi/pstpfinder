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

#ifndef COLORCHOOSER_H_
#define COLORCHOOSER_H_

#include <vector>

#include <gtkmm.h>
#include <gdkmm.h>

#include "Results.h"

namespace PstpFinder
{
  using namespace Gtk;
  using namespace std;

  class ColorsChooser : public Dialog
  {
    public:
      ColorsChooser();
      void set_colors(const vector<Results::Color>& pocketColors);
      vector<Results::Color> get_colors() const;
      virtual int run();

    private:
      HScale hScalePocketChooser;
      ColorSelection colorSelection;

      vector<Results::Color> colors, oldColors;
      size_t oldColorIndex;

      void signalHScalePocketChooserValueChanged();
      void signalResponse(int responseID);
  };
} /* namespace PstpFinder */
#endif /* COLORCHOOSER_H_ */

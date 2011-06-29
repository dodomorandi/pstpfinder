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

#ifndef _RESULTS_H
#define _RESULTS_H

#include "Pittpi.h"
#include "NewAnalysis.h"
#include <gtkmm.h>

using namespace Gtk;

class NewAnalysis;

namespace Gromacs
{
  struct PocketResidue
  {
    Residue* residue;
    vector<const Pocket*> pockets;

    PocketResidue(const Residue& residue) :
      residue(const_cast<Residue* const>(&residue)) { ; }

    static bool sortByResidueIndex(const PocketResidue& a,
                                   const PocketResidue& b)
    {
      return (a.residue->index < b.residue->index);
    }
  };

  class Results: public Window
  {
    public:
      Results(NewAnalysis& parent, const Pittpi& pittpi,
              const Gromacs& gromacs);
      void init();
    private:
      Notebook notebook;
      DrawingArea drawResultsGraph;
      TextView textViewData, textViewDetails;
      ScrolledWindow scrollData, scrollDetails;

      Gromacs gromacs;
      Pittpi pittpi;
      NewAnalysis& parent;
      vector<PocketResidue> residues;
      float maxPocketLength;
      vector<Gdk::Color> colors;

      bool removeFromParent(GdkEventAny* event);
      bool drawResultsGraphExposeEvent(GdkEventExpose* event);

      void fillResidues();
      static inline Gdk::Color rainbow(double value);
  };
};

#endif

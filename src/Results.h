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

namespace PstpFinder
{
  class NewAnalysis;

  struct PocketResidue
  {
    const Residue<SasPdbAtom>* residue;
    vector<const Pocket*> pockets;

    PocketResidue(const Residue<SasPdbAtom>& residue) :
      residue(&residue) { ; }

    PocketResidue(PocketResidue&& pocketResidue) :
      residue(std::move(pocketResidue.residue))
    {
      pockets = std::move(pocketResidue.pockets);
    }

    static bool sortByResidueIndex(const PocketResidue& a,
                                   const PocketResidue& b)
    {
      return (a.residue->index < b.residue->index);
    }
  };

  class Results: public Window
  {
    public:
#if GTKMM_MAJOR == 2
      typedef Gdk::Color Color;
#else
      typedef Gdk::RGBA Color;
#endif

      Results(NewAnalysis& parent, const shared_ptr<Pittpi>& pittpi,
              const Gromacs& gromacs);
      void init() throw();
    private:
      enum enumModifier
      {
        NOTHING,
        LABEL_X,
        LABEL_Y,
        POCKET_BAR
      };

      // FIXME: When the compiler will accept static initialization lists,
      // FIXME: this must be changed in std::array<std::string, n> and
      // FIXME: initialized here.
      const std::string statusBarMessages[6];

      Notebook notebook;
      DrawingArea drawResultsGraph;
      TextView textViewData, textViewDetails;
      ScrolledWindow scrollData, scrollDetails;
      VBox drawResultsVBox, vboxPocketInformation, vboxMain;
      HPaned panedMain;
      Statusbar drawResultsStatusBar;
      HBox hboxPocketCenter, hboxPocketStart,
           hboxPocketEnd, hboxPocketWidth;
      Label labelPocketCenter, labelPocketStart,
            labelPocketEnd, labelPocketWidth;
      Entry entryPocketCenter, entryPocketStart,
            entryPocketEnd, entryPocketWidth;
      Label labelPocketResidues;
      TextView textPocketResidues;
      Glib::RefPtr<UIManager> uiManager;
      Glib::RefPtr<ActionGroup> actionGroup;

      Gromacs gromacs;
      shared_ptr<Pittpi> pittpi;
      NewAnalysis& parent;
      list<PocketResidue> residues;
      float maxPocketLength;
      vector<Color> colors;
      const Pocket* selectedPocket;
      const Pocket* hoveringOnPocket;
      bool fixedSelection;
      float labelYMultiplier;
      float labelXMultiplier;

      const int graphLineWidth;
      const int graphBorder;
      const int graphOffsetStart;
      int graphLeftBorder;
      int graphBottomBorder;
      int graphHeaderHeight;
      enumModifier graphModifier;

      bool removeFromParent(GdkEventAny* event);
#if GTKMM_MAJOR == 3
      bool drawResultsGraphDrawEvent(
          const Cairo::RefPtr<Cairo::Context>& context) throw ();
#else
      bool drawResultsGraphExposeEvent(GdkEventExpose* event) throw();
#endif
      bool drawResultsGraphScrollEvent(GdkEventScroll* event) throw();
      bool drawResultsGraphMotionEvent(GdkEventMotion* event) throw();
      bool drawResultsGraphButtonPressEvent(GdkEventButton* event) throw();

      void fillResidues();
      static inline Color rainbow(double value);
      void updateInformation();
      void buildMenu();
      void runColorsChooserDialog();
  };
};

#endif

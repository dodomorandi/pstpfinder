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

#ifndef GTKMMWRAPPER_H_
#define GTKMMWRAPPER_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <memory>
#include <gtkmm.h>

namespace PstpFinder
{
  template<typename T>
  class GtkmmWrapper
  {
    public:
      template<typename ... Params>
      GtkmmWrapper(Params&&... parameters):
#if GTKMM_MAJOR == 3
        gtkmmObject(new GtkmmType())
      {
        *gtkmmObject = T::create(parameters...);
      }
#else
        gtkmmObject(new GtkmmType{parameters...})
      {
        ;
      }
#endif

      operator Glib::RefPtr<T>&()
      {
#if GTKMM_MAJOR == 3
        return *gtkmmObject;
#else
        return Glib::RefPtr<T>(&*gtkmmObject);
#endif
      }

      operator T&()
      {
#if GTKMM_MAJOR == 3
        return **gtkmmObject;
#else
        return *gtkmmObject;
#endif
      }

      T* operator ->()
      {
#if GTKMM_MAJOR == 3
        return ((*gtkmmObject).operator ->());
#else
        return &*gtkmmObject;
#endif
      }

    private:
#if GTKMM_MAJOR == 3
      typedef Glib::RefPtr<T> GtkmmType;
#else
      typedef T GtkmmType;
#endif
      std::unique_ptr<GtkmmType> gtkmmObject;
  };
} /* namespace Gromacs */
#endif /* GTKMMWRAPPER_H_ */

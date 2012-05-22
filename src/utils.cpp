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

#include "utils.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#else
#error "Some API implementation is missing for your system."\
       "Please contact program developer"
#endif

#include <cerrno>
#include <string>
// FIXME: Damn!! I want regex!!!!
// #include <regex>

using namespace std;

namespace PstpFinder
{
  bool exists(const string& filename)
  {
#ifdef HAVE_UNISTD_H
    int retval = access(filename.c_str(), R_OK);
    if(retval == -1 and errno != EACCES)
      return false;
    else
      return true;
#endif
  }

  string file_extension(const string& filename)
  {
    /*
    regex extensionRE("\\.[[:w:]]{1,3}$");
    smatch match;
    if(not regex_match(filename, match, extensionRE))
      return "";
    return match.str();
    */
    string lastFour(filename.substr(filename.size() - 4));
    auto dotPosition(lastFour.rfind("."));
    if(dotPosition == string::npos or dotPosition == 3)
      return "";
    else
      return lastFour.substr(dotPosition);
  }

  string change_extension(string filename, const string& new_extension)
  {
    /*
    regex extensionRE("\\.[[:w:]]{1,3}$");
    return regex_replace(filename, extensionRE, new_extension);
    */
    string extension(file_extension(filename));
    return filename.replace(end(filename) - extension.size(), end(filename),
                            new_extension);
  }
}

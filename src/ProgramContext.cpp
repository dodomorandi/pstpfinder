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

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#if GMXVER == 50

#include "ProgramContext.h"
#include <string>

const char*
PstpFinder::ProgramContext::programName() const
{
    return "pstpfinder";
}

const char*
PstpFinder::ProgramContext::displayName() const
{
    return "PSTP-finder";
}

const char*
PstpFinder::ProgramContext::fullBinaryPath() const
{
    /* Makes any differences? */
    return "pstpfinder";
}

const char*
PstpFinder::ProgramContext::defaultLibraryDataPath() const
{
    return GROMACS_PATH "/share/gromacs/top";
}

const char*
PstpFinder::ProgramContext::commandLine() const
{
    /* Makes any differences? */
    return "pstpfinder";
}

#endif /* GMXVER == 50 */

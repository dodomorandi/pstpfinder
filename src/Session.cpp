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

#include "Session.h"
#include <fstream>

using namespace std;

namespace PstpFinder
{
  Session::Session() :
      ready(false),
      sessionFileName()
  {
    ;
  }

  Session::Session(const string & fileName) :
      ready(true),
      sessionFileName(fileName),
      sessionFile(fileName.c_str(), ios::in | ios::binary),
      sasStream(fileName.c_str(), ios::in | ios::binary),
      pdbStream(fileName.c_str(), ios::in | ios::binary)
  {
    readSession(fileName);
  }

  string
  Session::getTrajectoryFileName() const
  {
    if(not ready)
      throw;
    return trajectoryFileName;
  }

  string
  Session::getTopologyFileName() const
  {
    if(not ready)
      throw;
    return topologyFileName;
  }

  unsigned long
  Session::getBeginTime() const
  {
    if(not ready)
      throw;
    return beginTime;
  }

  unsigned long
  Session::getEndTime() const
  {
    if(not ready)
      throw;
    return endTime;
  }

  double
  Session::getRadius() const
  {
    if(not ready)
      throw;
    return radius;
  }

  double
  Session::getPocketThreshold() const
  {
    if(not ready)
      throw;
    return pocketThreshold;
  }

  MetaStream&
  Session::getSasStream()
  {
    if(not ready)
      throw;
    return sasMetaStream;
  }

  unsigned long
  Session::getSasSize() const
  {
    if(not ready)
      throw;
    return sasDataEnd - sasDataStart;
  }

  MetaStream&
  Session::getPdbStream()
  {
    if(not ready)
      throw;
    return pdbMetaStream;
  }

  unsigned long
  Session::getPdbSize() const
  {
    if(not ready)
      throw;
    return pdbDataEnd - pdbDataStart;
  }

  void
  Session::readSession(const string & fileName)
  {
    std::locale oldLocale;
    std::locale::global(std::locale("C"));

    unsigned int dataUInt;
    sessionFile.seekg(0);
    getline(sessionFile, trajectoryFileName);
    getline(sessionFile, topologyFileName);
    sessionFile >> beginTime;
    sessionFile >> endTime;
    sessionFile >> radius;
    sessionFile >> pocketThreshold;

    sessionFile >> dataUInt;
    if(sessionFile.peek() == '\n')
      (void) (sessionFile.get());
    sasDataStart = sessionFile.tellg();
    sessionFile.seekg(dataUInt, ios::cur);
    sasDataEnd = sessionFile.tellg();
    sasMetaStream = MetaStream(sasStream, sasDataStart, sasDataEnd);

    sessionFile >> dataUInt;
    if(sessionFile.peek() == '\n')
      (void) (sessionFile.get());
    pdbDataStart = sessionFile.tellg();
    sessionFile.seekg(dataUInt, ios::cur);
    pdbDataEnd = sessionFile.tellg();
    pdbMetaStream = MetaStream(pdbStream, pdbDataStart, pdbDataEnd);

    std::locale::global(oldLocale);
  }
/* namespace PstpFinder */
}


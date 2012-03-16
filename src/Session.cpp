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
#include <boost/filesystem.hpp>

using namespace std;

namespace PstpFinder
{
  Session::Session() :
      ready(false), rawSasSession(false), sessionFileName()
  {
    ;
  }

  Session::Session(const string & fileName) :
      ready(true),
      rawSasSession(false),
      sessionFileName(fileName),
      sessionFile(fileName.c_str(), ios::in | ios::binary),
      sasStream(fileName.c_str(), ios::in | ios::binary),
      pdbStream(fileName.c_str(), ios::in | ios::binary)
  {
    readSession(fileName);
  }

  Session::~Session()
  {
    if(ready)
    {
      delete sasMetaStream;
      sasStream.close();
      if(not rawSasSession)
      {
        delete pdbMetaStream;
        pdbStream.close();
      }
      sessionFile.close();
    }
  }

  string
  Session::getTrajectoryFileName() const
  {
    if(not ready or rawSasSession)
      throw;
    return trajectoryFileName;
  }

  string
  Session::getTopologyFileName() const
  {
    if(not ready or rawSasSession)
      throw;
    return topologyFileName;
  }

  unsigned long
  Session::getBeginTime() const
  {
    if(not ready or rawSasSession)
      throw;
    return beginTime;
  }

  unsigned long
  Session::getEndTime() const
  {
    if(not ready or rawSasSession)
      throw;
    return endTime;
  }

  double
  Session::getRadius() const
  {
    if(not ready or rawSasSession)
      throw;
    return radius;
  }

  double
  Session::getPocketThreshold() const
  {
    if(not ready or rawSasSession)
      throw;
    return pocketThreshold;
  }

  MetaStream<ifstream>&
  Session::getSasStream()
  {
    if(not ready)
      throw;
    return *sasMetaStream;
  }

  unsigned long
  Session::getSasSize() const
  {
    if(not ready)
      throw;
    return sasDataEnd - sasDataStart;
  }

  const bool
  Session::isRawSasSession() const
  {
    return(rawSasSession);
  }

  MetaStream<ifstream>&
  Session::getPdbStream()
  {
    if(not ready or rawSasSession)
      throw;
    return *pdbMetaStream;
  }

  unsigned long
  Session::getPdbSize() const
  {
    if(not ready or rawSasSession)
      throw;
    return pdbDataEnd - pdbDataStart;
  }

  Session&
  Session::operator =(const Session& session)
  {
    if(&session == this)
      return *this;

    this->~Session();
    if(session.sessionFileName == "")
      new (this) Session();
    else
    {
      new (this) Session(session.sessionFileName);
      sasMetaStream->seekg(session.sasMetaStream->tellg());
      if(not session.rawSasSession)
        pdbMetaStream->seekg(session.pdbMetaStream->tellg());
    }

    return *this;
  }

  void
  Session::readSession(const string & fileName)
  {
    std::locale oldLocale;
    std::locale::global(std::locale("C"));

    unsigned int dataUInt;
    sessionFile.seekg(0);
    getline(sessionFile, trajectoryFileName);
    if(boost::filesystem::exists(boost::filesystem::path(trajectoryFileName)))
    {
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
      sasMetaStream = new MetaStream<ifstream>(sasStream, sasDataStart, sasDataEnd);

      sessionFile >> dataUInt;
      if(sessionFile.peek() == '\n')
        (void) (sessionFile.get());
      pdbDataStart = sessionFile.tellg();
      sessionFile.seekg(dataUInt, ios::cur);
      pdbDataEnd = sessionFile.tellg();
      pdbMetaStream = new MetaStream<ifstream>(pdbStream, pdbDataStart, pdbDataEnd);
    }
    else
    {
      rawSasSession = true;
      sessionFile.seekg(0);
      sasDataStart = sessionFile.tellg();
      sessionFile.seekg(0, ios_base::end);
      sasDataEnd = sessionFile.tellg();
      sasMetaStream = new MetaStream<ifstream>(sasStream, sasDataStart, sasDataEnd);
    }

    std::locale::global(oldLocale);
  }
/* namespace PstpFinder */
}


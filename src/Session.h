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

#ifndef SESSION_H_
#define SESSION_H_

#include "MetaStream.h"
#include <string>
#include <iostream>
#include <fstream>

using namespace std;

namespace PstpFinder
{
  class Session
  {
    public:
      Session();
      Session(const string& fileName);
      ~Session();
      string getTrajectoryFileName() const;
      string getTopologyFileName() const;
      unsigned long getBeginTime() const;
      unsigned long getEndTime() const;
      double getRadius() const;
      double getPocketThreshold() const;
      MetaStream<ifstream>& getSasStream();
      unsigned long getSasSize() const;
      const bool isRawSasSession() const;
      MetaStream<ifstream>& getPdbStream();
      unsigned long getPdbSize() const;
      const bool isPittpiAvailable() const;
      unsigned long getPittpiSize() const;
      MetaStream<ifstream>& getPittpiStream();

      Session& operator =(const Session& session);

    private:
      const bool ready;
      bool rawSasSession;
      bool pittpiAvailable;
      const string sessionFileName;
      ifstream sessionFile;
      ifstream sasStream;
      ifstream pdbStream;
      ifstream pittpiStream;
      string trajectoryFileName;
      string topologyFileName;
      unsigned long beginTime;
      unsigned long endTime;
      double radius;
      double pocketThreshold;
      streampos sasDataStart;
      streampos sasDataEnd;
      streampos pdbDataStart;
      streampos pdbDataEnd;
      streampos pittpiDataStart;
      streampos pittpiDataEnd;
      MetaStream<ifstream>* sasMetaStream;
      MetaStream<ifstream>* pdbMetaStream;
      MetaStream<ifstream>* pittpiMetaStream;

      void
      readSession(const string& fileName);
  };

} /* namespace PstpFinder */
#endif /* SESSION_H_ */

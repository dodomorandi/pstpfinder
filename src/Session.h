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
#include <memory>
#include <type_traits>

#include <boost/filesystem.hpp>

using namespace std;

namespace PstpFinder
{
  template<typename T>
  class Session
  {
    public:
      Session();
      Session(const string& fileName);
      Session& operator =(const Session&) = delete;
      string getTrajectoryFileName() const;
      string getTopologyFileName() const;
      unsigned long getBeginTime() const;
      unsigned long getEndTime() const;
      double getRadius() const;
      double getPocketThreshold() const;
      MetaStream<T>& getSasStream();
      unsigned long getSasSize() const;
      const bool isRawSasSession() const;
      MetaStream<T>& getPdbStream();
      unsigned long getPdbSize() const;
      const bool isPittpiAvailable() const;
      unsigned long getPittpiSize() const;
      MetaStream<T>& getPittpiStream();

    private:
      const bool ready;
      bool rawSasSession;
      bool pittpiAvailable;
      const string sessionFileName;
      T sessionFile;
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
      unique_ptr<MetaStream<T>> sasMetaStream;
      unique_ptr<MetaStream<T>> pdbMetaStream;
      unique_ptr<MetaStream<T>> pittpiMetaStream;

      void readSession(const string& fileName);
      inline void assertRegularType() const;
  };

  template<typename T>
  Session<T>::Session() :
      ready(false), rawSasSession(false), sessionFileName()
  {
    assertRegularType();
  }

  template<typename T>
  Session<T>::Session(const string & fileName) :
      ready(true),
      rawSasSession(false),
      sessionFileName(fileName),
      sessionFile(fileName.c_str(), ios::in | ios::binary)
  {
    assertRegularType();
    readSession(fileName);
  }

  template<typename T>
  string
  Session<T>::getTrajectoryFileName() const
  {
    assert(ready);
    assert(not rawSasSession);
    return trajectoryFileName;
  }

  template<typename T>
  string
  Session<T>::getTopologyFileName() const
  {
    assert(ready);
    assert(not rawSasSession);
    return topologyFileName;
  }

  template<typename T>
  unsigned long
  Session<T>::getBeginTime() const
  {
    assert(ready);
    assert(not rawSasSession);
    return beginTime;
  }

  template<typename T>
  unsigned long
  Session<T>::getEndTime() const
  {
    assert(ready);
    assert(not rawSasSession);
    return endTime;
  }

  template<typename T>
  double
  Session<T>::getRadius() const
  {
    assert(ready);
    assert(not rawSasSession);
    return radius;
  }

  template<typename T>
  double
  Session<T>::getPocketThreshold() const
  {
    assert(ready);
    assert(not rawSasSession);
    return pocketThreshold;
  }

  template<typename T>
  MetaStream<T>&
  Session<T>::getSasStream()
  {
    assert(ready);
    return *sasMetaStream;
  }

  template<typename T>
  unsigned long
  Session<T>::getSasSize() const
  {
    assert(ready);
    return sasDataEnd - sasDataStart;
  }

  template<typename T>
  const bool
  Session<T>::isRawSasSession() const
  {
    return(rawSasSession);
  }

  template<typename T>
  MetaStream<T>&
  Session<T>::getPdbStream()
  {
    assert(ready);
    assert(not rawSasSession);
    return *pdbMetaStream;
  }

  template<typename T>
  unsigned long
  Session<T>::getPdbSize() const
  {
    assert(ready);
    assert(not rawSasSession);
    return pdbDataEnd - pdbDataStart;
  }

  template<typename T>
  const bool
  Session<T>::isPittpiAvailable() const
  {
    return pittpiAvailable;
  }

  template<typename T>
  MetaStream<T>&
  Session<T>::getPittpiStream()
  {
    assert(ready);
    assert(not rawSasSession);
    assert(pittpiAvailable);
    return *pittpiMetaStream;
  }

  template<typename T>
  unsigned long
  Session<T>::getPittpiSize() const
  {
    assert(ready);
    assert(not rawSasSession);
    assert(pittpiAvailable);
    return pittpiDataEnd - pittpiDataStart;
  }

  template<typename T>
  void
  Session<T>::readSession(const string & fileName)
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
      sasMetaStream = unique_ptr<MetaStream<T>>(
            new MetaStream<T>(fileName,
                              ios_base::in | ios_base::out | ios_base::binary,
                              enumStreamType::STREAMTYPE_FIXED, sasDataStart,
                              sasDataEnd));
      if(sessionFile.peek() == '\n')
        (void) (sessionFile.get());

      sessionFile >> dataUInt;
      if(sessionFile.peek() == '\n')
        (void) (sessionFile.get());
      pdbDataStart = sessionFile.tellg();
      sessionFile.seekg(dataUInt, ios::cur);
      pdbDataEnd = sessionFile.tellg();
      pdbMetaStream = unique_ptr<MetaStream<T>>(
            new MetaStream<T>(fileName,
                              ios_base::in | ios_base::out | ios_base::binary,
                              enumStreamType::STREAMTYPE_FIXED, pdbDataStart,
                              pdbDataEnd));
      if(sessionFile.peek() == '\n')
        (void) (sessionFile.get());

      sessionFile.peek();
      if(not sessionFile.eof())
      {
        pittpiAvailable = true;
        sessionFile >> dataUInt;
        if(sessionFile.peek() == '\n')
          (void) (sessionFile.get());
        pittpiDataStart = sessionFile.tellg();
        sessionFile.seekg(dataUInt, ios::cur);
        pittpiDataEnd = sessionFile.tellg();
        pittpiMetaStream = unique_ptr<MetaStream<T>>(
              new MetaStream<T>(fileName,
                                ios_base::in | ios_base::out | ios_base::binary,
                                enumStreamType::STREAMTYPE_FIXED,
                                pittpiDataStart, pittpiDataEnd));
      }
      else
        pittpiAvailable = false;
    }
    else
    {
      rawSasSession = true;
      pittpiAvailable = false;
      sessionFile.seekg(0);
      sasDataStart = sessionFile.tellg();
      sessionFile.seekg(0, ios_base::end);
      sasDataEnd = sessionFile.tellg();
      sasMetaStream = unique_ptr<MetaStream<T>>(
            new MetaStream<T>(fileName,
                              ios_base::in | ios_base::out | ios_base::binary,
                              enumStreamType::STREAMTYPE_FIXED, sasDataStart,
                              sasDataEnd));
    }

    std::locale::global(oldLocale);
  }

  template<typename T>
  inline void
  Session<T>::assertRegularType() const
  {
    static_assert(is_base_of<base_stream(basic_istream), T>::value or
                  is_base_of<base_stream(basic_ostream), T>::value,
                  "T must have basic_istream or basic_ostream as base class");
  }
} /* namespace PstpFinder */
#endif /* SESSION_H_ */

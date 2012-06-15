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

#define SESSION_VERSION 1
namespace PstpFinder
{
  // Session forward declarations for Gromacs.h (and maybe others)
  template<typename T, typename = void> class Session;
}

#include "MetaStream.h"
#include "Gromacs.h"
#include "utils.h"
#include "Serializer.h"

#include <string>
#include <iostream>
#include <sstream>
#include <memory>
#include <type_traits>
#include <bitset>
#include <initializer_list>
#include <tuple>
#include <utility>

using namespace std;

namespace PstpFinder
{
  enum class SessionParameter
  {
    TRAJECTORY,
    TOPOLOGY,
    BEGIN,
    END,
    RADIUS,
    THRESHOLD
  };

  // FIXME: I'd like to use a union, but std::string has non trivial
  // constructor, copy constructor and destructor.
  struct SessionParameterValue
  {
    string str;
    unsigned long ulong;
    double dbl;
  };

  inline tuple<SessionParameter, SessionParameterValue>
  make_sessionParameter(SessionParameter parameter, const string& value)
  {
    SessionParameterValue spv;
    spv.str = value;
    return make_tuple(parameter, spv);
  }

  inline tuple<SessionParameter, SessionParameterValue>
  make_sessionParameter(SessionParameter parameter, unsigned long value)
  {
    SessionParameterValue spv;
    spv.ulong = value;
    return make_tuple(parameter, spv);
  }

  inline tuple<SessionParameter, SessionParameterValue>
  make_sessionParameter(SessionParameter parameter, double value)
  {
    SessionParameterValue spv;
    spv.dbl = value;
    return make_tuple(parameter, spv);
  }

  template<typename T>
  class Session_Base
  {
    public:
      typedef MetaStream<T> stream_type;
      Session_Base& operator =(const Session_Base&) = delete;
      string getTrajectoryFileName() const;
      string getTopologyFileName() const;
      unsigned long getBeginTime() const;
      unsigned long getEndTime() const;
      double getRadius() const;
      double getPocketThreshold() const;
      stream_type& getSasStream();
      unsigned long getSasSize() const;
      bool sasComplete() const;
      stream_type& getPdbStream();
      unsigned long getPdbSize() const;
      bool pdbComplete() const;
      const bool isPittpiAvailable() const;
      unsigned long getPittpiSize() const;
      stream_type& getPittpiStream();
      bool pittpiComplete() const;
      void abort();

      void eventSasStreamClosing();
      void eventPdbStreamClosing();
      void eventPittpiStreamClosing();

    protected:
      friend class MetaStream<T>;
      string trajectoryFileName;
      string topologyFileName;
      unsigned long beginTime;
      unsigned long endTime;
      double radius;
      double pocketThreshold;
      bitset<6> parameterSet;

      Session_Base();
      Session_Base(const string& fileName);
      Session_Base(const string& fileName,
          initializer_list<tuple<SessionParameter,
                                 SessionParameterValue>>&& parameters);
      void readSession();
      void prepareForWrite();
      inline void assertRegularType() const;
      inline void assertBaseIStream() const;
      inline void assertBaseOStream() const;
      inline void assertBaseIOStream() const;

    private:
      struct MetaData
      {
        unique_ptr<stream_type> stream;
        streampos info;
        streampos start;
        streampos end;
        bool complete;

        MetaData() : info(-1), start(-1), end(-1), complete(false) {};
      };

      const bool ready;
      unsigned short version;
      const string sessionFileName;
      unique_ptr<T> sessionFile;
      unique_ptr<Serializer<T>> serializer;
      MetaData metaSas;
      MetaData metaPdb;
      MetaData metaPittpi;
  };

  template<typename T>
  Session_Base<T>::Session_Base() :
      ready(false), version(0), sessionFileName()
  {
    assertRegularType();
  }

  template<typename T>
  Session_Base<T>::Session_Base(const string& fileName) :
      ready(true),
      version(0),
      sessionFileName(fileName)
  {
    assertRegularType();

    if(is_base_of<base_stream(basic_istream, T), T>::value)
      sessionFile = unique_ptr<T>(new T(fileName,
                                        ios_base::in bitor ios_base::out bitor
                                        ios_base::binary));
    else
      sessionFile = unique_ptr<T>(new T(fileName,
                                        ios_base::out bitor ios_base::binary));
    serializer = unique_ptr<Serializer<T>>(new Serializer<T>(*sessionFile));
  }

  template<typename T>
  Session_Base<T>::Session_Base(
        const string& fileName,
        initializer_list<tuple<SessionParameter,
                               SessionParameterValue>>&& parameters) :
      ready(true),
      version(0),
      sessionFileName(fileName)
  {
    assertRegularType();

    if(is_base_of<base_stream(basic_istream, T), T>::value)
      sessionFile = unique_ptr<T>(new T(fileName,
                                        ios_base::in bitor ios_base::out bitor
                                        ios_base::binary));
    else
      sessionFile = unique_ptr<T>(new T(fileName,
                                        ios_base::out bitor ios_base::binary));
    serializer = unique_ptr<Serializer<T>>(new Serializer<T>(*sessionFile));

    for(auto& parameter : parameters)
    {
      switch(get<0>(parameter))
      {
        case SessionParameter::TRAJECTORY:
          trajectoryFileName = get<1>(parameter).str;
          parameterSet |= 1;
          break;
        case SessionParameter::TOPOLOGY:
          topologyFileName = get<1>(parameter).str;
          parameterSet |= 2;
          break;
        case SessionParameter::BEGIN:
          beginTime = get<1>(parameter).ulong;
          parameterSet |= 4;
          break;
        case SessionParameter::END:
          endTime = get<1>(parameter).ulong;
          parameterSet |= 8;
          break;
        case SessionParameter::RADIUS:
          radius = get<1>(parameter).dbl;
          parameterSet |= 16;
          break;
        case SessionParameter::THRESHOLD:
          pocketThreshold = get<1>(parameter).dbl;
          parameterSet |= 32;
          break;
      }
    }
  }

  template<typename T>
  string
  Session_Base<T>::getTrajectoryFileName() const
  {
    assert(ready);
    return trajectoryFileName;
  }

  template<typename T>
  string
  Session_Base<T>::getTopologyFileName() const
  {
    assert(ready);
    return topologyFileName;
  }

  template<typename T>
  unsigned long
  Session_Base<T>::getBeginTime() const
  {
    assert(ready);
    return beginTime;
  }

  template<typename T>
  unsigned long
  Session_Base<T>::getEndTime() const
  {
    assert(ready);
    return endTime;
  }

  template<typename T>
  double
  Session_Base<T>::getRadius() const
  {
    assert(ready);
    return radius;
  }

  template<typename T>
  double
  Session_Base<T>::getPocketThreshold() const
  {
    assert(ready);
    return pocketThreshold;
  }

  template<typename T>
  typename Session_Base<T>::stream_type&
  Session_Base<T>::getSasStream()
  {
    assert(ready);
    assert(metaSas.stream);
    return *metaSas.stream;
  }

  template<typename T>
  unsigned long
  Session_Base<T>::getSasSize() const
  {
    assert(ready);
    if(metaSas.stream)
      return metaSas.end - metaSas.start;
    else
      return 0;
  }

  template<typename T>
  bool
  Session_Base<T>::sasComplete() const
  {
    assert(ready);
    return metaSas.stream and metaSas.complete;
  }

  template<typename T>
  typename Session_Base<T>::stream_type&
  Session_Base<T>::getPdbStream()
  {
    assert(ready);
    assert(metaPdb.stream);
    return *metaPdb.stream;
  }

  template<typename T>
  unsigned long
  Session_Base<T>::getPdbSize() const
  {
    assert(ready);
    if(metaPdb.stream)
      return metaPdb.end - metaPdb.start;
    else
      return 0;
  }

  template<typename T>
  bool
  Session_Base<T>::pdbComplete() const
  {
    assert(ready);
    return metaPdb.stream and metaPdb.complete;
  }

  template<typename T>
  const bool
  Session_Base<T>::isPittpiAvailable() const
  {
    return version > 1;
  }

  template<typename T>
  typename Session_Base<T>::stream_type&
  Session_Base<T>::getPittpiStream()
  {
    assert(ready);
    assert(version > 1);
    assert(metaPittpi.stream);
    return *metaPittpi.stream;
  }

  template<typename T>
  unsigned long
  Session_Base<T>::getPittpiSize() const
  {
    assert(ready);
    assert(version > 1);
    if(metaPittpi.stream)
      return metaPittpi.end - metaPittpi.start;
    else
      return 0;
  }

  template<typename T>
  bool
  Session_Base<T>::pittpiComplete() const
  {
    assert(ready);
    return metaPittpi.stream and metaPittpi.complete;
  }

  template<typename T>
  void
  Session_Base<T>::readSession()
  {
    unsigned int dataUInt;
    sessionFile->seekg(0);
    sessionFile->peek();
    if(sessionFile->eof())
      return;

    *serializer >> version;
    *serializer >> trajectoryFileName;
    *serializer >> topologyFileName;
    *serializer >> beginTime;
    *serializer >> endTime;
    *serializer >> radius;
    *serializer >> pocketThreshold;

    parameterSet.set();

    metaSas.info = sessionFile->tellg();
    *serializer >> dataUInt;
    metaSas.start = sessionFile->tellg();

    if(dataUInt > 0)
    {
      metaSas.complete = true;
      sessionFile->seekg(dataUInt, ios_base::cur);
      metaSas.end = sessionFile->tellg();
      metaSas.stream = unique_ptr<stream_type>(
        new stream_type(sessionFileName,
                        ios_base::in | ios_base::out | ios_base::binary,
                        enumStreamType::STREAMTYPE_FIXED, metaSas.start,
                        metaSas.end));
    }
    else
    {
      sessionFile->seekg(0, ios_base::end);
      metaSas.end = 0;
      metaSas.stream = unique_ptr<stream_type>(
        new stream_type(sessionFileName,
                        ios_base::in | ios_base::out | ios_base::binary,
                        enumStreamType::STREAMTYPE_ADJUST, metaSas.start));
      return;
    }

    if(sessionFile->eof())
      return;

    metaPdb.info = sessionFile->tellg();

    sessionFile->peek();
    if(sessionFile->eof())
      return;

    *serializer >> dataUInt;
    metaPdb.start = sessionFile->tellg();
    if(dataUInt > 0)
    {
      metaPdb.complete = true;
      sessionFile->seekg(dataUInt, ios_base::cur);
      metaPdb.end = sessionFile->tellg();
      metaPdb.stream = unique_ptr<stream_type>(
        new stream_type(sessionFileName,
                        ios_base::in | ios_base::out | ios_base::binary,
                        enumStreamType::STREAMTYPE_FIXED, metaPdb.start,
                        metaPdb.end));
    }
    else
    {
      sessionFile->seekg(0, ios_base::end);
      metaPdb.end = 0;
      metaPdb.stream = unique_ptr<stream_type>(
        new stream_type(sessionFileName,
                        ios_base::in | ios_base::out | ios_base::binary,
                        enumStreamType::STREAMTYPE_ADJUST, metaPdb.start));
      return;
    }

    if(version > 1)
    {
      if(sessionFile->eof())
        return;

      metaPittpi.info = sessionFile->tellg();

      sessionFile->peek();
      if(sessionFile->eof())
        return;

      *serializer >> dataUInt;;
      metaPittpi.start = sessionFile->tellg();
      if(dataUInt > 0)
      {
        metaPittpi.complete = true;
        sessionFile->seekg(dataUInt, ios_base::cur);
        metaPittpi.end = sessionFile->tellg();
        metaPittpi.stream = unique_ptr<stream_type>(
          new stream_type(sessionFileName,
                          ios_base::in | ios_base::out | ios_base::binary,
                          enumStreamType::STREAMTYPE_FIXED, metaPittpi.start,
                          metaPittpi.end));
      }
      else
      {
        sessionFile->seekg(0, ios_base::end);
        metaPittpi.end = 0;
        metaPittpi.stream = unique_ptr<stream_type>(
          new stream_type(sessionFileName,
                          ios_base::in | ios_base::out | ios_base::binary,
                          enumStreamType::STREAMTYPE_ADJUST,
                          metaPittpi.start));
        return;
      }
    }
  }

  template<typename T>
  void
  Session_Base<T>::prepareForWrite()
  {
    assert(parameterSet.all()); // Has radius and pocketThreshold set

    unsigned int dataUInt;
    switch(version)
    {
      case 0:   // Session have not been read or session is empty
        version = SESSION_VERSION;
        *serializer << version;
        *serializer << trajectoryFileName;
        *serializer << topologyFileName;
        *serializer << beginTime;
        *serializer << endTime;
        *serializer << radius;
        *serializer << pocketThreshold;

        metaSas.info = sessionFile->tellp();
        metaPdb.info = -1;
        metaPittpi.info = -1;
        dataUInt = 0;
        *serializer << dataUInt;
        metaSas.start = sessionFile->tellp();
        metaSas.complete = false;
        metaSas.stream = unique_ptr<stream_type>(
          new stream_type(sessionFileName,
                          ios_base::in | ios_base::out | ios_base::binary,
                          enumStreamType::STREAMTYPE_ADJUST, metaSas.start));

        metaSas.stream->callbackClose = bind(
              &Session_Base<T>::eventSasStreamClosing, ref(*this));
        break;
      case 1:  // SAS + PDB
        if(metaSas.end == 0)
        {
          metaSas.complete = false;
          metaSas.stream = unique_ptr<stream_type>(
            new stream_type(sessionFileName,
                            ios_base::in | ios_base::out | ios_base::binary,
                            enumStreamType::STREAMTYPE_ADJUST,
                            metaSas.start));

          metaSas.stream->callbackClose = bind(
              &Session_Base<T>::eventSasStreamClosing, ref(*this));
        }
        else if(metaPdb.end == 0)
        {
          metaSas.complete = true;
          metaPdb.complete = false;
          metaPdb.stream = unique_ptr<stream_type>(
            new stream_type(sessionFileName,
                            ios_base::in | ios_base::out | ios_base::binary,
                            enumStreamType::STREAMTYPE_ADJUST,
                            metaPdb.start));

          metaPdb.stream->callbackClose = bind(
              &Session_Base<T>::eventPdbStreamClosing, ref(*this));
        }
        break;
    }

    sessionFile->flush();
  }

  template<typename T>
  void
  Session_Base<T>::abort()
  {
    if(metaSas.stream and metaSas.stream->callbackClose)
      metaSas.stream->callbackClose = function<void()>();
    if(metaPdb.stream and metaPdb.stream->callbackClose)
      metaPdb.stream->callbackClose = function<void()>();
    if(metaPittpi.stream and metaPittpi.stream->callbackClose)
      metaPittpi.stream->callbackClose = function<void()>();
  }

  template<typename T>
  inline void
  Session_Base<T>::assertRegularType() const
  {
    static_assert(is_base_of<base_stream(basic_istream, T), T>::value or
                  is_base_of<base_stream(basic_ostream, T), T>::value,
                  "T must have basic_istream or basic_ostream as base class");
  }

  template<typename T>
  inline void
  Session_Base<T>::assertBaseIStream() const
  {
    static_assert(is_base_of<base_stream(basic_istream, T), T>::value,
                  "T must have basic_istream as base class");
  }

  template<typename T>
  inline void
  Session_Base<T>::assertBaseOStream() const
  {
    static_assert(is_base_of<base_stream(basic_ostream, T), T>::value,
                  "T must have basic_ostream as base class");
  }

  template<typename T>
  inline void
  Session_Base<T>::assertBaseIOStream() const
  {
    static_assert(is_base_of<base_stream(basic_istream, T), T>::value and
                  is_base_of<base_stream(basic_ostream, T), T>::value,
                  "T must have basic_istream and basic_ostream as base class");
  }

  template<typename T>
  void
  Session_Base<T>::eventSasStreamClosing()
  {
    if(not metaSas.stream or not metaSas.stream->is_open())
      return;

    metaSas.stream->seekp(0, ios_base::end);
    streampos endOfSas(metaSas.stream->tellp());
    metaSas.end = metaSas.start + endOfSas;
    metaSas.complete = true;
    sessionFile->seekp(metaSas.info);
    unsigned int dataUInt(endOfSas);
    *serializer << dataUInt;

    sessionFile->seekp(metaSas.end);
    metaPdb.info = sessionFile->tellp();
    dataUInt = 0;
    *serializer << dataUInt;
    metaPdb.complete = false;
    metaPdb.start = sessionFile->tellp();
    metaPdb.stream = unique_ptr<stream_type>(
        new stream_type(sessionFileName,
                        ios_base::in | ios_base::out | ios_base::binary,
                        enumStreamType::STREAMTYPE_ADJUST, metaPdb.start));

    metaPdb.stream->callbackClose = bind(
          &Session_Base<T>::eventPdbStreamClosing, ref(*this));
  }

  template<typename T>
  void
  Session_Base<T>::eventPdbStreamClosing()
  {
    if(not metaPdb.stream or not metaPdb.stream->is_open())
      return;

    metaPdb.stream->seekp(0, ios_base::end);
    streampos endOfPdb(metaPdb.stream->tellp());
    metaPdb.end = metaPdb.start + endOfPdb;
    metaPdb.complete = true;
    sessionFile->seekp(metaPdb.info);
    unsigned int dataUInt(endOfPdb);
    *serializer << dataUInt;

    if(version > 1)
    {
      sessionFile->seekp(metaPdb.end);
      metaPittpi.info = sessionFile->tellp();
      dataUInt = 0;
      *serializer << dataUInt;
      metaPittpi.complete = false;
      metaPittpi.start = sessionFile->tellp();
      metaPittpi.stream = unique_ptr<stream_type>(
          new stream_type(sessionFileName,
                          ios_base::in | ios_base::out | ios_base::binary,
                          enumStreamType::STREAMTYPE_ADJUST,
                          metaPittpi.start));

      metaPittpi.stream->callbackClose = bind(
            &Session_Base<T>::eventPittpiStreamClosing, ref(*this));
    }
    else
      sessionFile->close();
  }

  template<typename T>
  void
  Session_Base<T>::eventPittpiStreamClosing()
  {
    // TODO: closing a Pittpi Stream
  }

  template<typename T>
  class Session<T, typename enable_if<
          is_base_of<base_stream(basic_istream, T), T>::value and
          not is_base_of<base_stream(basic_ostream, T), T>::value>::type>
      : public Session_Base<T>
  {
    public:
      Session() : Base() {}
      Session(const string& fileName) : Base(fileName)
      {
        Base::assertBaseIStream();
        Base::readSession();
      }
      Session& operator =(const Session&) = delete;
    private:
      typedef Session_Base<T> Base;
  };

  template<typename T>
  class Session<T, typename enable_if<
          not is_base_of<base_stream(basic_istream, T), T>::value and
          is_base_of<base_stream(basic_ostream, T), T>::value>::type>
      : public Session_Base<T>
  {
    public:
      Session() : Base() {}
      Session(const string& fileName, Gromacs& gromacs, double radius,
              double pocketThreshold) :
          Base(fileName,
                { make_sessionParameter(SessionParameter::TRAJECTORY,
                                        gromacs.getTrajectoryFile()),
                  make_sessionParameter(SessionParameter::TOPOLOGY,
                                        gromacs.getTopologyFile()),
                  make_sessionParameter(
                      SessionParameter::BEGIN,
                      static_cast<unsigned long>(gromacs.getBegin())),
                  make_sessionParameter(
                      SessionParameter::END,
                      static_cast<unsigned long>(gromacs.getEnd())),
                  make_sessionParameter(SessionParameter::RADIUS, radius),
                  make_sessionParameter(SessionParameter::THRESHOLD,
                                        pocketThreshold) })
      {
        Base::assertBaseOStream();
        Base::prepareForWrite();
      }
      Session& operator =(const Session&) = delete;
    private:
      typedef Session_Base<T> Base;
  };

  template<typename T>
  class Session<T, typename enable_if<
          is_base_of<base_stream(basic_istream, T), T>::value and
          is_base_of<base_stream(basic_ostream, T), T>::value>::type>
      : public Session_Base<T>
  {
    public:
      Session() : Base() {}
      Session(const string& fileName) : Base(fileName)
      {
        Base::assertBaseIOStream();
        Base::readSession();
        Base::prepareForWrite();
      }
      Session& operator =(const Session&) = delete;

    private:
      typedef Session_Base<T> Base;
  };
} /* namespace PstpFinder */
#endif /* SESSION_H_ */

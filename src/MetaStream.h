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

#ifndef METASTREAM_H_
#define METASTREAM_H_

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <iosfwd>
#include <boost/iostreams/categories.hpp>

using namespace std;
namespace io = boost::iostreams;

namespace PstpFinder
{
  template<typename T>
  class MetaStream
  {
    public:
      typedef char char_type;
      typedef io::source_tag category;

      MetaStream() :
            stream(nullStream), valid(false), eofTrigger(false) {}

      MetaStream(T& modifiableStream, streampos begin = -1,
                   streampos end = -1) :
            stream(modifiableStream), valid(true), eofTrigger(false)
      {
        assert_basic_istream();

        if(begin == -1)
          streamBegin = stream.tellg();
        else
          streamBegin = begin;

        if(end == -1)
        {
          stream.seekg(0, ios_base::end);
          streamEnd = stream.tellg();
        }
        else
          streamEnd = end;

        stream.seekg(streamBegin, ios_base::beg);
      }

      MetaStream(const MetaStream& metaStream) :
            stream(metaStream.stream),
            valid(metaStream.valid),
            eofTrigger(false)
      {
        streamBegin = metaStream.streamBegin;
        streamEnd = metaStream.streamEnd;
      }

      // basic_istream related functions. Needed a static_assert for these.
      template<typename U>
      MetaStream&
      operator >>(U& out)
      {
        assert_basic_istream();

        if(not valid or eofTrigger)
          throw;

        if(eof())
          eofTrigger = true;
        streamoff currentPosition = stream.tellg();
        streampos finalPosition = currentPosition + sizeof(T);

        if(finalPosition > streamEnd)
        {
          stringstream tempStream(stringstream::in | stringstream::out);
          streamsize tempStreamSize = streamEnd - currentPosition;
          char *remainingData = new char[tempStreamSize];
          stream.read(remainingData, tempStreamSize);
          tempStream.write(remainingData, tempStreamSize);
          tempStream >> out;
          delete[] remainingData;
        }
        else
          stream >> out;
      }

      streamsize
      read(char* data, streamsize length)
      {
        assert_basic_istream();

        streamsize size;
        if(not valid or eofTrigger)
          throw;

        streamoff currentPosition = stream.tellg();

        if(currentPosition + length > streamEnd)
        {
          stream.read(data, streamEnd - currentPosition);
          size = streamEnd - currentPosition;
        }
        else
        {
          stream.read(data, length);
          size = length;
        }

        return size;
      }

      MetaStream&
      seekg(streamsize pos)
      {
        assert_basic_istream();

        stream.seekg(streamBegin + pos);
        eofTrigger = false;
        return *this;
      }

      MetaStream&
      seekg(streamsize off, ios_base::seek_dir way)
      {
        assert_basic_istream();

        if(way == ios_base::beg)
          stream.seekg(streamBegin + off);
        else if(way == ios_base::end)
          stream.seekg(streamEnd - off);
        else
          stream.seekg(off, ios_base::cur);

        eofTrigger = false;
        return *this;
      }

      streamsize
      tellg() const
      {
        assert_basic_istream();
        return stream.tellg() - streamBegin;
      }

      bool
      eof() const
      {
        assert_basic_istream();
        streampos position = stream.tellg();
        if(not eofTrigger and position < streamEnd and position >= streamBegin)
          return false;
        else
          return true;
      }

    private:
      T& stream;
      T nullStream;
      streampos streamBegin;
      streampos streamEnd;
      const bool valid;
      bool eofTrigger;

      void assert_basic_istream() const
      {
        static_assert(
            is_base_of<
            basic_istream<typename T::char_type, typename T::traits_type>,
            T>::value,
            "class derives from basic_istream");
      }
  };
} /* namespace PstpFinder */
#endif /* METASTREAM_H_ */

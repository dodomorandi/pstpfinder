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

#include "MetaStream.h"
#include <sstream>

using namespace std;

namespace PstpFinder
{
  MetaStream::MetaStream() :
      inputStream(nullStream), valid(false), eofTrigger(false)
  {
    ;
  }

  MetaStream::MetaStream(ifstream& modifiableStream, streampos begin,
                         streampos end) :
      inputStream(modifiableStream), valid(true), eofTrigger(false)
  {
    if(begin == -1)
      streamBegin = inputStream.tellg();
    else
      streamBegin = begin;

    if(end == -1)
    {
      inputStream.seekg(0, ios_base::end);
      streamEnd = inputStream.tellg();
    }
    else
      streamEnd = end;

    inputStream.seekg(streamBegin, ios_base::beg);
  }

  MetaStream::MetaStream(const MetaStream& metaStream) :
      inputStream(metaStream.inputStream),
      valid(metaStream.valid),
      eofTrigger(false)
  {
    streamBegin = metaStream.streamBegin;
    streamEnd = metaStream.streamEnd;
  }

  bool
  MetaStream::eof() const
  {
    streampos position = inputStream.tellg();
    if(not eofTrigger and position < streamEnd and position >= streamBegin)
      return false;
    else
      return true;
  }

  template<typename T>
    void
    MetaStream::getFromStream(T& out)
    {
      if(not valid or eofTrigger)
        throw;

      if(eof())
        eofTrigger = true;
      streamoff currentPosition = inputStream.tellg();
      streampos finalPosition = currentPosition + sizeof(T);

      if(finalPosition > streamEnd)
      {
        stringstream tempStream(stringstream::in | stringstream::out);
        streamsize tempStreamSize = streamEnd - currentPosition;
        char *remainingData = new char[tempStreamSize];
        inputStream.read(remainingData, tempStreamSize);
        tempStream.write(remainingData, tempStreamSize);
        tempStream >> out;
        delete[] remainingData;
      }
      else
        inputStream >> out;
    }

  MetaStream&
  MetaStream::operator >>(string& out)
  {
    getFromStream<string>(out);
    return *this;
  }

  MetaStream&
  MetaStream::operator >>(bool& out)
  {
    getFromStream<bool>(out);
    return *this;
  }

  MetaStream&
  MetaStream::operator >>(char& out)
  {
    getFromStream<char>(out);
    return *this;
  }

  MetaStream&
  MetaStream::operator >>(unsigned char& out)
  {
    getFromStream<unsigned char>(out);
    return *this;
  }

  MetaStream&
  MetaStream::operator >>(short& out)
  {
    getFromStream<short>(out);
    return *this;
  }

  MetaStream&
  MetaStream::operator >>(unsigned short& out)
  {
    getFromStream<unsigned short>(out);
    return *this;
  }

  MetaStream&
  MetaStream::operator >>(int& out)
  {
    getFromStream<int>(out);
    return *this;
  }

  MetaStream&
  MetaStream::operator >>(unsigned int& out)
  {
    getFromStream<unsigned int>(out);
    return *this;
  }

  MetaStream&
  MetaStream::operator >>(long& out)
  {
    getFromStream<long>(out);
    return *this;
  }

  MetaStream&
  MetaStream::operator >>(unsigned long& out)
  {
    getFromStream<unsigned long>(out);
    return *this;
  }

  MetaStream&
  MetaStream::operator >>(void*& out)
  {
    getFromStream<void*>(out);
    return *this;
  }

  MetaStream&
  MetaStream::operator >>(float& out)
  {
    getFromStream<float>(out);
    return *this;
  }

  MetaStream&
  MetaStream::operator >>(double& out)
  {
    getFromStream<double>(out);
    return *this;
  }

  MetaStream&
  MetaStream::operator >>(long double& out)
  {
    getFromStream<long double>(out);
    return *this;
  }

  MetaStream&
  MetaStream::operator =(const MetaStream& metaStream)
  {
    if(&metaStream == this)
      return *this;

    this->~MetaStream();
    new (this) MetaStream(metaStream.inputStream, metaStream.streamBegin,
                          metaStream.streamEnd);
    return *this;
  }

  MetaStream&
  MetaStream::seekg(streamsize pos)
  {
    inputStream.seekg(streamBegin + pos);
    eofTrigger = false;
    return *this;
  }

  streamsize
  MetaStream::read(char* data, streamsize length)
  {
    streamsize size;
    if(not valid or eofTrigger)
      throw;

    streamoff currentPosition = inputStream.tellg();

    if(currentPosition + length > streamEnd)
    {
      inputStream.read(data, streamEnd - currentPosition);
      size = streamEnd - currentPosition;
    }
    else
    {
      inputStream.read(data, length);
      size = length;
    }
    if(eof())
      eofTrigger = true;

    return size;
  }

  MetaStream&
  MetaStream::seekg(streamsize off, ios_base::seek_dir way)
  {
    if(way == ios_base::beg)
      inputStream.seekg(streamBegin + off);
    else if(way == ios_base::end)
      inputStream.seekg(streamEnd - off);
    else
      inputStream.seekg(off, ios_base::cur);

    eofTrigger = false;
    return *this;
  }

  streamsize
  MetaStream::tellg() const
  {
    return inputStream.tellg() - streamBegin;
  }
} /* namespace Gromacs */

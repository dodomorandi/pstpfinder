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

using namespace std;

namespace PstpFinder
{
  MetaStream::MetaStream() :
      inputStream(nullStream), valid(false)
  {
    ;
  }

  MetaStream::MetaStream(ifstream& modifiableStream, streampos begin,
                         streampos end) :
      inputStream(modifiableStream), valid(true)
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
      valid(metaStream.valid)
  {
    streamBegin = metaStream.streamBegin;
    streamEnd = metaStream.streamEnd;
    currentPosition = metaStream.currentPosition;
  }

  bool
  MetaStream::eof() const
  {
    if(currentPosition < streamEnd and currentPosition >= streamBegin)
      return true;
    else
      return false;
  }

  template<typename T>
    void
    MetaStream::getFromStream(T& out)
    {
      if(not valid or eof())
        throw;

      inputStream >> out;
      currentPosition += sizeof(T);

      if(currentPosition > streamEnd)
      {
        inputStream.seekg(streamEnd);
        currentPosition = streamEnd;
      }
      else if(currentPosition < streamBegin)
      {
        inputStream.seekg(streamBegin);
        currentPosition = streamBegin;
      }
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
    return *this;
  }

  streamsize
  MetaStream::read(char* data, streamsize length)
  {
    streamsize size;

    if(currentPosition + length > streamEnd)
    {
      inputStream.read(data, streamEnd - currentPosition);
      size = streamEnd - currentPosition;
      currentPosition = streamEnd;
    }
    else
    {
      inputStream.read(data, length);
      size = length;
      currentPosition += length;
    }

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
      inputStream.seekg(off);

    return *this;
  }

  streamsize
  MetaStream::tellg() const
  {
    return inputStream.tellg() - streamBegin;
  }
} /* namespace Gromacs */

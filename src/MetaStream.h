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

#include "utils.h"

#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <type_traits>
#include <functional>

namespace PstpFinder
{
  enum class enumStreamType
  {
    STREAMTYPE_FIXED,
    STREAMTYPE_ADJUST
  };

  template<typename T, typename = void>
  class MetaStream;

  template<typename T>
  class MetaStream_Base : public T
  {
    public:
      typedef T stream_type;
      typedef typename T::off_type off_type;
      typedef typename T::char_type char_type;
      typedef typename T::traits_type traits_type;
      typedef typename T::int_type int_type;

      std::function<void(std::size_t)> callbackBeforeExtending;
      std::function<void()> callbackAfterExtending;
      std::function<void()> callbackClose;
      std::function<void()> callbackDestroy;

      template<template<typename, typename> class Stream>
      using base_stream = base_stream_to<Stream, T>;

      MetaStream_Base() :
          T(), streamType(enumStreamType::STREAMTYPE_FIXED), valid(false)
      {
      }

      /* FIXME: use iostreams move sematics, pleeeeease!!!
       * Damn it! It's 2012/03/27, and move and swap operations are still not
       * implemented inside ios_base!! Grrrr... I can't move streams, therefore
       * I need another solution...
       */
      // MetaStream_Base(T&& stream, enumStreamType streamType);

      MetaStream_Base(const std::string& filename,
                        std::ios_base::ios_base::openmode flags,
                        enumStreamType streamType) :
        T(filename, flags),
        streamType(streamType),
        valid(true)
      {
      }

      MetaStream_Base(T&& stream,
        enumStreamType streamType) :
        T(move(stream)), streamType(streamType),
        valid(true)
      {
      }

      ~MetaStream_Base()
      {
        if(T::is_open())
          close();
        if(callbackDestroy)
        {
          try
          {
            callbackDestroy();
          }
          catch(const std::bad_function_call& e)
          {
            std::cerr << "Warning: " << e.what() << std::endl;
          }
        }
      }

      // basic_istream related functions.
      template<typename U>
      MetaStream_Base&
      operator >>(U& out)
      {
        assert_basic_istream();
        assert(valid);

        std::streamoff currentPosition = T::tellg();
        off_type finalPosition;
        {
          U tempObject(out);
          T::operator>>(tempObject);
          finalPosition = T::tellg();
          T::seekg(currentPosition);
        }

        if(finalPosition > streamEnd)
        {
          std::stringstream tempStream(
              std::stringstream::in | std::stringstream::out);
          std::streamsize tempStreamSize = streamEnd - currentPosition;
          char *remainingData = new char[tempStreamSize];
          T::read(remainingData, tempStreamSize);
          tempStream.write(remainingData, tempStreamSize);
          tempStream >> out;
          delete[] remainingData;
          T::setstate(std::ios_base::eofbit);
        }
        else
          T::operator>>(out);

        return *this;
      }

      std::streamsize
      read(char_type* data, std::streamsize length)
      {
        assert_basic_istream();
        assert(valid);

        std::streamsize size;
        std::streamoff currentPosition = T::tellg();

        if(currentPosition + length > streamEnd)
        {
          T::read(data, streamEnd - currentPosition);
          size = streamEnd - currentPosition;
          T::setstate(std::ios_base::eofbit);
        }
        else
        {
          T::read(data, length);
          size = length;
        }

        return size;
      }

      MetaStream_Base&
      seekg(std::streamsize pos)
      {
        assert_basic_istream();

        T::seekg(streamBegin + pos);
        return *this;
      }

      MetaStream_Base&
      seekg(std::streamsize off, std::ios_base::seek_dir way)
      {
        assert_basic_istream();

        if(way == std::ios_base::beg)
          T::seekg(streamBegin + off);
        else if(way == std::ios_base::end)
          T::seekg(streamEnd - off);
        else
          T::seekg(off, std::ios_base::cur);

        return *this;
      }

      std::streamsize
      tellg()
      {
        assert_basic_istream();
        return T::tellg() - streamBegin;
      }

      int_type
      peek()
      {
        assert_basic_istream();
        off_t position { T::tellg() };
        if(position >= streamEnd)
        {
          T::setstate(std::ios_base::eofbit);
          return std::char_traits<char_type>::eof();
        }
        else
          return T::peek();
      }

      // basic_ostream related functions.
      template<typename U>
      MetaStream_Base&
      operator <<(U in)
      {
        return dumpStream(in);
      }

      MetaStream_Base&
      operator <<(base_stream<std::basic_ostream>&(*pf)
                  (base_stream<std::basic_ostream>&))
      {
        return dumpStream(pf);
      }

      MetaStream_Base&
      operator <<(std::ios_base&(*pf)(std::ios_base&))
      {
        return dumpStream(pf);
      }

      MetaStream_Base&
      operator <<(std::basic_ios<char_type, traits_type>&(*pf)
                  (std::basic_ios<char_type, traits_type>&))
      {
        return dumpStream(pf);
      }

      template<typename U>
      inline MetaStream_Base&
      dumpStream(U in)
      {
        assert_basic_ostream();
        assert(valid);

        std::streamoff currentPosition = T::tellp();
        off_type finalPosition = currentPosition;
        static std::basic_stringstream<char_type> tempBuffer;
        {
          tempBuffer << in;
          std::streamoff tempOffset(tempBuffer.tellp());
          if(tempOffset > 0)
          {
            currentPosition += tempOffset;
            tempBuffer.str("");
          }
        }

        if(finalPosition > streamEnd)
        {
          if(streamType == enumStreamType::STREAMTYPE_FIXED)
            throw;

          if(callbackBeforeExtending)
          {
            try
            {
              callbackBeforeExtending(finalPosition - streamEnd);
            }
            catch(const std::bad_function_call& e)
            {
              std::cerr << "Warning: " << e.what() << std::endl;
            }
          }

          static_cast<T&>(*this) << in;
          streamEnd = T::tellp();

          if(callbackAfterExtending)
          {
            try
            {
              callbackAfterExtending();
            }
            catch(const std::bad_function_call& e)
            {
              std::cerr << "Warning: " << e.what() << std::endl;
            }
          }
        }
        else
          static_cast<T&>(*this) << in;

        return *this;
      }

      MetaStream_Base&
      write(const char_type* data, std::streamsize length)
      {
        assert_basic_ostream();
        assert(valid);

        std::streamoff currentPosition = T::tellp();
        if(currentPosition + length > streamEnd)
        {
          std::streamsize writtenBytes;
          if(streamEnd - currentPosition >= 0)
              writtenBytes = streamEnd - currentPosition;
          else
              writtenBytes = 0;
          std::streamsize remainingBytes { length - writtenBytes };
          T::write(data, writtenBytes);

          if(streamType == enumStreamType::STREAMTYPE_FIXED)
            throw;

          if(callbackBeforeExtending)
          {
            try
            {
              callbackBeforeExtending(remainingBytes);
            }
            catch(const std::bad_function_call& e)
            {
              std::cerr << "Warning: " << e.what() << std::endl;
            }
          }

          T::write(data + writtenBytes, remainingBytes);
          streamEnd += remainingBytes;

          if(callbackAfterExtending)
          {
            try
            {
              callbackAfterExtending();
            }
            catch(const std::bad_function_call& e)
            {
              std::cerr << "Warning: " << e.what() << std::endl;
            }
          }
        }
        else
          T::write(data, length);

        return *this;
      }

      MetaStream_Base&
      seekp(std::streamsize pos)
      {
        assert_basic_ostream();

        T::seekp(streamBegin + pos);
        return *this;
      }

      MetaStream_Base&
      seekp(std::streamsize off, std::ios_base::seek_dir way)
      {
        assert_basic_ostream();

        if(way == std::ios_base::beg)
          T::seekp(streamBegin + off);
        else if(way == std::ios_base::end)
          T::seekp(streamEnd - off);
        else
          T::seekp(off, std::ios_base::cur);

        return *this;
      }

      std::streamsize
      tellp()
      {
        assert_basic_ostream();
        return T::tellp() - streamBegin;
      }

      void close()
      {
        if(T::is_open())
        {
          if(callbackClose)
          {
            try
            {
              callbackClose();
            }
            catch(const std::bad_function_call& e)
            {
              std::cerr << "Warning: " << e.what() << std::endl;
            }
          }
        }
        T::close();
      }

    protected:
      off_type streamBegin;
      off_type streamEnd;
      const enumStreamType streamType;

    private:
      const bool valid;

      inline void
      assert_basic_istream() const
      {
        static_assert(is_stream_base_of<std::basic_istream, T>::value,
                      "class derives from basic_istream");
      }

      inline void
      assert_basic_ostream() const
      {
        static_assert(is_stream_base_of<std::basic_ostream, T>::value,
                      "class derives from basic_ostream");
      }

      inline void
      assert_basic_stream() const
      {
        static_assert(is_stream_base_of<std::basic_istream, T>::value or
          is_stream_base_of<std::basic_ostream, T>::value,
          "class derives from basic_istream or basic_ostream");
      }
  };

  template<typename T>
  class MetaStream<T,
        typename std::enable_if<
            is_stream_base_of<std::basic_istream, T>::value and
    not is_stream_base_of<std::basic_ostream, T>::value>::type>
    : public MetaStream_Base<T>
  {
    public:
      typedef typename T::off_type off_type;

      MetaStream() : Base() {}
      MetaStream(const MetaStream& metaStream) : Base(metaStream) {}

      MetaStream(T&& stream, enumStreamType streamType =
                   enumStreamType::STREAMTYPE_FIXED,
                 off_type begin = -1, off_type end = -1) :
            Base(move(stream), streamType)
      {
        init(begin, end, streamType);
      }

      MetaStream(const std::string& filename, std::ios_base::openmode flags,
                 enumStreamType streamType = enumStreamType::STREAMTYPE_FIXED,
                 off_type begin = -1, off_type end = -1) :
            Base(filename, flags, streamType)
      {
        init(begin, end, streamType);
      }

    private:
      typedef MetaStream_Base<T> Base;
      void init(off_type begin, off_type end, enumStreamType streamType)
      {
        // FIXME: Unused streamType during initialization
        (void) streamType;
        assert(streamType == enumStreamType::STREAMTYPE_FIXED);

        if(begin == -1)
          Base::streamBegin = T::tellg();
        else
          Base::streamBegin = begin;

        if(end == -1)
        {
          T::seekg(0, std::ios_base::end);
          Base::streamEnd = T::tellg();
        }
        else
          Base::streamEnd = end;

        if(Base::streamEnd < Base::streamBegin)
          Base::streamEnd = Base::streamBegin;

        T::seekg(Base::streamBegin, std::ios_base::beg);
      }
  };

  template<typename T>
  class MetaStream<T,
        typename std::enable_if<
            not is_stream_base_of<std::basic_istream, T>::value and
    is_stream_base_of<std::basic_ostream, T>::value>::type>
    : public MetaStream_Base<T>
  {
    public:
      typedef typename T::off_type off_type;

      MetaStream() : Base() {}
      MetaStream(const MetaStream& metaStream) : Base(metaStream) {}

      MetaStream(T&& stream, enumStreamType streamType =
                   enumStreamType::STREAMTYPE_ADJUST,
                 off_type begin = -1, off_type end = -1) :
            Base(move(stream), streamType)
      {
        init(begin, end, streamType);
      }

      MetaStream(const std::string& filename, std::ios_base::openmode flags,
                 enumStreamType streamType = enumStreamType::STREAMTYPE_ADJUST,
                 off_type begin = -1, off_type end = -1) :
            Base(filename, flags, streamType)
      {
        init(begin, end, streamType);
      }

      static MetaStream_Base<T>&
      endl(MetaStream_Base<T>& stream)
      {
        std::endl(static_cast<T&>(stream));
        return stream;
      }

    private:
      typedef MetaStream_Base<T> Base;
      void init(off_type begin, off_type end, enumStreamType streamType)
      {
        // FIXME: Unused streamType during initialization
        (void) streamType;
        assert(
            not (streamType == enumStreamType::STREAMTYPE_FIXED and
                end == -1));

        if(begin == -1)
          Base::streamBegin = T::tellp();
        else
          Base::streamBegin = begin;

        if(end == -1)
        {
          T::seekp(0, std::ios_base::end);
          Base::streamEnd = T::tellp();
        }
        else
          Base::streamEnd = end;

        if(Base::streamEnd < Base::streamBegin)
          Base::streamEnd = Base::streamBegin;

        T::seekp(Base::streamBegin, std::ios_base::beg);
      }
  };

  template<typename T>
  class MetaStream<T,
        typename std::enable_if<
          is_stream_base_of<std::basic_istream,T>::value and
          is_stream_base_of<std::basic_ostream, T>::value>::type>
    : public MetaStream_Base<T>
  {
    public:
      typedef typename T::off_type off_type;

      MetaStream() : Base() {}
      MetaStream(const MetaStream& metaStream) : Base(metaStream) {}

      MetaStream(T&& stream, enumStreamType streamType =
                   enumStreamType::STREAMTYPE_ADJUST,
                 off_type begin = -1, off_type end = -1) :
            Base(move(stream), streamType)
      {
        init(begin, end, streamType);
      }

      MetaStream(const std::string& filename, std::ios_base::openmode flags,
                 enumStreamType streamType = enumStreamType::STREAMTYPE_ADJUST,
                 off_type begin = -1, off_type end = -1) :
            Base(filename, flags, streamType)
      {
        init(begin, end, streamType);
      }

    private:
      typedef MetaStream_Base<T> Base;
      void init(off_type begin, off_type end, enumStreamType streamType)
      {
        // TODO: Handling of different kind of streamType
        (void) streamType;

        if(begin == -1)
        {
          off_type tellg { T::tellg() };
          off_type tellp { T::tellp() };

          assert(tellg != -1 or tellp != -1);
          if(tellg == -1)
            Base::streamBegin = tellp;
          else if(tellp == -1)
            Base::streamBegin = tellg;
          else
            Base::streamBegin = std::min(tellg, tellp);
        }
        else
          Base::streamBegin = begin;

        if(end == -1)
        {
          T::seekg(0, std::ios_base::end);
          Base::streamEnd = T::tellg();
        }
        else
          Base::streamEnd = end;

        if(Base::streamEnd < Base::streamBegin)
          Base::streamEnd = Base::streamBegin;

        T::seekg(Base::streamBegin, std::ios_base::beg);
      }
  };

} /* namespace PstpFinder */
#endif /* METASTREAM_H_ */

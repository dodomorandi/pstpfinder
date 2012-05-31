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

using namespace std;

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
      typedef typename T::off_type off_type;
      typedef typename T::char_type char_type;
      typedef typename T::traits_type traits_type;
      typedef typename T::int_type int_type;

      function<void(size_t)> callbackBeforeExtending;
      function<void()> callbackAfterExtending;
      function<void()> callbackClose;
      function<void()> callbackDestroy;

      /*
       * FIXME
       * Can't use this... g++ 4.6.1 still doesn't support template aliases.
       * Waiting for 4.7 on Ubuntu/Mint
       *
      template<typename Stream>
        using base_stream = Stream<typename T::char_type,
        typename T::traits_type>;
       *
       * For now the implementation have been done using a #define (see above)
      */

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

      MetaStream_Base(const string& filename,
                        ios_base::ios_base::openmode flags,
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
          catch(const bad_function_call& e)
          {
            cerr << "Warning: " << e.what() << endl;
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

        streamoff currentPosition = T::tellg();
        off_type finalPosition;
        {
          U tempObject(out);
          T::operator>>(tempObject);
          finalPosition = T::tellg();
          T::seekg(currentPosition);
        }

        if(finalPosition > streamEnd)
        {
          stringstream tempStream(stringstream::in | stringstream::out);
          streamsize tempStreamSize = streamEnd - currentPosition;
          char *remainingData = new char[tempStreamSize];
          T::read(remainingData, tempStreamSize);
          tempStream.write(remainingData, tempStreamSize);
          tempStream >> out;
          delete[] remainingData;
          T::setstate(ios_base::eofbit);
        }
        else
          T::operator>>(out);

        return *this;
      }

      streamsize
      read(char_type* data, streamsize length)
      {
        assert_basic_istream();
        assert(valid);

        streamsize size;
        streamoff currentPosition = T::tellg();

        if(currentPosition + length > streamEnd)
        {
          T::read(data, streamEnd - currentPosition);
          size = streamEnd - currentPosition;
          T::setstate(ios_base::eofbit);
        }
        else
        {
          T::read(data, length);
          size = length;
        }

        return size;
      }

      MetaStream_Base&
      seekg(streamsize pos)
      {
        assert_basic_istream();

        T::seekg(streamBegin + pos);
        return *this;
      }

      MetaStream_Base&
      seekg(streamsize off, ios_base::seek_dir way)
      {
        assert_basic_istream();

        if(way == ios_base::beg)
          T::seekg(streamBegin + off);
        else if(way == ios_base::end)
          T::seekg(streamEnd - off);
        else
          T::seekg(off, ios_base::cur);

        return *this;
      }

      streamsize
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
          T::setstate(ios_base::eofbit);
          return char_traits<char_type>::eof();
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
      operator <<(base_stream(basic_ostream, T)&(*pf)
                  (base_stream(basic_ostream, T)&))
      {
        return dumpStream(pf);
      }

      MetaStream_Base&
      operator <<(ios_base&(*pf)(ios_base&))
      {
        return dumpStream(pf);
      }

      MetaStream_Base&
      operator <<(basic_ios<char_type, traits_type>&(*pf)
                  (basic_ios<char_type, traits_type>&))
      {
        return dumpStream(pf);
      }

      template<typename U>
      inline MetaStream_Base&
      dumpStream(U in)
      {
        assert_basic_ostream();
        assert(valid);

        streamoff currentPosition = T::tellp();
        off_type finalPosition = currentPosition;
        static basic_stringstream<char_type> tempBuffer;
        {
          tempBuffer << in;
          streamoff tempOffset(tempBuffer.tellp());
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
            catch(const bad_function_call& e)
            {
              cerr << "Warning: " << e.what() << endl;
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
            catch(const bad_function_call& e)
            {
              cerr << "Warning: " << e.what() << endl;
            }
          }
        }
        else
          static_cast<T&>(*this) << in;

        return *this;
      }

      MetaStream_Base&
      write(const char_type* data, streamsize length)
      {
        assert_basic_ostream();
        assert(valid);

        streamoff currentPosition = T::tellp();
        if(currentPosition + length > streamEnd)
        {
          streamsize writtenBytes;
          if(streamEnd - currentPosition >= 0)
              writtenBytes = streamEnd - currentPosition;
          else
              writtenBytes = 0;
          streamsize remainingBytes { length - writtenBytes };
          T::write(data, writtenBytes);

          if(streamType == enumStreamType::STREAMTYPE_FIXED)
            throw;

          if(callbackBeforeExtending)
          {
            try
            {
              callbackBeforeExtending(remainingBytes);
            }
            catch(const bad_function_call& e)
            {
              cerr << "Warning: " << e.what() << endl;
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
            catch(const bad_function_call& e)
            {
              cerr << "Warning: " << e.what() << endl;
            }
          }
        }
        else
          T::write(data, length);

        return *this;
      }

      MetaStream_Base&
      seekp(streamsize pos)
      {
        assert_basic_ostream();

        T::seekp(streamBegin + pos);
        return *this;
      }

      MetaStream_Base&
      seekp(streamsize off, ios_base::seek_dir way)
      {
        assert_basic_ostream();

        if(way == ios_base::beg)
          T::seekp(streamBegin + off);
        else if(way == ios_base::end)
          T::seekp(streamEnd - off);
        else
          T::seekp(off, ios_base::cur);

        return *this;
      }

      streamsize
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
            catch(const bad_function_call& e)
            {
              cerr << "Warning: " << e.what() << endl;
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
        static_assert(is_base_of<base_stream(basic_istream, T), T>::value,
                      "class derives from basic_istream");
      }

      inline void
      assert_basic_ostream() const
      {
        static_assert(is_base_of<base_stream(basic_ostream, T), T>::value,
                      "class derives from basic_ostream");
      }

      inline void
      assert_basic_stream() const
      {
        static_assert(is_base_of<base_stream(basic_istream, T), T>::value or
                      is_base_of<base_stream(basic_ostream, T), T>::value,
                      "class derives from basic_istream or basic_ostream");
      }
  };

  template<typename T>
  class MetaStream<T, typename enable_if<
            is_base_of<base_stream(basic_istream, T), T>::value and
            not is_base_of<base_stream(basic_ostream, T), T>::value>::type>
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

      MetaStream(const string& filename, ios_base::openmode flags,
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
        assert(streamType == enumStreamType::STREAMTYPE_FIXED);

        if(begin == -1)
          Base::streamBegin = T::tellg();
        else
          Base::streamBegin = begin;

        if(end == -1)
        {
          T::seekg(0, ios_base::end);
          Base::streamEnd = T::tellg();
        }
        else
          Base::streamEnd = end;

        if(Base::streamEnd < Base::streamBegin)
          Base::streamEnd = Base::streamBegin;

        T::seekg(Base::streamBegin, ios_base::beg);
      }
  };

  template<typename T>
  class MetaStream<T, typename enable_if<
            not is_base_of<base_stream(basic_istream, T), T>::value and
            is_base_of<base_stream(basic_ostream, T), T>::value>::type>
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

      MetaStream(const string& filename, ios_base::openmode flags,
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
        assert(
            not (streamType == enumStreamType::STREAMTYPE_FIXED and
                end == -1));

        if(begin == -1)
          Base::streamBegin = T::tellp();
        else
          Base::streamBegin = begin;

        if(end == -1)
        {
          T::seekp(0, ios_base::end);
          Base::streamEnd = T::tellp();
        }
        else
          Base::streamEnd = end;

        if(Base::streamEnd < Base::streamBegin)
          Base::streamEnd = Base::streamBegin;

        T::seekp(Base::streamBegin, ios_base::beg);
      }
  };

  template<typename T>
  class MetaStream<T, typename enable_if<
            is_base_of<base_stream(basic_istream, T), T>::value and
            is_base_of<base_stream(basic_ostream, T), T>::value>::type>
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

      MetaStream(const string& filename, ios_base::openmode flags,
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
            Base::streamBegin = min(tellg, tellp);
        }
        else
          Base::streamBegin = begin;

        if(end == -1)
        {
          T::seekg(0, ios_base::end);
          Base::streamEnd = T::tellg();
        }
        else
          Base::streamEnd = end;

        if(Base::streamEnd < Base::streamBegin)
          Base::streamEnd = Base::streamBegin;

        T::seekg(Base::streamBegin, ios_base::beg);
      }
  };

} /* namespace PstpFinder */
#endif /* METASTREAM_H_ */

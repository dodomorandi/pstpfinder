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

#ifndef SERIALIZER_H_
#define SERIALIZER_H_

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "utils.h"

#include <iostream>
#include <cassert>
#include <array>
#include <type_traits>
#include <functional>
#include <algorithm>

using namespace std;

namespace PstpFinder
{
  template<typename Type, typename = void>
  struct has_serialize_member_function : false_type {};

  template<typename Type>
  class has_serialize_member_function<Type,
    typename enable_if<is_class<Type>::value>::type>
  {
      class yes { char m;};
      class no { yes m[2];};
      struct Test
      {
        void serialize(){}
      };
      struct Base : public Type, public Test {};
      template <typename T, T t>  class Helper{};
      template <typename U>
      static no test(U*, Helper<void (Test::*)(), &U::serialize>* = 0);
      static yes test(...);
   public:
      static const bool value = sizeof(yes) == sizeof(test((Base*)(0)));
  };

  template<typename Stream, typename = void>
  class Serializer;

  template<typename Stream>
  class Serializer<
      Stream,
      typename enable_if<
          is_base_of<base_stream(basic_istream, Stream), Stream>::value>::type>
  {
    public:
      typedef typename Stream::char_type char_type;
      Serializer(Stream& stream) :
        stream(stream) {}

      template<typename Serializable>
      typename enable_if<is_arithmetic<Serializable>::value, Serializer&>::type
      operator &(Serializable& serializable)
      {
        array<char_type, sizeof(Serializable)> buffer;
        stream.read(buffer.data(), sizeof(Serializable));
        char_type* pointer((char_type*)&serializable);

#ifdef PSTPFINDER_BIG_ENDIAN // Default little endian
        reverse(begin(buffer), end(buffer));
#endif

        for(auto byte(begin(buffer)); byte < end(buffer); byte++, pointer++)
          *pointer = *byte;

        return *this;
      }

      template<typename Serializable>
      typename enable_if<is_class<Serializable>::value, Serializer&>::type
      operator &(Serializable& serializable)
      {
        size_t sz;
        *this & sz;

        char_type* buffer = new char_type[sz];
        stream.read(buffer, sz);
        stringstream sstream;
        sstream.write(buffer, sz);
        sstream >> serializable;
        delete[] buffer;

        return *this;
      }

      template<typename Output>
      typename enable_if<
            has_serialize_member_function<Output>::value,
            Serializer&>::type
      operator >>(Output& output)
      {
        output.serialize(*this);
        return *this;
      }

      template<typename Output>
      typename enable_if<
            not has_serialize_member_function<Output>::value,
            Serializer&>::type
      operator >>(Output& output)
      {
        *this & output;
        return *this;
      }

    private:
      Stream& stream;
  };

  template<typename Stream>
  class Serializer<
      Stream,
      typename enable_if<
          is_base_of<base_stream(basic_ostream, Stream), Stream>::value>::type>
  {
    public:
      typedef typename Stream::char_type char_type;
      Serializer(Stream& stream) :
        stream(stream) {}

      template<typename Serializable>
      typename enable_if<is_arithmetic<Serializable>::value, Serializer&>::type
      operator &(Serializable& serializable)
      {
        array<char_type, sizeof(Serializable)> buffer;
        char_type* pointer((char_type*)&serializable);

        for(auto byte(begin(buffer)); byte < end(buffer); byte++, pointer++)
          *byte = *pointer;

#ifdef PSTPFINDER_BIG_ENDIAN // Default little endian
        reverse(begin(buffer), end(buffer));
#endif
        stream.write(buffer.data(), sizeof(Serializable));

        return *this;
      }

      template<typename Serializable>
      typename enable_if<is_class<Serializable>::value, Serializer&>::type
      operator &(Serializable& serializable)
      {
        stringstream sstream;
        sstream << serializable;

        string buffer(sstream.str());
        size_t sz(buffer.size());
        *this & sz;
        stream << buffer;

        return *this;
      }

      template<typename Input>
      typename enable_if<
          has_serialize_member_function<Input>::value,
          Serializer&>::type
      operator <<(Input& input)
      {
        input.serialize(*this);
        return *this;
      }

      template<typename Input>
      typename enable_if<
          not has_serialize_member_function<Input>::value,
          Serializer&>::type
      operator <<(Input& input)
      {
        *this & input;
        return *this;
      }

    private:
      Stream& stream;
  };
} /* namespace PstpFinder */
#endif /* SERIALIZER_H_ */
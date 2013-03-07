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
#include "container_traits.h"

#include <iostream>
#include <cassert>
#include <array>
#include <type_traits>
#include <functional>
#include <algorithm>
#include <vector>
#include <sstream>

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
  class SerializerHelper
  {
    public:
      template<typename Serializable>
      static size_t getSerializedSize(Serializable serializable)
      {
        stringstream buffer;
        Serializer<stringstream> serializer(buffer);
        serializer << serializable;

        return static_cast<const size_t>(buffer.tellp());
      }

    protected:
      typedef typename Stream::char_type char_type;
      Stream& stream;

      SerializerHelper(Stream& stream) :
        stream(stream) {}

      template<typename Serializable>
      typename enable_if<is_arithmetic<Serializable>::value,
        const SerializerHelper&>::type
      serializeData(Serializable serializable) const
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
      typename enable_if<is_same<Serializable, string>::value,
        const SerializerHelper&>::type
      serializeData(const Serializable& serializable) const
      {
        serializeData(serializable.size());
        stream << serializable;

        return *this;
      }

      template<typename Serializable>
      typename enable_if<has_serialize_member_function<Serializable>::value,
        const SerializerHelper&>::type
      serializeData(Serializable& serializable) const
      {
        serializable.serialize(*static_cast<const Serializer<Stream>*>(this));
        return *this;
      }

      template<typename Serializable>
      typename enable_if<is_std_container<Serializable>::value,
        const SerializerHelper&>::type
      serializeData(Serializable& serializableVector) const
      {
        serializeContainer(serializableVector);

        return *this;
      }

      template<typename Serializable>
      typename enable_if<is_arithmetic<Serializable>::value,
        const SerializerHelper&>::type
      deserializeData(Serializable& serializable) const
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
      typename enable_if<is_same<Serializable, string>::value,
        const SerializerHelper&>::type
      deserializeData(Serializable& serializable) const
      {
        size_t sz;
        deserializeData(sz);

        serializable.resize(sz);
        char_type* buffer = new char_type[sz];
        stream.read(buffer, sz);
        serializable.assign(buffer, sz);
        return *this;
      }

      template<typename Serializable>
      typename enable_if<has_serialize_member_function<Serializable>::value,
        const SerializerHelper&>::type
      deserializeData(Serializable& serializable) const
      {
        serializable.serialize(*static_cast<const Serializer<Stream>*>(this));

        return *this;
      }

      template<typename Serializable>
      typename enable_if<is_std_container<Serializable>::value,
        const SerializerHelper&>::type
      deserializeData(Serializable& serializableVector) const
      {
        deserializeContainer(serializableVector);

        return *this;
      }

    private:
      template<typename Container>
      inline void
      serializeContainer(Container& container) const
      {
        serializeData(container.size());
        for(auto& serializable : container)
          serializeData(serializable);
      }

      template<typename Container>
      inline void
      deserializeContainer(Container& container) const
      {
        typedef typename Container::value_type containerType;
        size_t containerSize;
        deserializeData(containerSize);

        container.reserve(containerSize);
        for(size_t i(0); i < containerSize; i++)
        {
          containerType serializable;
          deserializeData(serializable);
          container.push_back(move(serializable));
        }
      }
  };

  template<typename Stream>
  class Serializer<
      Stream,
      typename enable_if<
          is_base_of<base_stream(basic_istream, Stream), Stream>::value and
          not is_base_of<base_stream(basic_ostream, Stream), Stream>::value>
        ::type> : public SerializerHelper<Stream>
  {
    public:
      typedef typename Stream::char_type char_type;
      Serializer(Stream& stream) :
        SerializerHelper<Stream>(stream) {}

      template<typename Serializable>
      typename enable_if<is_arithmetic<Serializable>::value, const Serializer&>::type
      operator &(Serializable& serializable) const
      {
        SerializerHelper<Stream>::deserializeData(serializable);
        return *this;
      }

      template<typename Serializable>
      typename enable_if<is_class<Serializable>::value, const Serializer&>::type
      operator &(Serializable& serializable) const
      {
        SerializerHelper<Stream>::deserializeData(serializable);
        return *this;
      }

      template<typename Output>
      typename enable_if<
            has_serialize_member_function<Output>::value,
            const Serializer&>::type
      operator >>(Output& output) const
      {
        output.serialize(*this);
        return *this;
      }

      template<typename Output>
      typename enable_if<
            not has_serialize_member_function<Output>::value,
            const Serializer&>::type
      operator >>(Output& output) const
      {
        *this & output;
        return *this;
      }
  };

  template<typename Stream>
  class Serializer<
      Stream,
      typename enable_if<
          not is_base_of<base_stream(basic_istream, Stream), Stream>::value and
          is_base_of<base_stream(basic_ostream, Stream), Stream>::value>
        ::type> : public SerializerHelper<Stream>
  {
    public:
      typedef typename Stream::char_type char_type;
      Serializer(Stream& stream) :
        SerializerHelper<Stream>(stream) {}

      template<typename Serializable>
      typename enable_if<is_arithmetic<Serializable>::value,
        const Serializer&>::type
      operator &(Serializable serializable) const
      {
        SerializerHelper<Stream>::serializeData(serializable);
        return *this;
      }

      template<typename Serializable>
      typename enable_if<is_class<Serializable>::value,
        const Serializer&>::type
      operator &(Serializable serializable) const
      {
        SerializerHelper<Stream>::serializeData(serializable);
        return *this;
      }

      template<typename Input>
      typename enable_if<
          has_serialize_member_function<Input>::value,
          const Serializer&>::type
      operator <<(Input input) const
      {
        input.serialize(*this);
        return *this;
      }

      template<typename Input>
      typename enable_if<
          not has_serialize_member_function<Input>::value,
          const Serializer&>::type
      operator <<(Input input) const
      {
        *this & input;
        return *this;
      }
  };

  template<typename Stream>
  class Serializer<
      Stream,
      typename enable_if<
          is_base_of<base_stream(basic_istream, Stream), Stream>::value and
          is_base_of<base_stream(basic_ostream, Stream), Stream>::value>
        ::type> : public SerializerHelper<Stream>
  {
    public:
      typedef typename Stream::char_type char_type;
      Serializer(Stream& stream) :
        SerializerHelper<Stream>(stream) {}

      template<typename Serializable>
      typename enable_if<is_arithmetic<Serializable>::value, const Serializer&>::type
      operator &(Serializable& serializable) const
      {
        if(mode == Mode::INPUT)
          SerializerHelper<Stream>::deserializeData(serializable);
        else
          SerializerHelper<Stream>::serializeData(serializable);
        return *this;
      }

      template<typename Serializable>
      typename enable_if<is_class<Serializable>::value, const Serializer&>::type
      operator &(Serializable& serializable) const
      {
        if(mode == Mode::INPUT)
          SerializerHelper<Stream>::deserializeData(serializable);
        else
          SerializerHelper<Stream>::serializeData(serializable);
        return *this;
      }

      template<typename Input>
      typename enable_if<
          has_serialize_member_function<Input>::value,
          const Serializer&>::type
      operator <<(Input input) const
      {
        mode = Mode::OUTPUT;
        input.serialize(*this);
        return *this;
      }

      template<typename Input>
      typename enable_if<
          not has_serialize_member_function<Input>::value,
          const Serializer&>::type
      operator <<(Input input) const
      {
        mode = Mode::OUTPUT;
        *this & input;
        return *this;
      }

      template<typename Output>
      typename enable_if<
            has_serialize_member_function<Output>::value,
            const Serializer&>::type
      operator >>(Output& output) const
      {
        mode = Mode::INPUT;
        output.serialize(*this);
        return *this;
      }

      template<typename Output>
      typename enable_if<
            not has_serialize_member_function<Output>::value,
            const Serializer&>::type
      operator >>(Output& output) const
      {
        mode = Mode::INPUT;
        *this & output;
        return *this;
      }

    private:
      enum class Mode
      {
          INPUT,
          OUTPUT
      };

      mutable Mode mode;
  };
} /* namespace PstpFinder */
#endif /* SERIALIZER_H_ */

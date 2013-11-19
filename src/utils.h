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

#ifndef UTILS_H_
#define UTILS_H_

#include <string>
#include <sstream>
#include <array>
#include <type_traits>

namespace PstpFinder
{
  #define base_stream(stream, T) stream<typename T::char_type, \
                                      typename T::traits_type>

  template<typename T>
  using remove_all = typename std::remove_cv<
      typename std::remove_reference<typename std::remove_pointer<
        typename std::remove_all_extents<T>::type>::type>::type>::type;

  template<template<typename, typename> class Stream, typename T,
            bool = std::is_base_of<std::ios_base, remove_all<T>>::value>
    struct is_stream_base_of_helper : public std::false_type {};

  template<template<typename, typename> class Stream,
            typename T>
  struct is_stream_base_of_helper<Stream, T, true>
  {
      typedef remove_all<T> __T;
      typedef typename std::is_base_of<Stream<typename __T::char_type,
                                           typename __T::traits_type>,
                                    __T>::type type;
      static constexpr bool value = type::value;
  };

  template<template<typename, typename> class Stream, typename T>
  struct is_stream_base_of : public is_stream_base_of_helper<Stream, T> {};

  bool exists(const std::string& filename);
  std::string file_extension(const std::string& filename);
  std::string change_extension(std::string filename,
                               const std::string& new_extension);


  template<std::size_t N>
  std::string
  array2string(const std::array<char, N> arr)
  {
    std::string str;
    str.reserve(N);
    for(const char& c : arr)
      if(c != '\0')
        str.push_back(c);
      else
        break;

    str.shrink_to_fit();
    return str;
  }
}

#endif /* UTILS_H_ */

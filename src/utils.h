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
#include <type_traits>

using namespace std;

namespace PstpFinder
{
  #define base_stream(stream, T) stream<typename T::char_type, \
                                      typename T::traits_type>

  template<typename T>
  using remove_all = typename remove_cv<
      typename remove_reference<typename remove_pointer<
        typename remove_all_extents<T>::type>::type>::type>::type;

  template<template<typename, typename> class Stream, typename T,
            bool = is_base_of<ios_base, remove_all<T>>::value>
    struct is_stream_base_of_helper : public false_type {};

  template<template<typename, typename> class Stream,
            typename T>
  struct is_stream_base_of_helper<Stream, T, true>
  {
      typedef remove_all<T> __T;
      typedef typename is_base_of<Stream<typename __T::char_type,
                                           typename __T::traits_type>,
                                    __T>::type type;
      static constexpr bool value = type::value;
  };

  template<template<typename, typename> class Stream, typename T>
  struct is_stream_base_of : public is_stream_base_of_helper<Stream, T> {};

  bool exists(const string& filename);
  string file_extension(const string& filename);
  string change_extension(string filename, const string& new_extension);
}
#endif /* UTILS_H_ */

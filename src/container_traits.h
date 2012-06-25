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

#ifndef CONTAINER_TRAITS_H_
#define CONTAINER_TRAITS_H_

#include <type_traits>
#include <array>
#include <vector>
#include <deque>
#include <forward_list>
#include <list>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <stack>
#include <queue>

namespace std
{
  template<typename Type>
  struct is_std_container : false_type {};

  template<typename T, size_t N>
  struct is_std_container<array<T, N>> : true_type {};

  template<typename T, typename Alloc>
  struct is_std_container<vector<T, Alloc>> : true_type {};

  template<typename T, typename Alloc>
  struct is_std_container<deque<T, Alloc>> : true_type {};

  template<typename T, typename Alloc>
  struct is_std_container<forward_list<T, Alloc>> : true_type {};

  template<typename T, typename Alloc>
  struct is_std_container<list<T, Alloc>> : true_type {};

  template<typename Key, typename Compare, typename Alloc>
  struct is_std_container<set<Key, Compare, Alloc>> : true_type {};

  template<typename Key, typename Compare, typename Alloc>
  struct is_std_container<multiset<Key, Compare, Alloc>> : true_type {};

  template<typename Key, typename T, typename Compare, typename Alloc>
  struct is_std_container<map<Key, T, Compare, Alloc>> : true_type {};

  template<typename Key, typename T, typename Compare, typename Alloc>
  struct is_std_container<multimap<Key, T, Compare, Alloc>> : true_type {};

  template<typename Key, typename Compare, typename Alloc>
  struct is_std_container<unordered_set<Key, Compare, Alloc>> : true_type {};

  template<typename Key, typename Compare, typename Alloc>
  struct is_std_container<unordered_multiset<Key, Compare, Alloc>> : true_type {};

  template<typename Key, typename T, typename Compare, typename Alloc>
  struct is_std_container<unordered_map<Key, T, Compare, Alloc>> : true_type {};

  template<typename Key, typename T, typename Compare, typename Alloc>
  struct is_std_container<unordered_multimap<Key, T, Compare, Alloc>> : true_type {};

  template<typename T, typename Sequence>
  struct is_std_container<stack<T, Sequence>> : true_type {};

  template<typename T, typename Sequence>
  struct is_std_container<queue<T, Sequence>> : true_type {};

  template<typename T, typename Sequence>
  struct is_std_container<priority_queue<T, Sequence>> : true_type {};
}

#endif /* CONTAINER_TRAITS_H_ */

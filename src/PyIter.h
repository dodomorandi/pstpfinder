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

#ifndef PYITER_H
#define PYITER_H

#include "Python.h"

#include <iterator>

// Utility for Python Object Iteration

namespace PstpFinder
{
namespace py
{

class python_iterator : public
    std::iterator<std::forward_iterator_tag, PyObject*>
{
  public:
    /* DefaultConstructible */
    python_iterator();
    /* Steals object reference */
    python_iterator(PyObject* obj);

    /* MoveConstructible */
    python_iterator(python_iterator&& iter) noexcept;

    /* CopyConstructible */
    python_iterator(const python_iterator& iter) noexcept;

    /* Move Assignable */
    python_iterator& operator=(python_iterator && iter) noexcept;

    /* Copy Assignable */
    python_iterator& operator=(const python_iterator& iter) noexcept;

    /* Destructible */
    ~python_iterator() noexcept;

    /* Iterator */
    reference operator*();
    python_iterator& operator++();

    /* Equality Comparable */
    bool operator==(const python_iterator& iter) noexcept;

    /* Input Iterator */
    bool operator!=(const python_iterator& iter) noexcept;
    python_iterator operator++(int);

  private:
    PyObject* it;
    PyObject* item;
};

class python_iterable
{
  public:
    python_iterable(PyObject* obj);

    python_iterable(python_iterable&& iterable) noexcept;
    python_iterable(const python_iterable& iterable) noexcept;

    ~python_iterable() noexcept;

    python_iterator begin() noexcept;
    python_iterator end() noexcept;

  private:
    PyObject* obj;
};

} /* namespace py */
} /* namespace PstpFinder */ 

#endif /* PYITER_H */

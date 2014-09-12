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

#include "config.h"
#ifdef HAVE_PYMOD_SADIC

#include "PyIter.h"

namespace PstpFinder
{
python_iterator::python_iterator() : it(nullptr), item(nullptr)
{}

python_iterator::python_iterator(PyObject* obj)
{
    if(PyIter_Check(obj))
        it = obj;
    else
    {
        it = PyObject_GetIter(obj);
        if(it == nullptr)
            throw;
    }

    item = PyIter_Next(it);
}

python_iterator::python_iterator(python_iterator&& iter) noexcept : it(
    std::move(iter.it)), item(std::move(iter.item))
{
}

python_iterator::python_iterator(const python_iterator& iter) noexcept : it(
    iter.it), item(iter.item)
{
    Py_XINCREF(item);
}

python_iterator&
python_iterator::operator=(python_iterator && iter) noexcept
{
    Py_XDECREF(item);

    item = nullptr;
    it = nullptr;

    std::swap(it, iter.it);
    std::swap(item, iter.item);

    return *this;
}

python_iterator&
python_iterator::operator=(const python_iterator& iter) noexcept
{
    it = iter.it;

    item = iter.item;
    Py_XINCREF(item);

    return *this;
}

python_iterator::~python_iterator() noexcept
{
    Py_XDECREF(item);
}

python_iterator::reference
python_iterator::operator*()
{
    if(item == nullptr)
        throw;

    return item;
}

python_iterator&
python_iterator::operator++()
{
    if(it == nullptr or item == nullptr)
        throw;

    Py_DECREF(item);
    item = PyIter_Next(it);

    return *this;
}

bool
python_iterator::operator==(const python_iterator& iter) noexcept
{
    return item == iter.item;
}

bool
python_iterator::operator!=(const python_iterator& iter) noexcept
{
    return item != iter.item;
}

python_iterator
python_iterator::operator++(int)
{
    python_iterator new_iter(*this);
    return ++new_iter;
}

/* ------------------------------------------------------------------------- */

python_iterable::python_iterable(PyObject* obj)
{
    this->obj = obj;
    Py_INCREF(this->obj);
}

python_iterable::python_iterable(python_iterable&& iterable) noexcept
{
    this->obj = iterable.obj;
}

python_iterable::python_iterable(const python_iterable& iterable) noexcept
{
    this->obj = iterable.obj;
    Py_INCREF(this->obj);
}

python_iterable::~python_iterable() noexcept
{
    Py_DECREF(this->obj);
}

python_iterator
python_iterable::begin() noexcept
{
    return python_iterator(obj);
}

python_iterator
python_iterable::end() noexcept
{
    return python_iterator();
}

} /* namespace PstpFinder */

#endif /* HAVE_PYMOD_SADIC */

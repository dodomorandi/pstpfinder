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

#ifndef PYUTILS_H
#define PYUTILS_H

#ifdef HAVE_PYTHON

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <Python.h>

namespace PstpFinder
{
namespace py
{
    template<typename... Args>
    inline PyObject*
    callMethod(PyObject* obj, const char* method, const char* arg1, Args&&... args)
    {
#if HAVE_PYTHON >= 3
        return PyObject_CallMethod(obj, method, arg1, std::forward<Args>(args)...);
#else
        return PyObject_CallMethod(obj, const_cast<char*>(method), const_cast<char*>(arg1), std::forward<Args>(args)...);
#endif
    }

    inline const char*
    toCString(PyObject* object)
    {
#if HAVE_PYTHON >= 3
        return reinterpret_cast<const char*>(PyUnicode_1BYTE_DATA(object));
#else
        return const_cast<const char*>(PyString_AsString(object));
#endif
    }
} /* namespace py */

} /* namespace PstpFinder */

#endif /* HAVE_PYTHON */

#endif /* PYUTILS_H */


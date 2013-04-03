AC_INIT([pstpfinder], [0.4], [enki@dodomorandi.no-ip.biz])
AM_INIT_AUTOMAKE([-Wall -Werror])
LT_INIT
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([src/pstpfinder.cpp])
AC_CONFIG_HEADERS([config.h])

# We don't want default AC_PROG_* behaviour
OLD_CXXFLAGS="$CXXFLAGS"
OLD_CPPFLAGS="$CPPFLAGS"
OLD_CFLAGS="$CFLAGS"

# Checks for programs.
AC_PROG_CXX
AX_CXX_CHECK_STD

# Restoring output variables
CXXFLAGS="$OLD_CXXFLAGS"
CPPFLAGS="$OLD_CPPFLAGS"
CFLAGS="$OLD_CFLAGS"

# Check for endianness
AX_ENDIAN

if test $CXX_SUPPORTS_CXX11 -eq 1; then
   CXXFLAGS="$CXXFLAGS -std=c++11"
elif test $CXX_SUPPORTS_CXX0X -eq 1; then
   CXXFLAGS="$CXXFLAGS -std=c++0x"
else
   AC_ERROR([cannot compile PSTP-Finder. A compiler compatible with C++11 features is required.])
fi

# Checks for libraries.
AX_GROMACS
AX_GROMACS_ANALYSIS
PKG_CHECK_MODULES([GTKMM30], [gtkmm-3.0], [
   CFLAGS="$CFLAGS $GTKMM30_CFLAGS"
   CXXFLAGS="$CXXFLAGS $GTKMM30_CFLAGS"
   LIBS="$LIBS $GTKMM30_LIBS"
   AC_DEFINE([GTKMM_MAJOR], [3], "GtkMM major version")
], [PKG_CHECK_MODULES([GTKMM24], [gtkmm-2.4], [
      CFLAGS="$CFLAGS $GTKMM24_CFLAGS"
      CXXFLAGS="$CXXFLAGS $GTKMM24_CFLAGS"
      LIBS="$LIBS $GTKMM24_LIBS"
      AC_DEFINE([GTKMM_MAJOR], [2], "GtkMM major version")
   ])
])

PKG_CHECK_MODULES([PYTHON], [python >= 2.6], [
   CFLAGS="$CFLAGS $PYTHON_CFLAGS"
   CXXFLAGS="$CXXFLAGS $PYTHON_CFLAGS"
   LIBS="$LIBS $PYTHON_LIBS"
   AC_DEFINE([HAVE_PYTHON], [1], "Python")
   HAVE_PYTHON=1
], [
   AC_DEFINE([HAVE_PYTHON], [0], "Python")
   HAVE_PYTHON=0
])

if test $HAVE_PYTHON -eq 1; then
   AX_PYTHON_MODULE([sadic])
   if test "$HAVE_PYMOD_SADIC" = "yes"; then
      AC_DEFINE([HAVE_PYMOD_SADIC], [1], "Sadic Python module")
      HAVE_PYMOD_SADIC=1
   else
      AC_DEFINE([HAVE_PYMOD_SADIC], [0], "Sadic Python module")
      HAVE_PYMOD_SADIC=0
   fi
else
   AC_DEFINE([HAVE_PYMOD_SADIC], [0], "Sadic Python module")
   HAVE_PYMOD_SADIC=0
fi

BOOST_REQUIRE([1.37])
BOOST_SYSTEM
if test $HAVE_PYTHON -eq 1 -a $HAVE_PYMOD_SADIC -eq 1; then
   BOOST_PYTHON
fi
LIBS="$LIBS $BOOST_SYSTEM_LIBS $BOOST_PYTHON_LIBS"
CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"

# Checks for header files.
AC_CHECK_HEADER([GL/gl.h])
AC_CHECK_HEADER([GL/glu.h])
AC_CHECK_HEADER([GL/glfw.h])

# Checks for typedefs, structures, and compiler characteristics.
AX_CONF_ARGS
AX_ENABLE_WARNINGS

# Checks for library functions.
AC_CHECK_LIB([GL], [glClear])
AC_CHECK_LIB([GLU], [gluBeginCurve])
AC_CHECK_LIB([glfw], [glfwInit])

# Checks if everything for GL context is ok
if test x$ac_cv_header_GL_gl_h = xyes -a x$ac_cv_header_GL_glu_h = xyes -a \
    x$ac_cv_header_GL_glfw_h = xyes -a x$ac_cv_lib_GL_glClear = xyes -a \
    x$ac_cv_lib_GLU_gluBeginCurve = xyes -a x$ac_cv_lib_glfw_glfwInit  = xyes
then
  GL_CONTEXT_OK=1
else
  GL_CONTEXT_OK=0
fi

# Define EXTRAS for bin_PROGRAMS
if test $GL_CONTEXT_OK -eq 1; then
  EXTRAS=marching_cubes
fi
AC_SUBST([EXTRAS])

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT

# Write out status
echo
echo -n "Compiling marching_cubes: "
if test $GL_CONTEXT_OK -eq 1; then
   echo "yes"
else
   echo "no (missing GL/GLU/GLFW package)"
fi

echo -n "Python support: "
if test $HAVE_PYTHON -eq 1; then
   echo "yes"
   
   echo -n "Boost Python support: "
   if test "x$BOOST_PYTHON_LIBS" != "x"; then
      echo "yes"
   else
      echo "no. Sadic support won't be built."
   fi
else
   echo "no. Sadic support can't be checked and it won't be built."
fi

echo -n "Sadic support: "
if test $HAVE_PYTHON -eq 1 -a "x$BOOST_PYTHON_LIBS" != "x" -a "$HAVE_PYMOD_SADIC" = "yes"; then
   echo "yes"
elif test $HAVE_PYTHON -eq 1 -a "x$BOOST_PYTHON_LIBS" != "x"; then
   echo "no"
   echo
   echo "You can download Sadic sources from http://www.sbl.unisi.it/?pg=download-software"
else
   echo "no"
fi

echo "Debug information: $ax_cv_enable_debug"
echo -n "Optimization level: "
if test $ax_cv_enable_debug = yes; then
  if test $ax_cv_enable_optimization = yes; then
    if test $ax_cv_optimization_level = default; then
       echo "3"
    else
       echo "$ax_cv_optimization_level"
    fi
  else
    echo "0"
  fi
else
  if test $ax_cv_enable_optimization \!= no; then
    if test $ax_cv_optimization_level = default; then
       echo "3"
    else
       echo "$ax_cv_optimization_level"
    fi
  else
    echo "0"
  fi
fi
echo

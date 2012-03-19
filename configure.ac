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
AC_SUBST([CPPFLAGS])
AC_SUBST([CXXFLAGS])
AC_SUBST([CFLAGS])

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
   AC_SUBST([CFLAGS])
   AC_SUBST([CXXFLAGS])
   AC_SUBST([LIBS])
   AC_DEFINE([GTKMM_MAJOR], [3], "GtkMM major version")
], [PKG_CHECK_MODULES([GTKMM24], [gtkmm-2.4], [
      CFLAGS="$CFLAGS $GTKMM24_CFLAGS"
      CXXFLAGS="$CXXFLAGS $GTKMM24_CFLAGS"
      LIBS="$LIBS $GTKMM24_LIBS"
      AC_SUBST([CFLAGS])
      AC_SUBST([CXXFLAGS])
      AC_SUBST([LIBS])
      AC_DEFINE([GTKMM_MAJOR], [2], "GtkMM major version")
   ])
])

PKG_CHECK_MODULES([PYTHON], [python >= 2.6], [
   CFLAGS="$CFLAGS $PYTHON_CFLAGS"
   CXXFLAGS="$CXXFLAGS $PYTHON_CFLAGS"
   LIBS="$LIBS $PYTHON_LIBS"
   AC_SUBST([CFLAGS])
   AC_SUBST([CXXFLAGS])
   AC_SUBST([LIBS])
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
   else
      AC_DEFINE([HAVE_PYMOD_SADIC], [0], "Sadic Python module")
   fi
else
   AC_DEFINE([HAVE_PYMOD_SADIC], [0], "Sadic Python module")
fi

BOOST_REQUIRE([1.37])
BOOST_IOSTREAMS
BOOST_SERIALIZATION
BOOST_FILESYSTEM
BOOST_SYSTEM
if test $HAVE_PYTHON -eq 1; then
   BOOST_PYTHON
fi
LIBS="$LIBS $BOOST_SERIALIZATION_LIBS $BOOST_IOSTREAMS_LIBS $BOOST_FILESYSTEM_LIBS $BOOST_SYSTEM_LIBS $BOOST_PYTHON_LIBS"
AC_SUBST([LIBS])

# Checks for header files.

CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
AC_SUBST([CPPFLAGS])

# Checks for typedefs, structures, and compiler characteristics.
AX_CONF_ARGS
AX_ENABLE_WARNINGS

# Checks for library functions.

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT

# Write out status
echo
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
echo

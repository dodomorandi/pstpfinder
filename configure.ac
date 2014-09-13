AC_INIT([pstpfinder], [0.4.3.2], [enki@dodomorandi.no-ip.biz])
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

# Check for endianness
AX_ENDIAN

if test $CXX_SUPPORTS_CXX11 -eq 1; then
   CXXFLAGS="$CXXFLAGS -std=c++11"
elif test $CXX_SUPPORTS_CXX0X -eq 1; then
   CXXFLAGS="$CXXFLAGS -std=c++0x"
else
   AC_MSG_ERROR([cannot compile PSTP-Finder. A compiler compatible with C++11 features is required.])
fi

# Checks for libraries.
AX_GROMACS
PKG_CHECK_MODULES([GTKMM30], [gtkmm-3.0], [
   CFLAGS="$CFLAGS $GTKMM30_CFLAGS"
   CXXFLAGS="$CXXFLAGS $GTKMM30_CFLAGS"
   LIBS="$LIBS $GTKMM30_LIBS"
   AC_SUBST([CFLAGS])
   AC_SUBST([CXXFLAGS])
   AC_SUBST([LIBS])
   AC_DEFINE([GTKMM_MAJOR], [3], [GtkMM major version])
], [PKG_CHECK_MODULES([GTKMM24], [gtkmm-2.4], [
      CFLAGS="$CFLAGS $GTKMM24_CFLAGS"
      CXXFLAGS="$CXXFLAGS $GTKMM24_CFLAGS"
      LIBS="$LIBS $GTKMM24_LIBS"
      AC_SUBST([CFLAGS])
      AC_SUBST([CXXFLAGS])
      AC_SUBST([LIBS])
      AC_DEFINE([GTKMM_MAJOR], [2], [GtkMM major version])
   ])
])

PKG_CHECK_MODULES([PYTHON2], [python2 >= 2.6], [
   AC_DEFINE([HAVE_PYTHON2], [1], "Python2")
   HAVE_PYTHON2=1
], [
   HAVE_PYTHON2=0
])

PKG_CHECK_MODULES([PYTHON3], [python3 >= 3.0], [
   AC_DEFINE([HAVE_PYTHON3], [1], "Python3")
   HAVE_PYTHON3=1
], [ 
   AC_DEFINE([HAVE_PYTHON3], [0], "Python3")
   HAVE_PYTHON3=0
])

if test $HAVE_PYTHON3 -eq 1; then
   AX_PYTHON_MODULE([sadic], [3])
   if test "$HAVE_PYMOD_SADIC" = "yes"; then
      CFLAGS="$CFLAGS $PYTHON3_CFLAGS"
      CXXFLAGS="$CXXFLAGS $PYTHON3_CFLAGS"
      LIBS="$LIBS $PYTHON3_LIBS"
      AC_SUBST([CFLAGS])
      AC_SUBST([CXXFLAGS])
      AC_SUBST([LIBS])
      AC_DEFINE([HAVE_PYMOD_SADIC], [1], "Sadic Python3 module")
   fi
elif test $HAVE_PYTHON2 -eq 1; then
   AX_PYTHON_MODULE([sadic], [2])
   if test "$HAVE_PYMOD_SADIC" = "yes"; then
      CFLAGS="$CFLAGS $PYTHON2_CFLAGS"
      CXXFLAGS="$CXXFLAGS $PYTHON2_CFLAGS"
      LIBS="$LIBS $PYTHON2_LIBS"
      AC_SUBST([CFLAGS])
      AC_SUBST([CXXFLAGS])
      AC_SUBST([LIBS])
      AC_DEFINE([HAVE_PYMOD_SADIC], [1], "Sadic Python2 module")
   fi
fi

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
if test $HAVE_PYTHON2 -eq 1 -o $HAVE_PYTHON3 -eq 1; then
   echo "yes"
else
   echo "no. Sadic support can't be checked and it won't be built."
fi

echo -n "Sadic support: "
if test "$HAVE_PYMOD_SADIC" = "yes"; then
   echo "yes"
elif test $HAVE_PYTHON2 -eq 1 -o $HAVE_PYTHON3 -eq 1; then
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

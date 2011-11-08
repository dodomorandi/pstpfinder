AC_INIT([pstpfinder], [0.1], [enki@dodomorandi.no-ip.biz])
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
PKG_CHECK_MODULES([GTKMM], [gtkmm-2.4], [
   CPPFLAGS="$CPPFLAGS $GTKMM_CFLAGS"
   LDFLAGS="$LDFLAGS $GTKMM_LIBS"
   AC_SUBST([CPPFLAGS])
   AC_SUBST([LDFLAGS])
])

BOOST_REQUIRE([1.34])
BOOST_IOSTREAMS
BOOST_SERIALIZATION
BOOST_THREADS
BOOST_FILESYSTEM
LIBS="$LIBS $BOOST_SERIALIZATION_LIBS $BOOST_THREAD_LIBS $BOOST_IOSTREAMS_LIBS $BOOST_FILESYSTEM_LIBS"
AC_SUBST([LIBS])

# Checks for header files.

CPPFLAGS=" -I../include/ $CPPFLAGS $BOOST_CPPFLAGS"
AC_SUBST([CPPFLAGS])

# Checks for typedefs, structures, and compiler characteristics.
AX_CONF_ARGS
AX_ENABLE_WARNINGS

# Checks for library functions.

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT

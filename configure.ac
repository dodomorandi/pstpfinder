AC_INIT([trapof], [0.1], [enki@dodomorandi.no-ip.biz])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
LT_INIT
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([src/trapof.cpp])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX

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
BOOST_THREADS
LIBS="$LIBS $BOOST_THREAD_LIBS"
AC_SUBST([LIBS])

# Checks for header files.

CPPFLAGS=" -I../include/ $CPPFLAGS $BOOST_CPPFLAGS"
AC_SUBST([CPPFLAGS])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT

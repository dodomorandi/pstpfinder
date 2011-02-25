AC_INIT([trapof], [0.1], [enki@dodomorandi.no-ip.biz])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([src/main.cpp])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX

# Checks for libraries.
AX_GROMACS
PKG_CHECK_MODULES([GTKMM], [gtkmm-2.4], [
   CPPFLAGS="$CPPFLAGS $GTKMM_CFLAGS"
   LDFLAGS="$LDFLAGS $GTKMM_LDFLAGS"
   AC_SUBST([CPPFLAGS])
   AC_SUBST([LDFLAGS])
])

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT

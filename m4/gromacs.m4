AC_DEFUN([AX_GROMACS],[
AC_ARG_VAR([GROMACS_INCLUDE_PATH], [Define the path for GROMACS headers])
AC_ARG_VAR([GROMACS_LIBRARY_PATH], [Define the path for GROMACS library])

_AX_HEADER_GROMACS([trnio.h],[$GROMACS_INCLUDE_PATH])
_AX_LIB_GROMACS([$GROMACS_LIBRARY_PATH])

CPPFLAGS="$CPPFLAGS -I$ax_cv_gromacs_headers"
if test "x$ax_cv_gromacs_library" != "x"; then
   LDFLAGS="$LDFLAGS -L$ax_cv_gromacs_library"
fi

AC_SUBST([CPPFLAGS])
AC_SUBST([LDFLAGS])

AC_CHECK_TYPE([struct output_env], [AC_DEFINE([GMX45], [1], [Gromacs major.minor version is 4.5])],,[[#include "oenv.h"]])

])

AC_DEFUN([AX_GROMACS_ANALYSIS],[

if test "x$ax_cv_gromacs_headers" = "x" ; then
  AX_GROMACS
fi

AC_CHECK_LIB([gmxana], [nsc_dclm_pbc],, AC_MSG_FAILURE([Missing GROMACS Analaysis library libgmxana. Is everything alright with your compilation?]))
])

AC_DEFUN([_AX_HEADER_GROMACS],[
OLD_CPPFLAGS="$CPPFLAGS"
for directory in $2 /usr/include/gromacs /usr/local/include/gromacs; do
   CPPFLAGS="$OLD_CFLAGS -I$directory"
   AC_CHECK_HEADER([$1], [ax_cv_gromacs_headers=$directory; break], [header_cache="ac_cv_header_$1"; unset ${header_cache/./_}])
done
if test "x$ax_cv_gromacs_headers" = "x" ; then
   AC_MSG_FAILURE([GROMACS headers needed to compile the program. Install it or set the environment variable GROMACS_INCLUDE_PATH])
fi
CPPFLAGS="$OLD_CPPFLAGS"
])

AC_DEFUN([_AX_LIB_GROMACS],[
AC_CHECK_LIB([gmx], [read_tpx],, [
   unset ac_cv_lib_gmx_read_tpx
   OLD_LDFLAGS="$LDFLAGS"
   found=0
   for directory in $1 $1/lib "" $ax_cv_gromacs_headers/../lib /lib /usr/lib /usr/local/lib; do
      if test "x$directory" != "x"; then
         LDFLAGS="$OLD_LDFLAGS -L$directory"
      else
         LDFLAGS="$OLD_LDFLAGS"
      fi
      AC_CHECK_LIB([gmx], [read_tpx], [ax_cv_gromacs_library=$directory; LIBS="$LIBS -lgmx"; found=1; break])
      unset ac_cv_lib_gmx_read_tpx
   done
   if test $found -eq 0 ; then
      AC_MSG_FAILURE([GROMACS libraries needed to compile the program. Install it or set the environment variable GROMACS_LIBRARY_PATH])
   fi
   LDFLAGS="$OLD_LDFLAGS"
])
])

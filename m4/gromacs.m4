AC_DEFUN([AX_GROMACS],[
AC_ARG_VAR([GROMACS_INCLUDE_PATH], [Define the path for GROMACS headers])
AC_ARG_VAR([GROMACS_LIBRARY_PATH], [Define the path for GROMACS library])

OLD_CPPFLAGS=$CPPFLAGS
if test x$GROMACS_LIBRARY_PATH != x; then
   LDFLAGS="-L$GROMACS_LIBRARY_PATH -Wl,-rpath,$GROMACS_LIBRARY_PATH -Wl,-rpath-link,$GROMACS_LIBRARY_PATH $LDFLAGS"
fi

path_found=0
for INCLUDE_PATH in $GROMACS_INCLUDE_PATH "" /usr/local/include; do
   if test x$INCLUDE_PATH != x; then
      CPPFLAGS="$OLD_CPPFLAGS -I$INCLUDE_PATH"
   fi

   AC_CHECK_HEADER([gromacs/trnio.h],[path_found=1; break])   
done

if test $path_found -eq 0; then
   AC_MSG_FAILURE([GROMACS headers needed to compile the program. Install it or set the environment variable GROMACS_INCLUDE_PATH])
fi

AC_CHECK_TYPE([struct output_env], [AC_DEFINE([GMX45], [1], [Gromacs major.minor version is 4.5])],[
   echo "Sorry, GROMACS versions before 4.5 still not implemented."
   exit 1
],[[#include <gromacs/oenv.h>]])

AC_CHECK_LIB([gmx],[read_tpx],,
   [AC_MSG_FAILURE([[GROMACS libraries needed to compile the program. ],
     [Install it or set the environment variable GROMACS_LIBRARY_PATH]])])

ax_cv_gromacs=ax_cv_gromacs

AC_SUBST([CPPFLAGS])
AC_SUBST([LDFLAGS])
])

AC_DEFUN([AX_GROMACS_ANALYSIS],[
if test x$ax_cv_gromacs = x; then
   AX_GROMACS()
fi

AC_CHECK_LIB([gmxana], [nsc_dclm_pbc],, AC_MSG_FAILURE([Missing GROMACS Analaysis library libgmxana. Is everything alright with your compilation?]))
])

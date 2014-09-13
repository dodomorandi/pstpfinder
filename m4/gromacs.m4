AC_DEFUN([AX_GROMACS],[
AC_ARG_VAR([GROMACS_PATH], [Define the prefix for GROMACS installation. If omitted '/usr/local' and '/usr' directories will be used.])
AC_ARG_VAR([GROMACS_INCLUDE_PATH], [Define the path for GROMACS headers. Default if '$GROMACS_PATH/include'.])
AC_ARG_VAR([GROMACS_LIBRARY_PATH], [Define the path for GROMACS library. Default if '$GROMACS_PATH/lib[32/64]'.])

OLD_CPPFLAGS=$CPPFLAGS
path_found=0
for INCLUDE_PATH in $GROMACS_INCLUDE_PATH $GROMACS_PATH/include "" /usr/local/include; do
   if test x$INCLUDE_PATH != x; then
      CPPFLAGS="$OLD_CPPFLAGS -I$INCLUDE_PATH"
   fi

   unset ac_cv_header_gromacs_trnio_h
   unset ac_cv_header_gromacs_fileio_trnio_h
   AC_CHECK_HEADER([gromacs/trnio.h],[path_found=1; gmx_ver=45; break])   

   AC_CHECK_HEADER([gromacs/fileio/trnio.h],[path_found=1; gmx_ver=50; break])
done

if test $path_found -eq 0; then
   AC_MSG_FAILURE([GROMACS headers needed to compile the program. Install it or set the environment variable GROMACS_PATH or GROMACS_INCLUDE_PATH])
fi

if test $gmx_ver -eq 45; then
   AC_CHECK_TYPE([struct output_env], [
      AC_DEFINE([GMXVER], [45], [Gromacs version is >= 4.5])
   ],[
      echo "Sorry, GROMACS versions before 4.5 still not implemented."
      exit 1
   ],[[#include <gromacs/oenv.h>]])
else
   AC_DEFINE([GMXVER], [50], [Gromacs version is >= 5.0])
fi

OLD_LDFLAGS=$LDFLAGS
path_found=0
for LIBRARY_PATH in $GROMACS_LIBRARY_PATH ${GROMACS_PATH}/lib ${GROMACS_PATH}/lib64 ${GROMACS_PATH}/lib32 ""; do
   if test x$LIBRARY_PATH != x; then
      LDFLAGS="$OLD_LDFLAGS -L$LIBRARY_PATH -Wl,-rpath,$LIBRARY_PATH -Wl,-rpath-link,$LIBRARY_PATH"
   fi
   
   if test $gmx_ver -le 45; then
      unset ac_cv_lib_gmx_read_tpx
      AC_CHECK_LIB([gmx],[read_tpx],[
         path_found=1
         LIBS="-lgmx $LIBS"
         AC_DEFINE([HAVE_LIBGMX],[1],"Have gmx library")
         break
      ])
   else
      unset ac_cv_lib_gromacs_read_tpx
      AC_CHECK_LIB([gromacs],[read_tpx],[
         path_found=1
         LIBS="-lgromacs $LIBS"
         AC_DEFINE([HAVE_LIBGROMACS],[1], "Have gromacs library")
         break
      ])
   fi
done

if test $path_found -eq 0; then
   AC_MSG_FAILURE([[GROMACS libraries needed to compile the program. ],
     [Install it or set the environment variable GROMACS_PATH or GROMACS_LIBRARY_PATH]])
fi

if test $gmx_ver -le 45; then
   AC_CHECK_LIB([gmxana], [nsc_dclm_pbc],, AC_MSG_FAILURE([Missing GROMACS Analaysis library libgmxana. Is everything alright with your compilation?]))
fi

AC_SUBST([CPPFLAGS])
AC_SUBST([LDFLAGS])
])

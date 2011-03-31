AC_DEFUN([AX_ENABLE_WARNINGS],[
   CPPFLAGS="-Wall ${CPPFLAGS}"
   AC_SUBST([CPPFLAGS])])

AC_DEFUN([AX_CONF_ARGS],[
   AC_ARG_ENABLE([debug],
      [AS_HELP_STRING([--enable-debug[=[yes|no]]],
         [enable debug or set to a state (default is no)])],
      [
         case x$enableval in
         xyes|xno)
            ax_cv_enable_debug=$enableval
            ;;
         x)
            ax_cv_enable_debug=yes
            ;;
         *)
            ax_cv_enable_debug=no
            ;;
         esac
      ],
      [ax_cv_enable_debug=no])
   
   AC_ARG_ENABLE([optimization],
      [AS_HELP_STRING([--disable-optimization[=[yes|no]]],
         [disable optimization (normally it would depend on --enable-debug)])],
      [
         case x$enableval in
         xyes | xno)
            ax_cv_disable_optimization=$enableval
            ;;
         x)
            ax_cv_disable_optimization=yes
            ;;
         *)
            ax_cv_disable_optimization=unspec
            ;;
         esac
      ],
      [ax_cv_disable_optimization=unspec])
   
   if test ${ax_cv_enable_debug} = yes; then
      CPPFLAGS="-g ${CPPFLAGS}"
      if test ${ax_cv_disable_optimization} = no; then
        CPPFLAGS="-O3 ${CPPFLAGS}"
      fi
   else
      if test ${ax_cv_disable_optimization} != yes; then
         CPPFLAGS="-O3 ${CPPFLAGS}"
      fi
   fi
   AC_SUBST([CPPFLAGS])
])

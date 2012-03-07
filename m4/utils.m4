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
      [AS_HELP_STRING([--enable-optimization[=[yes|no]]],
         [enable optimization (normally it would depend on --enable-debug)])],
      [
         case x$enableval in
         xyes | xy | xno | xn)
            ax_cv_enable_optimization=$enableval
            ;;
         x)
            ax_cv_enable_optimization=yes
            ;;
         *)
            ax_cv_enable_optimization=default
            ;;
         esac
      ],
      [ax_cv_enable_optimization=default]) 

   AC_ARG_WITH([optimization-level],
      [AS_HELP_STRING([--with-optimization-level=[[0|1|2|3]]],
         [defines the level of optimization (default is 3 if optimization is enabled, ignored if disabled)])],
      [
         case x$withval in
         x0 | x1 | x2 | x3)
            ax_cv_optimization_level=$withval
            ;;
         *)
            ax_cv_optimization_level=default
            ;;
         esac
      ],
      [ax_cv_optimization_level=default])
  
   if test x${ax_cv_enable_debug} = xyes; then
      ax_cv_old_debug_flags="-g"
      CPPFLAGS="-g ${CPPFLAGS}"
      if test `echo x${ax_cv_enable_optimization}|cut -c -2` = xy; then
         if test x${ax_cv_optimization_level} = xdefault; then
            CPPFLAGS="-O3 ${CPPFLAGS}"
            ax_cv_old_debug_flags="-O3 ${ax_cv_old_debug_flags}"
         else
            CPPFLAGS="-O${ax_cv_optimization_level} ${CPPFLAGS}"
            ax_cv_old_debug_flags="-O${ax_cv_optimization_level} ${ax_cv_old_debug_flags}"
         fi
      else
         CPPFLAGS="-O0 ${CPPFLAGS}"
         ax_cv_old_debug_flags="-O0 ${ax_cv_old_debug_flags}"
      fi
   else
      if test `echo x${ax_cv_enable_optimization}|cut -c -2` != xn; then
         if test x${ax_cv_optimization_level} = xdefault; then
            CPPFLAGS="-O3 ${CPPFLAGS}"
            ax_cv_old_debug_flags="-O3"
         else
            CPPFLAGS="-O${ax_cv_optimization_level} ${CPPFLAGS}"
            ax_cv_old_debug_flags="-O${ax_cv_optimization_level}"
         fi
      else
         CPPFLAGS="-O0 ${CPPFLAGS}"
         ax_cv_old_debug_flags="-O0"
      fi
   fi

   AC_SUBST([CPPFLAGS])
])

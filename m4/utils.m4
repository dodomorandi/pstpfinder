AC_DEFUN([AX_ENDIAN],[
   AC_LANG_PUSH([C])
   AC_MSG_CHECKING([system endianness])
   AC_RUN_IFELSE(
       [AC_LANG_PROGRAM([],[
             int main(){
                 unsigned short cafe = 0xCAFE;
                 unsigned char *p = (unsigned char*)&cafe;
                 if(*p == 0xFE)
                   return 0; // Little endian
                 else
                   return 1; // Big endian
             }
          ])],[
       AC_MSG_RESULT([little endian])
       AC_DEFINE([PSTPFINDER_LITTLE_ENDIAN],[],[Little endian])
   ],[
       AC_MSG_RESULT([big endian])
       AC_DEFINE([PSTPFINDER_BIG_ENDIAN],[],[Big endian])
   ])
      
   AC_LANG_POP([C])
])

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
      CPPFLAGS="-g ${CPPFLAGS}"
      if test `echo x${ax_cv_enable_optimization}|cut -c -2` = xy; then
         if test x${ax_cv_optimization_level} = xdefault; then
            CPPFLAGS="-O3 ${CPPFLAGS}"
         else
            CPPFLAGS="-O${ax_cv_optimization_level} ${CPPFLAGS}"
         fi
      else
         CPPFLAGS="-O0 ${CPPFLAGS}"
      fi
   else
      if test `echo x${ax_cv_enable_optimization}|cut -c -2` != xn; then
         if test x${ax_cv_optimization_level} = xdefault; then
            CPPFLAGS="-O3 -DNDEBUG ${CPPFLAGS}"
         else
            CPPFLAGS="-O${ax_cv_optimization_level} -DNDEBUG ${CPPFLAGS}"
         fi
      else
         CPPFLAGS="-O0 -DNDEBUG ${CPPFLAGS}"
      fi
   fi

   AC_SUBST([CPPFLAGS])
])

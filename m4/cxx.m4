# AX_CXX_CHECK_STD()
# ------------------
# Check for C++ standards support. Defines CXX_SUPPORTS_CXX11 and
# CXX_SUPPORTS_CXX0X, setting them to 1 or 0 if the compiler supports
# the options -std=c++11 and -std=c++0x.
# For now the check is made using initializer_list presence
AC_DEFUN([AX_CXX_CHECK_STD], [
  ax_is_cxx11_save_CPPFLAGS=CPPFLAGS
  CPPFLAGS="-std=c++11"
  AC_LANG_PUSH([C++])
  AC_MSG_CHECKING([if C++ compiler supports -std=c++11])
  AC_COMPILE_IFELSE(
    [AC_LANG_PROGRAM([#include <initializer_list>],
  	  [std::initializer_list<int> list;])],[
  	AC_MSG_RESULT([yes])
  	AC_SUBST([CXX_SUPPORTS_CXX11], [1])
  ],[
  	AC_MSG_RESULT([no])
  	AC_SUBST([CXX_SUPPORTS_CXX11], [0])
  ])
  	
  CPPFLAGS="-std=c++0x"
  AC_MSG_CHECKING([if C++ compiler supports -std=c++0x])
  AC_LANG_CONFTEST(
    [AC_LANG_SOURCE([int main(){return 0;}])])
  AC_COMPILE_IFELSE(
    [AC_LANG_PROGRAM([#include <initializer_list>],
     [std::initializer_list<int> list;])],[
    AC_MSG_RESULT([yes])
    AC_SUBST([CXX_SUPPORTS_CXX0X], [1])
  ],[
    AC_MSG_RESULT([no])
    AC_SUBST([CXX_SUPPORTS_CXX0X], [0])
  ])
  AC_LANG_POP([C++])
  CPPFLAGS=ax_is_cxx1_save_CPPFLAGS
])

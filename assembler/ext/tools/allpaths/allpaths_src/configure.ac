AC_INIT([AllPathsLG], [40819], [crdhelp@broadinstitute.org])
AM_INIT_AUTOMAKE
AM_PROG_LIBTOOL([disable-shared])
AC_PROG_CXX
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_HEADERS([src/config.h])

AC_LITTLE_ENDIAN_CHECK
AC_GXX_CHECK
AX_PTHREAD
AC_OPENMP
AC_OPENMP_CEHCK

CXXFLAGS="$CXXFLAGS -Wextra -Wall -ansi -pedantic -Wsign-promo -Woverloaded-virtual -Wendif-labels -Wno-unused -Wno-deprecated -Wno-long-long -Wno-parentheses -fno-nonansi-builtins -ftemplate-depth-50 -mieee-fp -fno-strict-aliasing -iquote . -ggdb3 -DNDEBUG"

OS_RELEASE="`uname -r`"
AC_SUBST(CXXFLAGS)
AC_SUBST(OS_RELEASE, [$OS_RELEASE])
AC_SUBST(VERSION, [$VERSION])
AC_OUTPUT

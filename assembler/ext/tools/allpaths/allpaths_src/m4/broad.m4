host_machine=`uname -m`
host_os=`uname -o`

AC_DEFUN([AC_LINUX_CHECK], [
    AC_MSG_CHECKING([for GNU/Linux platform.])
    AS_IF(
        [test "x$host_os" != "xGNU/Linux"], 
        [AC_MSG_BROAD_FAIL([Host must be GNU/Linux])], 
        []
)]) 

AC_DEFUN([AC_64BIT_CHECK], [
    AC_MSG_CHECKING([for 64-bit platform.])
    AS_IF([test "x$host_machine" != "xx86_64"],
        [AC_MSG_BROAD_FAIL([Host must a 64-bit x86 platform])], 
        [])
])

AC_DEFUN([AC_LITTLE_ENDIAN_CHECK], [
    AC_C_BIGENDIAN(AC_MSG_BROAD_FAIL([Host must be little endian]),
        [],
        AC_MSG_BROAD_FAIL([Host must be little endian]), 
        AC_MSG_BROAD_FAIL([Host must be little endian])) 
])

AC_DEFUN([AC_GXX_CHECK], [
    CXX_VER=`$CXX -dumpversion`
    AX_COMPARE_VERSION([$CXX_VER], [ge], [4.3.3],
    AC_MSG_RESULT([g++ version is >= 4.3.3... yes]),
    AC_MSG_BROAD_FAIL([You must compile this with g++ 4.3.3 or higher.]))])

AC_DEFUN([AC_OPENMP_CEHCK], [
    OLD_CXXFLAGS=$CXXFLAGS
    CXXFLAGS="$OPENMP_CFLAGS $CXXFLAGS"
    OLD_CFLAGS=$CFLAGS
    CFLAGS="$OPENMP_CFLAGS $CFLAGS"
    AC_MSG_CHECKING([validity of OpenMP configuation.])
    AC_LINK_IFELSE([#ifndef _OPENMP
     choke me
    #endif
    #include <omp.h>
    int main () { return omp_get_num_threads (); }],
    [],
    [AC_MSG_BROAD_FAIL([Your compiler must support OpenMP.])])
    CXXFLAGS=$OLD_CXXFLAGS
    CFLAGS=$OLD_CFLAGS
])

AC_DEFUN([AC_MSG_BROAD_FAIL], [
  AC_MSG_RESULT([Failed...])                    
  AC_MSG_NOTICE([Configure failed with the error message: $1])
  AC_MSG_NOTICE([])
  AC_MSG_NOTICE([Common error conditions in the build process are documented])
  AC_MSG_NOTICE([on our wiki page:])
  AC_MSG_NOTICE([http://www.broadinstitute.org/crd/wiki/index.php/Build_FAQ])
  AC_MSG_NOTICE([We also offer email support at crdhelp@broadinstitute.org])
  AC_MSG_FAILURE([Configure failed.])
])


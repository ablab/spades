#.rst:
# FindAltiVec
# --------
#
# Finds AltiVec support
#
# This module can be used to detect AltiVec support in a C compiler.  If
# the compiler supports AltiVec, the flags required to compile with
# AltiVec support are returned in variables for the different languages.
# The variables may be empty if the compiler does not need a special
# flag to support AltiVec.
#
# The following variables are set:
#
# ::
#
#    ALTIVEC_C_FLAGS - flags to add to the C compiler for AltiVec support
#    ALTIVEC_FOUND - true if AltiVec is detected
#
#=============================================================================

set(_ALTIVEC_REQUIRED_VARS)
set(CMAKE_REQUIRED_QUIET_SAVE ${CMAKE_REQUIRED_QUIET})
set(CMAKE_REQUIRED_QUIET ${ALTIVEC_FIND_QUIETLY})

# sample AltiVec source code to test
set(ALTIVEC_C_TEST_SOURCE
"
#include <altivec.h>
vector signed int vec_2sComp (vector signed int x)
{
    vector signed int one = (vector signed int) (1);
    x = vec_nor (x, x);
    x = vec_add (x, one);
    return x;
}
int main(void)
{
    vector signed int two = (vector signed int) (2);
    vector signed int result = vec_2sComp(two);
    return vec_extract(two,0);
}
")

# if these are set then do not try to find them again,
# by avoiding any try_compiles for the flags
if((DEFINED ALTIVEC_C_FLAGS) OR (DEFINED HAVE_ALTIVEC))
else()
  if(WIN32)
    set(ALTIVEC_C_FLAG_CANDIDATES
      #Empty, if compiler automatically accepts AltiVec
      " ")
  else()
    set(ALTIVEC_C_FLAG_CANDIDATES
      #Empty, if compiler automatically accepts AltiVec
      " "
      "-maltivec"
      "-faltivec"
    )
  endif()

  include(CheckCSourceCompiles)

  foreach(FLAG IN LISTS ALTIVEC_C_FLAG_CANDIDATES)
    set(SAFE_CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
    set(CMAKE_REQUIRED_FLAGS "${FLAG}")
    unset(HAVE_ALTIVEC CACHE)
    if(NOT CMAKE_REQUIRED_QUIET)
      message(STATUS "Try AltiVec C flag = [${FLAG}]")
    endif()
    check_c_source_compiles("${ALTIVEC_C_TEST_SOURCE}" HAVE_ALTIVEC)
    set(CMAKE_REQUIRED_FLAGS "${SAFE_CMAKE_REQUIRED_FLAGS}")
    if(HAVE_ALTIVEC)
      set(ALTIVEC_C_FLAGS_INTERNAL "${FLAG}")
      break()
    endif()
  endforeach()

  unset(ALTIVEC_C_FLAG_CANDIDATES)
  
  set(ALTIVEC_C_FLAGS "${ALTIVEC_C_FLAGS_INTERNAL}"
    CACHE STRING "C compiler flags for AltiVec intrinsics")
endif()

list(APPEND _ALTIVEC_REQUIRED_VARS ALTIVEC_C_FLAGS)

set(CMAKE_REQUIRED_QUIET ${CMAKE_REQUIRED_QUIET_SAVE})

if(_ALTIVEC_REQUIRED_VARS)
  include(FindPackageHandleStandardArgs)

  find_package_handle_standard_args(ALTIVEC
                                    REQUIRED_VARS ${_ALTIVEC_REQUIRED_VARS})

  mark_as_advanced(${_ALTIVEC_REQUIRED_VARS})

  unset(_ALTIVEC_REQUIRED_VARS)
else()
  message(SEND_ERROR "FindAltiVec requires C or CXX language to be enabled")
endif()


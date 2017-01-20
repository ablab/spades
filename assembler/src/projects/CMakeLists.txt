############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

add_subdirectory(spades)
add_subdirectory(hammer)
add_subdirectory(ionhammer)
add_subdirectory(dipspades)
add_subdirectory(corrector)
add_subdirectory(scaffold_correction)

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
  # Require at least gcc 4.8
  if (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.8)
    add_subdirectory(mph_test)
  endif()
else()
  add_subdirectory(mph_test)
endif()

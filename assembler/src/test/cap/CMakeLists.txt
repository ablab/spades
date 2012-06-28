############################################################################
# Copyright (c) 2011-2012 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(cap_test CXX)

include_directories(${CMAKE_SOURCE_DIR}/debruijn)
include_directories(${CMAKE_SOURCE_DIR}/cap)

add_executable(cap_test
 test.cpp)

target_link_libraries(cap_test ${COMMON_LIBRARIES} input)

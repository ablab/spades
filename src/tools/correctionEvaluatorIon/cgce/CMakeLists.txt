############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(cgce CXX)

include_directories(${CMAKE_SOURCE_DIR}/debruijn
                    ${CMAKE_SOURCE_DIR}/online_vis
                    ${CMAKE_SOURCE_DIR}/cap)

add_executable(cgce
  cgce.cpp)

target_link_libraries(cgce input mph_index yaml-cpp ${COMMON_LIBRARIES})

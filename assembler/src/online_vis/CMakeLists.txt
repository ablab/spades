############################################################################
# Copyright (c) 2011-2012 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(online_vis CXX)

add_executable(online_vis
 main.cpp)

target_link_libraries(online_vis ${COMMON_LIBRARIES} input)

include_directories(${CMAKE_SOURCE_DIR}/debruijn)

if (DEFINED static_build)
  set_target_properties(online_vis PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()


############################################################################
# Copyright (c) 2011-2012 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(online_vis CXX)

add_executable(online_vis main.cpp)

if (READLINE_FOUND)
 include_directories(${READLINE_INCLUDE_DIR})
endif()
include_directories(./drawing_commands)
include_directories(${CMAKE_SOURCE_DIR}/debruijn)

target_link_libraries(online_vis ${COMMON_LIBRARIES} ${READLINE_LIBRARY} input mph_index )

if (DEFINED static_build)
  set_target_properties(online_vis PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()


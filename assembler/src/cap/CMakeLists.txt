############################################################################
# Copyright (c) 2011-2012 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(cap CXX)

add_executable(cap
 tools.cpp)

target_link_libraries(cap ${COMMON_LIBRARIES} input)

include_directories(${CMAKE_SOURCE_DIR}/debruijn)

if (DEFINED static_build)
  set_target_properties(cap PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()


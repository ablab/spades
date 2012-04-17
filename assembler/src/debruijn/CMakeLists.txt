############################################################################
# Copyright (c) 2011-2012 Saint-Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(spades CXX)

add_executable(spades
 main.cpp)

target_link_libraries(spades ${COMMON_LIBRARIES} input)

if (DEFINED static_build)
  set_target_properties(spades PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()


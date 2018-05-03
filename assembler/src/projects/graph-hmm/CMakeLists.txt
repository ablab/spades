############################################################################
# Copyright (c) 2018 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(graph-hmm CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(graph-hmm
               main.cpp)

target_link_libraries(graph-hmm input utils pipeline ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
  set_target_properties(graph-hmm PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

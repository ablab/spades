############################################################################
# Copyright (c) 2018 Saint Petersburg State University
# Copyright (c) 2019 University of Warwick
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(cds-subgraphs CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})


add_executable(cds-subgraphs
               stop_condon_finder.cpp
               subgraph_extraction.cpp
               main.cpp)

add_executable(cds-mapping-stats
               stats.cpp)

target_link_libraries(cds-subgraphs graphio common_modules ${COMMON_LIBRARIES})
target_link_libraries(cds-mapping-stats graphio common_modules ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
  set_target_properties(cds-subgraphs PROPERTIES LINK_SEARCH_END_STATIC 1)
  set_target_properties(cds-mapping-stats PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS cds-subgraphs
        DESTINATION bin
        COMPONENT runtime)

install(TARGETS cds-mapping-stats
        DESTINATION bin
        COMPONENT runtime)

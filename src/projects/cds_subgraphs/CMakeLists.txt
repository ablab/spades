############################################################################
# Copyright (c) 2023-2024 SPAdes team
# Copyright (c) 2018-2022 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(cds-subgraphs CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})


add_executable(cds-subgraphs
               stop_condon_finder.cpp
               subgraph_extraction.cpp
               cds_subgraphs.cpp)

add_executable(cds-mapping-stats
               stats.cpp)

add_executable(mag-improve
               mag_improve.cpp)

add_executable(unitig-coverage
               unitig_coverage.cpp profile_storage.cpp)
             
target_link_libraries(cds-subgraphs toolchain common_modules ${COMMON_LIBRARIES})
target_link_libraries(cds-mapping-stats toolchain common_modules ${COMMON_LIBRARIES})
target_link_libraries(mag-improve toolchain common_modules ${COMMON_LIBRARIES})
target_link_libraries(unitig-coverage toolchain graphio common_modules ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
  set_target_properties(cds-subgraphs PROPERTIES LINK_SEARCH_END_STATIC 1)
  set_target_properties(cds-mapping-stats PROPERTIES LINK_SEARCH_END_STATIC 1)
  set_target_properties(mag-improve PROPERTIES LINK_SEARCH_END_STATIC 1)
  set_target_properties(unitig-coverage PROPERTIES LINK_SEARCH_END_STATIC 1)  
endif()

install(TARGETS cds-subgraphs
        DESTINATION bin
        COMPONENT cds-subgraphs)

install(TARGETS cds-mapping-stats
        DESTINATION bin
        COMPONENT cds-subgraphs)

install(TARGETS mag-improve
        DESTINATION bin
        COMPONENT cds-subgraphs)

############################################################################
# Copyright (c) 2018 Saint Petersburg State University
# Copyright (c) 2019 University of Warwick
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(cds_subgraphs CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})


add_executable(cds_subgraphs
               main.cpp)

target_link_libraries(cds_subgraphs graphio common_modules ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
  set_target_properties(cds_subgraphs PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS cds_subgraphs
        DESTINATION bin
        COMPONENT runtime)

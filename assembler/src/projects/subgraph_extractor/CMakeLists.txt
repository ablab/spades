############################################################################
# Copyright (c) 2018 Saint Petersburg State University
# Copyright (c) 2019 University of Warwick
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(subgraph-extractor CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})


add_executable(subgraph-extractor
               main.cpp)

target_link_libraries(subgraph-extractor graphio common_modules ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
  set_target_properties(subgraph-extractor PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS subgraph-extractor
        DESTINATION bin
        COMPONENT runtime)

############################################################################
# Copyright (c) 2021 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(containment_index CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(cont-index
               main.cpp barcode_index_construction.cpp long_edge_statistics.cpp reference_path_checker.cpp
               scaffold_graph_helper.cpp multiplex_gfa_reader.cpp vertex_info_getter.cpp vertex_resolver.cpp)

target_link_libraries(cont-index spades-stages graphio toolchain common_modules ${COMMON_LIBRARIES})

if (SPADES_STATIC_BUILD)
    set_target_properties(cont-index PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS cont-index
        DESTINATION bin
        COMPONENT runtime)

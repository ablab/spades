############################################################################
# Copyright (c) 2021-2022 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(containment_index CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(splitter
               main.cpp barcode_index_construction.cpp cloud_only_links.cpp graph_resolver.cpp graph_resolver_io.cpp
               long_edge_statistics.cpp path_cluster_extractor.cpp path_extractor.cpp reference_path_checker.cpp
               scaffold_graph_helper.cpp vertex_info_getter.cpp)

target_link_libraries(splitter spades-stages graphio toolchain common_modules ${COMMON_LIBRARIES} auxiliary_graphs)

if (SPADES_STATIC_BUILD)
    set_target_properties(splitter PROPERTIES LINK_SEARCH_END_STATIC 1)
endif()

install(TARGETS splitter
        DESTINATION bin
        COMPONENT runtime)

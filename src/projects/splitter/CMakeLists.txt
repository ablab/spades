############################################################################
# Copyright (c) 2021-2023 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(splitter CXX)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_executable(splitter
               main.cpp barcode_index_construction.cpp graph_resolver.cpp graph_resolver_io.cpp
               path_extractor.cpp scaffold_graph_helper.cpp)

target_link_libraries(splitter graphio toolchain common_modules ${COMMON_LIBRARIES} auxiliary_graphs)

install(TARGETS splitter
        DESTINATION bin
        COMPONENT splitter)

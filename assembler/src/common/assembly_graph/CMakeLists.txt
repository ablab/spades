############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(assembly_graph CXX)

add_library(assembly_graph STATIC
            core/id_distributor.cpp
            components/connected_component.cpp
            paths/bidirectional_path.cpp
            paths/bidirectional_path_io/io_support.cpp
            paths/bidirectional_path_io/bidirectional_path_output.cpp
            graph_support/scaff_supplementary.cpp graph_support/coverage_uniformity_analyzer.cpp
            ../modules/alignment/edge_index_refiller.cpp)

target_link_libraries(assembly_graph utils llvm-support)
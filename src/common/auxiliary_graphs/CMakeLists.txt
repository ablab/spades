############################################################################
# Copyright (c) 2023-2025 SPAdes team
# Copyright (c) 2019-2022 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(auxiliary_graphs CXX)

add_library(auxiliary_graphs STATIC
            contracted_graph/contracted_graph.cpp
            contracted_graph/contracted_graph_builder.cpp
            contracted_graph/contracted_graph_helper.cpp
            contracted_graph/graph_condensation.cpp
            contracted_graph/contracted_statistics.cpp
            scaffold_graph/scaffold_vertex.cpp
)

target_link_libraries(auxiliary_graphs assembly_graph)

############################################################################
# Copyright (c) 2023-2024 SPAdes team
# Copyright (c) 2015-2022 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

project(path_extend CXX)

add_library(path_extend STATIC
            loop_detector.cpp
            loop_traverser.cpp
            gap_analyzer.cpp
            path_extenders.cpp
            pe_resolver.cpp
            overlap_remover.cpp
            pipeline/launch_support.cpp
            pipeline/launcher.cpp
            pipeline/extenders_logic.cpp
            scaffolder2015/extension_chooser2015.cpp
            scaffolder2015/scaffold_graph.cpp
            scaffolder2015/scaffold_graph_constructor.cpp
            scaffolder2015/scaffold_graph_visualizer.cpp
            scaffolder2015/connection_condition2015.cpp
            scaffolder2015/path_polisher.cpp
            scaffolder2015/scaffold_vertex.cpp
            read_cloud_path_extend/read_cloud_connection_conditions.cpp
            read_cloud_path_extend/scaffold_graph_construction_pipeline.cpp
            read_cloud_path_extend/scaffold_graph_storage.cpp
            read_cloud_path_extend/validation/scaffold_graph_validation.cpp
            read_cloud_path_extend/validation/transition_extractor.cpp
            read_cloud_path_extend/scaffold_graph_extractor.cpp
            read_cloud_path_extend/intermediate_scaffolding/scaffold_graph_gap_closer.cpp
            read_cloud_path_extend/transitions/transitions.cpp
            read_cloud_path_extend/intermediate_scaffolding/scaffold_vertex_predicates.cpp
            read_cloud_path_extend/path_cluster_helper.cpp
            read_cloud_path_extend/intermediate_scaffolding/pe_extraction.cpp
            read_cloud_path_extend/intermediate_scaffolding/predicate_builders.cpp
            read_cloud_path_extend/path_scaffolder.cpp
            read_cloud_path_extend/containment_index_threshold_finder.cpp
            read_cloud_path_extend/tip_connection_predicate.cpp
            read_cloud_path_extend/read_cloud_polisher_support.cpp)

target_link_libraries(path_extend assembly_graph ssw configs)



############################################################################
# Copyright (c) 2015 Saint Petersburg State University
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
            path_scaffolder.cpp
            pipeline/launch_support.cpp
            pipeline/launcher.cpp
            pipeline/extenders_logic.cpp
            scaffolder2015/extension_chooser2015.cpp
            scaffolder2015/scaffold_graph_constructor.cpp
            scaffolder2015/scaffold_graph_visualizer.cpp
            scaffolder2015/connection_condition2015.cpp
            scaffolder2015/path_polisher.cpp
            read_cloud_path_extend/scaffold_graph_extractor.cpp
            read_cloud_path_extend/extender_support/candidate_selectors.cpp
            read_cloud_path_extend/extender_support/entry_collectors.cpp
            read_cloud_path_extend/fragment_statistics/secondary_stats_estimators.cpp
            read_cloud_path_extend/intermediate_scaffolding/contracted_graph_from_simple.cpp
            read_cloud_path_extend/intermediate_scaffolding/pair_entry_processors.cpp
            read_cloud_path_extend/intermediate_scaffolding/path_cluster_helper.cpp
            read_cloud_path_extend/intermediate_scaffolding/predicate_builders.cpp
            read_cloud_path_extend/intermediate_scaffolding/scaffold_graph_path_cleaner.cpp
            read_cloud_path_extend/intermediate_scaffolding/scaffold_graph_polisher.cpp
            read_cloud_path_extend/scaffold_graph_construction/construction_callers.cpp
            read_cloud_path_extend/scaffold_graph_construction/containment_index_threshold_finder.cpp
            read_cloud_path_extend/scaffold_graph_construction/extender_searcher.cpp
            read_cloud_path_extend/scaffold_graph_construction/read_cloud_connection_conditions.cpp
            read_cloud_path_extend/scaffold_graph_construction/scaffold_graph_construction_pipeline.cpp
            read_cloud_path_extend/scaffold_graph_construction/scaffold_graph_storage.cpp
            read_cloud_path_extend/scaffold_graph_construction/scaffold_graph_storage_constructor.cpp
            read_cloud_path_extend/statistics/component_validation.cpp
            read_cloud_path_extend/statistics/long_edge_dataset.cpp
            read_cloud_path_extend/statistics/perfect_transitions.cpp
            read_cloud_path_extend/statistics/split_index_statistics.cpp
            read_cloud_path_extend/validation/path_cluster_validation.cpp
            read_cloud_path_extend/validation/pipeline_validation.cpp
            read_cloud_path_extend/validation/reference_path_index.cpp
            read_cloud_path_extend/validation/scaffold_graph_validation.cpp
            read_cloud_path_extend/validation/transition_extractor.cpp
            read_cloud_path_extend/validation/transition_subgraph_validation.cpp)

target_link_libraries(path_extend assembly_graph auxiliary_graphs ssw configs)



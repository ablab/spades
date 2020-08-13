//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once
#include "common/assembly_graph/paths/bidirectional_path_container.hpp"
#include "common/pipeline/graph_pack.hpp"


struct PathThreadingParams {
    bool use_distances = true;          //try to fit the distance between edges
    double good_distance_coeff = 0.10;   //a path is considered as good if its distance differs by at most good_distance_coeff
    double best_of_good_coeff = 0.05;   //if there are multiple *good& paths, the best one has to have distance difference < best_of_good_coeff * second_best_distance_difference
    bool extend_unique_paths = false;    //extend ends if possible
    bool make_transitive_connections = true; //check next max_steps_forward edges
    size_t max_steps_forward = 5;
    size_t max_distance = 50000;        //stop Dijkstra path search when path length exceeds max_distance or there are more than 3000 paths in buffer
    bool use_scaffolds = true;          //use pre-contructed scaffolds if path was not found
};

struct PathWithEdgePostions {

    PathWithEdgePostions(): edge_set(), positions() {}

    PathWithEdgePostions(const std::vector<std::vector<debruijn_graph::EdgeId>>&& paths_, const std::vector<int>&& positions_):
        edge_set(paths_), positions(positions_) {}

    std::vector<std::vector<debruijn_graph::EdgeId>> edge_set;
    std::vector<int> positions;
};

struct SeqString {
    std::string name;
    std::string seq;
};

using PathWithEdgePostionsContainer = std::vector<PathWithEdgePostions>;

PathWithEdgePostionsContainer ParseInputPaths(const std::string& transcript_paths_file, const debruijn_graph::Graph& graph);

path_extend::PathContainer Launch(debruijn_graph::GraphPack const & gp, 
                                  PathThreadingParams params,
                                  PathWithEdgePostionsContainer const & input_paths,
                                  std::vector<SeqString> & contigs,
                                  std::vector<std::string> const & paths_names,
                                  path_extend::PathContainer const & scaffolds);

std::vector<SeqString> ReadContigs(std::string const & contigs_file);
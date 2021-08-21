//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once
#include "common/assembly_graph/paths/bidirectional_path_container.hpp"
#include "common/pipeline/graph_pack.hpp"

namespace helpers {

// #define GOOD_NAME "contig_350"

struct PathThreadingParams {
    double score_domination_coeff = 0.50;   //if there are multiple *good* paths, the best one has to have score < score_domination_coeff * second_best_score
    double good_distance_coeff = 0.10;   //a path is considered as good if its distance differs by at most good_distance_coeff
    double best_of_good_coeff = 0.05;   //if there are multiple *good* paths, the best one has to have distance difference < best_of_good_coeff * second_best_distance_difference
    bool extend_unique_paths = false;    //extend ends if possible
    bool use_agressive_filling = true;  //use dijkstra's algorithm for edges connecting
    size_t max_steps_forward = 1;
    size_t max_distance = 50000;        //stop Dijkstra path search when path length exceeds max_distance or there are more than 3000 paths in buffer
    bool use_scaffolds = false;         //use pre-contructed scaffolds if path was not found
};

struct PathWithEdgePostions {

    PathWithEdgePostions() = default;

    PathWithEdgePostions(std::string path_name_, std::vector<debruijn_graph::EdgeId> paths_, std::vector<long long> start_positions_, std::vector<long long> end_positions_)
        : path_name(std::move(path_name_))
        , edges(std::move(paths_))
        , start_positions(std::move(start_positions_))
        , end_positions(std::move(end_positions_))
    {
        VERIFY(edges.size() == start_positions.size() && edges.size() == end_positions.size());
    }

    std::string path_name;
    std::vector<debruijn_graph::EdgeId> edges;
    std::vector<long long> start_positions; // inclusive
    std::vector<long long> end_positions;   // exclusive
};

struct SeqString {
    std::string name;
    std::string info;
    std::string seq;
    bool corrected = false;
};

using PathWithEdgePostionsContainer = std::vector<PathWithEdgePostions>;

path_extend::PathContainer Launch(debruijn_graph::GraphPack const & gp, 
                                  PathThreadingParams params,
                                  PathWithEdgePostionsContainer const & input_paths,
                                  std::vector<SeqString> & contigs,
                                  path_extend::PathContainer const & scaffolds,
                                  size_t nthreads);

std::vector<SeqString> ReadContigs(std::string const & contigs_file);

} // namespace helpers

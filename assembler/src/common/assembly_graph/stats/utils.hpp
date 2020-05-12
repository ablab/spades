//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/core/graph.hpp"
#include <vector>
#include <cstdlib>

namespace debruijn_graph {
namespace stats {

template<class Graph>
double AvgCoverage(const Graph &g, const std::vector<typename Graph::EdgeId> &edges) {
    double total_cov = 0.;
    size_t total_length = 0;
    for (auto e : edges) {
        total_cov += g.coverage(e) * (double) g.length(e);
        total_length += g.length(e);
    }
    return total_cov / (double) total_length;
}

template<class Graph>
size_t Nx(Graph &g, double percent) {
    size_t sum_edge_length = 0;
    std::vector<size_t> lengths;
    for (typename Graph::EdgeId e : g.edges()) {
        lengths.push_back(g.length(e));
        sum_edge_length += g.length(e);
    }
    std::sort(lengths.begin(), lengths.end());
    double len_perc = (1.0 - percent * 0.01) * (double) (sum_edge_length);
    for (size_t i = 0; i < lengths.size(); i++) {
        if (lengths[i] >= len_perc)
            return lengths[i];
        else
            len_perc -= (double) lengths[i];
    }
    return 0;
}

}
}

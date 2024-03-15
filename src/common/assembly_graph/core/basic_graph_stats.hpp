//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <algorithm>
#include <vector>
#include <cstddef>

namespace omnigraph {

template<class Graph>
class AvgCoverageCounter {
private:
    const Graph& graph_;
    const size_t min_length_;
public:
    AvgCoverageCounter(const Graph &graph, size_t min_length = 0)
            : graph_(graph), min_length_(min_length) {}

    double Count() const {
        double cov = 0;
        size_t length = 0;
        for (typename Graph::EdgeId e : graph_.edges()) {
            if (graph_.length(e) < min_length_)
                continue;

            cov += graph_.coverage(e) * (double) graph_.length(e);
            length += graph_.length(e);
        }
        if (length == 0)
            return 0.;
        return cov / (double) length;
    }
};

template<class Graph>
class CumulativeLengthCounter {
private:
    const Graph& graph_;
    const size_t min_length_;
public:
    CumulativeLengthCounter(const Graph &graph,
                            size_t min_length = 0)
            : graph_(graph), min_length_(min_length) {}

    size_t Count(bool canonical_only = true) const {
        size_t length = 0;
        if (canonical_only) {
            // FIXME: provide helper
            for (typename Graph::EdgeId e : graph_.canonical_edges()) {
                if (graph_.length(e) < min_length_)
                    continue;

                length += graph_.length(e);
            }
        } else {
            for (typename Graph::EdgeId e : graph_.edges()) {
                if (graph_.length(e) < min_length_)
                    continue;

                length += graph_.length(e);
            }
        }

        return length;
    }
};

template<class Graph>
size_t CumulativeLength(const Graph& g,
                        const std::vector<typename Graph::EdgeId>& path) {
    size_t s = 0;
    for (auto e : path)
        s += g.length(e);

    return s;
}

template<class Graph>
double AvgCoverage(const Graph& g,
                   const std::vector<typename Graph::EdgeId>& path) {
    double unnormalized_coverage = 0;
    size_t path_length = 0;
    for (auto edge : path) {
        size_t length = g.length(edge);
        path_length += length;
        unnormalized_coverage += g.coverage(edge) * (double) length;
    }
    return unnormalized_coverage / (double) path_length;
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
        if (double(lengths[i]) >= len_perc)
            return lengths[i];
        else
            len_perc -= (double) lengths[i];
    }
    return 0;
}

}

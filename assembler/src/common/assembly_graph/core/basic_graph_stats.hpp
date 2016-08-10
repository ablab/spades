#pragma once

#include "utils/standard_base.hpp"
namespace omnigraph {

template<class Graph>
class AvgCovereageCounter {
private:
    const Graph &graph_;
    const size_t min_length_;
public:
    AvgCovereageCounter(const Graph &graph, size_t min_length = 0) :
            graph_(graph), min_length_(min_length) {
    }

    double Count() const {
        double cov = 0;
        size_t length = 0;
        for (auto it = graph_.ConstEdgeBegin(); !it.IsEnd(); ++it) {
            if (graph_.length(*it) >= min_length_) {
                cov += graph_.coverage(*it) * (double) graph_.length(*it);
                length += graph_.length(*it);
            }
        }
        if (length == 0)
            return 0.;
        return cov / (double) length;
    }
};

template<class Graph>
size_t CumulativeLength(const Graph& g,
                        const std::vector<typename Graph::EdgeId>& path) {
    size_t s = 0;
    for (auto it = path.begin(); it != path.end(); ++it)
        s += g.length(*it);

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
}
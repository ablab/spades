#pragma once

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
}
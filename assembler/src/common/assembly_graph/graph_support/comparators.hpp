#pragma once

#include <functional>

namespace omnigraph {

template<class Graph>
struct CoverageComparator {
private:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    std::reference_wrapper<const Graph> graph_;
public:
    CoverageComparator(const Graph &graph)
            : graph_(graph) {}

    /**
     * Standard comparator function as used in collections.
     */
    double operator()(EdgeId edge) const {
        const Graph &g = graph_;
        return double(g.kmer_multiplicity(edge)) / double(g.length(edge));
    }

    bool operator()(EdgeId edge1, EdgeId edge2) const {
        const Graph &g = graph_;

        uint64_t lhs = g.kmer_multiplicity(edge1) * g.length(edge2),
                 rhs = g.kmer_multiplicity(edge2) * g.length(edge1);
        if (lhs < rhs)
            return true;
        else if (lhs == rhs)
            return edge1 < edge2;

        return false;
    }
};

template<class Graph>
struct LengthComparator {
  private:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    std::reference_wrapper<const Graph> graph_;
  public:
    LengthComparator(const Graph &graph)
            : graph_(graph) {}

    /**
     * Standard comparator function as used in collections.
     */
    size_t operator()(EdgeId edge) const {
        const Graph &g = graph_;
        return g.length(edge);
    }

    bool operator()(EdgeId edge1, EdgeId edge2) const {
        const Graph &g = graph_;

        size_t l1 = g.length(edge1), l2 = g.length(edge2);

        if (l1 < l2)
            return true;
        else if (l1 == l2)
            return edge1 < edge2;

        return false;
    }
};
}

#pragma once
namespace omnigraph {

template<class Graph>
struct CoverageComparator {
private:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    const Graph &graph_;
public:
    CoverageComparator(const Graph &graph)
            : graph_(graph) {
    }

    /**
     * Standard comparator function as used in collections.
     */
    bool operator()(EdgeId edge1, EdgeId edge2) const {
        if (math::eq(graph_.coverage(edge1), graph_.coverage(edge2))) {
            return edge1 < edge2;
        }
        return math::ls(graph_.coverage(edge1), graph_.coverage(edge2));
    }
};

/**
 * This class defines which edge is more likely to be tip. In this case we just assume shorter edges
 * are more likely tips then longer ones.
 */
    template<class Graph>
    struct LengthComparator {
    private:
        typedef typename Graph::EdgeId EdgeId;
        typedef typename Graph::VertexId VertexId;
        const Graph &graph_;
    public:
        /**
         * TipComparator should never be created with default constructor but it is necessary on order for
         * code to compile.
         */
        //  TipComparator() {
        //    VERIFY(false);
        //  }
        /**
         * Construct TipComparator for given graph
         * @param graph graph for which comparator is created
         */
        LengthComparator(const Graph &graph)
                : graph_(graph) {
        }

        /**
         * Standard comparator function as used in collections.
         */
        bool operator()(EdgeId edge1, EdgeId edge2) const {
            if (graph_.length(edge1) == graph_.length(edge2)) {
                return edge1 < edge2;
            }
            return graph_.length(edge1) < graph_.length(edge2);
        }
    };
}
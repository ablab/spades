#pragma once

#include "assembly_graph/graph_support/basic_vertex_conditions.hpp"
#include "assembly_graph/graph_support/graph_processing_algorithm.hpp"
#include "assembly_graph/graph_support/parallel_processing.hpp"

namespace omnigraph {

template<class Graph>
class Cleaner : public PersistentProcessingAlgorithm<Graph, typename Graph::VertexId> {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef PersistentProcessingAlgorithm<Graph, VertexId> base;
    typedef IsolatedVertexCondition<Graph> ConditionT;

    Graph &g_;
    ConditionT isolated_condition_;

public:
    Cleaner(Graph &g, size_t chunk_cnt = 1) :
            base(g,
                 std::make_shared<ParallelInterestingElementFinder<Graph, VertexId>>(ConditionT(g), chunk_cnt),
                    /*canonical only*/true),
            g_(g), isolated_condition_(g) {
    }

protected:

    bool Process(VertexId v) {
        if (isolated_condition_.Check(v)) {
            g_.DeleteVertex(v);
            return true;
        } else {
            return false;
        }
    }
};

template<class Graph>
size_t CleanIsolatedVertices(Graph &g, size_t chunk_cnt = 1) {
    Cleaner<Graph> cleaner(g, chunk_cnt);
    return cleaner.Run();
}

}

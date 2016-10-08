#pragma once
#include "assembly_graph/graph_support/parallel_processing.hpp"
#include "assembly_graph/graph_support/basic_vertex_conditions.hpp"
namespace omnigraph {

/**
* Compressor compresses vertices with unique incoming and unique outgoing edge in linear time while
* simple one-by-one compressing has square complexity.
*/
template<class Graph>
class Compressor {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef CompressCondition<Graph> ConditionT;

    Graph &graph_;
    ConditionT compress_condition_;
    bool safe_merging_;

    bool GoUniqueWayForward(EdgeId &e) {
        VertexId u = graph_.EdgeEnd(e);
        if (!graph_.CheckUniqueOutgoingEdge(u)
            || !graph_.CheckUniqueIncomingEdge(u)) {
            return false;
        }
        e = graph_.GetUniqueOutgoingEdge(u);
        return true;
    }

    bool GoUniqueWayBackward(EdgeId &e) {
        VertexId u = graph_.EdgeStart(e);
        if (!graph_.CheckUniqueOutgoingEdge(u)
            || !graph_.CheckUniqueIncomingEdge(u)) {
            return false;
        }
        e = graph_.GetUniqueIncomingEdge(u);
        return true;
    }

    //do not use without checks:)
    EdgeId CompressWithoutChecks(VertexId v) {
        EdgeId e = graph_.GetUniqueOutgoingEdge(v);
        EdgeId start_edge = e;
        while (GoUniqueWayBackward(e) && e != start_edge
               && !graph_.RelatedVertices(graph_.EdgeStart(e),
                                          graph_.EdgeEnd(e))) {
        }
        vector <EdgeId> mergeList;
        start_edge = e;
        do {
            mergeList.push_back(e);
        } while (GoUniqueWayForward(e) && e != start_edge
                 && !graph_.RelatedVertices(graph_.EdgeStart(e),
                                            graph_.EdgeEnd(e)));
        EdgeId new_edge = graph_.MergePath(mergeList, safe_merging_);
        TRACE("Vertex compressed and is now part of edge "
              << graph_.str(new_edge));
        return new_edge;

    }

public:
    Compressor(Graph& graph, bool safe_merging = true) :
            graph_(graph),
            compress_condition_(graph),
            safe_merging_(safe_merging) {
    }

    /**
     * Method compresses longest possible path, containing given vertex.
     * @param vertex to be compressed as part of a path
     * @return true if vertex can be compressed and false otherwise
     */
    bool CompressVertex(VertexId v) {
        return CompressVertexEdgeId(v) != EdgeId(0);
    }

    EdgeId CompressVertexEdgeId(VertexId v) {
        TRACE("Processing vertex " << graph_.str(v) << " started");
        if (!compress_condition_.Check(v)) {
            return EdgeId(0);
        }
        TRACE("Vertex " << graph_.str(v) << " judged compressible");
        return CompressWithoutChecks(v);
    }

//    bool IsOfInterest(VertexId v) const {
//        return CanCompressVertex(v);
//    }

private:
    DECL_LOGGER("Compressor")
};

template<class Graph>
class CompressingProcessor : public PersistentProcessingAlgorithm<Graph, typename Graph::VertexId> {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef PersistentProcessingAlgorithm<Graph, VertexId> base;
    typedef CompressCondition<Graph> ConditionT;

    Compressor<Graph> compressor_;
public:
    CompressingProcessor(Graph &graph, size_t chunk_cnt = 1, bool safe_merging = true) :
            base(graph,
                 std::make_shared<ParallelInterestingElementFinder<Graph, VertexId>>(ConditionT(graph), chunk_cnt),
                    /*canonical only*/true),
            compressor_(graph, safe_merging) {
    }

protected:
    bool Process(VertexId v) override {
        return compressor_.CompressVertex(v);
    }
};

/**
* Method compresses all vertices which can be compressed.
*/
template<class Graph>
bool CompressAllVertices(Graph &g, bool safe_merging = true, size_t chunk_cnt = 1) {
    CompressingProcessor<Graph> compressor(g, chunk_cnt, safe_merging);
    return compressor.Run();
}
}

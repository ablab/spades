#pragma once
#include "graph_support/parallel_processing.hpp"
#include "graph_support/basic_vertex_conditions.hpp"
namespace omnigraph {

/**
* Compressor compresses vertices with unique incoming and unique outgoing edge in linear time while
* simple one-by-one compressing has square complexity.
*/
template<class Graph>
class Compressor : public PersistentProcessingAlgorithm<Graph, typename Graph::VertexId,
        ParallelInterestingElementFinder < Graph, typename Graph::VertexId>> {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef PersistentProcessingAlgorithm <Graph,
    VertexId, ParallelInterestingElementFinder<Graph, VertexId>> base;
    typedef CompressCondition <Graph> ConditionT;

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
        //		e = graph_.conjugate(e);
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

//	//todo use graph method!
//	bool CanCompressVertex(VertexId v) const {
//		if (!graph_.CheckUniqueOutgoingEdge(v)
//			|| !graph_.CheckUniqueIncomingEdge(v)) {
//			TRACE(
//					"Vertex "
//							<< graph_.str(v)
//							<< " judged NOT compressible. Proceeding to the next vertex");
//			TRACE("Processing vertex " << graph_.str(v) << " finished");
//			return false;
//		}
//		return true;
//	}
public:
    Compressor(Graph &graph, size_t chunk_cnt = 1, bool safe_merging = true) :
            base(graph,
                 ParallelInterestingElementFinder<Graph, VertexId>(graph,
                                                                   ConditionT(graph), chunk_cnt),
                    /*canonical only*/true),
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
        TRACE("Processing vertex " << graph_.str(v) << " started");
        if (!compress_condition_.Check(v)) {
            return false;
        }
        TRACE("Vertex " << graph_.str(v) << " judged compressible");
        CompressWithoutChecks(v);
        return true;
    }

    EdgeId CompressVertexEdgeId(VertexId v) {
        TRACE("Processing vertex " << graph_.str(v) << " started");
        if (!compress_condition_.Check(v)) {
            return EdgeId(0);
        }
        TRACE("Vertex " << graph_.str(v) << " judged compressible");
        return CompressWithoutChecks(v);
    }

//	bool IsOfInterest(VertexId v) const {
//	    return CanCompressVertex(v);
//	}

protected:
    bool Process(VertexId v) override {
        if (compress_condition_.Check(v)) {
            CompressWithoutChecks(v);
            return true;
        } else {
            return false;
        }
    }

private:
    DECL_LOGGER("Compressor")
};

/**
* Method compresses all vertices which can be compressed.
*/
template<class Graph>
bool CompressAllVertices(Graph &g, bool safe_merging = true, size_t chunk_cnt = 1) {
    Compressor<Graph> compressor(g, chunk_cnt, safe_merging);
    return compressor.Run();
}
}
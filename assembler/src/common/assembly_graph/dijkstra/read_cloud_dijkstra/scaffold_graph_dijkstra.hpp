#pragma once
#include "common/assembly_graph/dijkstra/dijkstra_helper.hpp"
#include "modules/path_extend/scaffolder2015/scaffold_graph.hpp"
#include "common/barcode_index/barcode_info_extractor.hpp"
#include "read_cloud_path_extend/scaffold_graph_construction/read_cloud_connection_conditions.hpp"

namespace omnigraph {
using path_extend::scaffold_graph::ScaffoldGraph;

template <class Graph, typename distance_t = size_t>
class SimpleScaffoldGraphLengthCalculator {
 protected:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
 public:
    distance_t GetLength(EdgeId ) const {
        return 1;
    }
};

template <class Graph, typename distance_t = size_t>
class DistanceBasedScaffoldGraphLengthCalculator {
 protected:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    const Graph &graph_;
 public:
    explicit DistanceBasedScaffoldGraphLengthCalculator(const Graph& graph) : graph_(graph) {}
    distance_t GetLength(EdgeId edge) const {
        return graph_.length(edge) + graph_.length(edge.getEnd());
    }
};

template <>
class ForwardNeighbourIterator<ScaffoldGraph> {
    typedef typename ScaffoldGraph::VertexId VertexId;
    typedef typename ScaffoldGraph::EdgeId EdgeId;
    typedef typename vector<ScaffoldGraph::EdgeId>::const_iterator edge_const_iterator;
    vector<EdgeId> out_edges_;
    edge_const_iterator current_;
 public:
    ForwardNeighbourIterator(const ScaffoldGraph &graph, VertexId vertex) :
        out_edges_(graph.OutgoingEdges(vertex)), current_(out_edges_.begin()) { }

    bool HasNext() {
        return current_ != out_edges_.end();
    }

    vertex_neighbour<ScaffoldGraph> Next() {
        TRACE("Before increment");
        TRACE(current_->getStart().int_id() << ", " << current_->getEnd().int_id());
        vertex_neighbour<ScaffoldGraph> res(current_->getEnd(), *current_);
        current_++;
//        TRACE("After increment");
//        TRACE(current_->getStart().int_id() << ", " << current_->getEnd().int_id());
        return res;
    }

    DECL_LOGGER("ScaffoldForwardNeighbourItetator");
};

template<>
class BackwardNeighbourIterator<ScaffoldGraph> {
    typedef typename ScaffoldGraph::VertexId VertexId;
    typedef typename ScaffoldGraph::EdgeId EdgeId;
    typedef typename vector<EdgeId>::const_iterator edge_const_iterator;

    vector<EdgeId> in_edges_;
    edge_const_iterator current_;
 public:
    BackwardNeighbourIterator(const ScaffoldGraph &graph, VertexId vertex) :
        in_edges_(graph.IncomingEdges(vertex)), current_(in_edges_.begin()) { }

    bool HasNext() {
        return current_ != in_edges_.end();
    }

    vertex_neighbour<ScaffoldGraph> Next() {
        vertex_neighbour<ScaffoldGraph> res(current_->getStart(), *current_);
        current_++;
        return res;
    }
};

template<class Graph, typename distance_t = size_t>
class ScaffoldBarcodedPathPutChecker {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    const Graph& g_;
    const VertexId first_;
    const VertexId second_;
    shared_ptr<path_extend::read_cloud::ScaffoldVertexPredicate> predicate_;

 public:
    ScaffoldBarcodedPathPutChecker(const Graph& g, const VertexId& first, const VertexId& second,
                                   shared_ptr<path_extend::read_cloud::ScaffoldVertexPredicate> predicate) :
        g_(g),
        first_(first),
        second_(second),
        predicate_(predicate) {
        TRACE("Construction");
        TRACE("First id: " << first_.int_id());
        TRACE("Second id: " << second_.int_id());
    }

    bool Check(VertexId vertex, EdgeId /*unused*/, distance_t distance) const {
        TRACE("Checking vertex " << g_.str(vertex));
        TRACE("Id: " << vertex.int_id());
        TRACE("First id: " << first_.int_id());
        TRACE("Second id: " << second_.int_id());
        bool target_reached = distance > 0 and (vertex == first_ or vertex == second_);
        if (target_reached) {
            TRACE("Target reached");
            return false;
        }
        TRACE("Checking");
        return predicate_->Check(vertex);
    }
    DECL_LOGGER("ScaffoldBarcodePutChecker");
};

template<class Graph, typename distance_t = size_t>
class StartPredicateProcessChecker {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    const Graph& g_;
    const VertexId start_;
    const func::TypedPredicate<VertexId>& predicate_;
 public:
    StartPredicateProcessChecker(const Graph &g_,
                                 const VertexId start_,
                                 const func::TypedPredicate<VertexId>& predicate_)
        : g_(g_), start_(start_), predicate_(predicate_) {}

    bool Check(VertexId vertex, distance_t /*distance*/) {
        return vertex == start_ or not predicate_(vertex);
    }
};

template<class Graph, typename distance_t = size_t>
class TrivialScaffoldPutChecker {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

 public:
    TrivialScaffoldPutChecker() {}

    bool Check(VertexId /*unused*/, EdgeId /*unused*/, distance_t /*unused*/) const {
        return true;
    }
    DECL_LOGGER("TrivialScaffoldPutChecker");
};

typedef ComposedDijkstraSettings<path_extend::scaffold_graph::ScaffoldGraph,
                                 DistanceBasedScaffoldGraphLengthCalculator<path_extend::scaffold_graph::ScaffoldGraph>,
                                 StartPredicateProcessChecker<path_extend::scaffold_graph::ScaffoldGraph>,
                                 TrivialScaffoldPutChecker<path_extend::scaffold_graph::ScaffoldGraph>,
                                 ForwardNeighbourIteratorFactory<path_extend::scaffold_graph::ScaffoldGraph> >
    PredicateBasedScaffoldDijkstraSettings;

typedef Dijkstra<path_extend::scaffold_graph::ScaffoldGraph, PredicateBasedScaffoldDijkstraSettings> PredicateBasedScaffoldDijkstra;

//forward scaffold dijkstra

typedef ComposedDijkstraSettings<path_extend::scaffold_graph::ScaffoldGraph,
                                 SimpleScaffoldGraphLengthCalculator<path_extend::scaffold_graph::ScaffoldGraph>,
                                 BoundedVertexTargetedProcessChecker<path_extend::scaffold_graph::ScaffoldGraph>,
                                 ScaffoldBarcodedPathPutChecker<path_extend::scaffold_graph::ScaffoldGraph>,
                                 ForwardNeighbourIteratorFactory<path_extend::scaffold_graph::ScaffoldGraph> >
    ForwardBoundedScaffoldDijkstraSettings;

typedef Dijkstra<path_extend::scaffold_graph::ScaffoldGraph, ForwardBoundedScaffoldDijkstraSettings> ForwardBoundedScaffoldDijkstra;

//backward scaffold dijkstra

typedef ComposedDijkstraSettings<path_extend::scaffold_graph::ScaffoldGraph,
                                 SimpleScaffoldGraphLengthCalculator<path_extend::scaffold_graph::ScaffoldGraph>,
                                 BoundedVertexTargetedProcessChecker<path_extend::scaffold_graph::ScaffoldGraph>,
                                 ScaffoldBarcodedPathPutChecker<path_extend::scaffold_graph::ScaffoldGraph>,
                                 BackwardNeighbourIteratorFactory<path_extend::scaffold_graph::ScaffoldGraph> >
    BackwardBoundedScaffoldDijkstraSettings;

typedef Dijkstra<path_extend::scaffold_graph::ScaffoldGraph, BackwardBoundedScaffoldDijkstraSettings> BackwardBoundedScaffoldDijkstra;


class ScaffoldDijkstraHelper {
 public:
    static BackwardBoundedScaffoldDijkstra CreateBackwardBoundedScaffoldDijkstra(
        const path_extend::scaffold_graph::ScaffoldGraph&graph,
        const ScaffoldGraph::ScaffoldGraphVertex first,
        const ScaffoldGraph::ScaffoldGraphVertex second,
        size_t length_bound,
        shared_ptr<path_extend::read_cloud::ScaffoldVertexPredicate> predicate,
        size_t max_vertex_number = -1ul){
        return BackwardBoundedScaffoldDijkstra(graph, BackwardBoundedScaffoldDijkstraSettings(
            SimpleScaffoldGraphLengthCalculator<path_extend::scaffold_graph::ScaffoldGraph>(),
            BoundedVertexTargetedProcessChecker<path_extend::scaffold_graph::ScaffoldGraph>(first, length_bound),
            ScaffoldBarcodedPathPutChecker<path_extend::scaffold_graph::ScaffoldGraph>(graph, first, second, predicate),
            BackwardNeighbourIteratorFactory<path_extend::scaffold_graph::ScaffoldGraph>(graph)),
                                               max_vertex_number);
    }

    static ForwardBoundedScaffoldDijkstra CreateForwardBoundedScaffoldDijkstra(
        const path_extend::scaffold_graph::ScaffoldGraph& graph,
        const ScaffoldGraph::ScaffoldGraphVertex& first,
        const ScaffoldGraph::ScaffoldGraphVertex& second,
        size_t length_bound,
        shared_ptr<path_extend::read_cloud::ScaffoldVertexPredicate> predicate,
        size_t max_vertex_number = -1ul){
        return ForwardBoundedScaffoldDijkstra(graph, ForwardBoundedScaffoldDijkstraSettings(
            SimpleScaffoldGraphLengthCalculator<path_extend::scaffold_graph::ScaffoldGraph>(),
            BoundedVertexTargetedProcessChecker<path_extend::scaffold_graph::ScaffoldGraph>(second, length_bound),
            ScaffoldBarcodedPathPutChecker<path_extend::scaffold_graph::ScaffoldGraph>(graph, first, second, predicate),
            ForwardNeighbourIteratorFactory<path_extend::scaffold_graph::ScaffoldGraph>(graph)),
                                              max_vertex_number);
    }

    static PredicateBasedScaffoldDijkstra CreatePredicateBasedScaffoldDijkstra(
        const path_extend::scaffold_graph::ScaffoldGraph &graph,
        const path_extend::scaffold_graph::ScaffoldGraph::ScaffoldGraphVertex &vertex,
        const func::TypedPredicate<path_extend::scaffold_graph::ScaffoldVertex>& predicate,
        size_t max_vertex_number = -1ul){
        return PredicateBasedScaffoldDijkstra(graph, PredicateBasedScaffoldDijkstraSettings(
            DistanceBasedScaffoldGraphLengthCalculator<path_extend::scaffold_graph::ScaffoldGraph>(graph),
            StartPredicateProcessChecker<path_extend::scaffold_graph::ScaffoldGraph>(graph, vertex, predicate),
            TrivialScaffoldPutChecker<path_extend::scaffold_graph::ScaffoldGraph>(),
            ForwardNeighbourIteratorFactory<path_extend::scaffold_graph::ScaffoldGraph>(graph)),
                                              max_vertex_number);
    }
};
}
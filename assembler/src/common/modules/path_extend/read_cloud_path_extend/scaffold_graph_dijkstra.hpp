#pragma once
#include "common/assembly_graph/dijkstra/dijkstra_helper.hpp"
#include "modules/path_extend/scaffolder2015/scaffold_graph.hpp"
#include "common/barcode_index/barcode_info_extractor.hpp"
#include "read_cloud_connection_conditions.hpp"

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
class DistanceBasedScaffoldGraphLengthCalculator: public LengthCalculator<Graph, distance_t> {
 protected:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
 public:
    explicit DistanceBasedScaffoldGraphLengthCalculator(const Graph& graph) : LengthCalculator<Graph, distance_t>(graph) {}
    distance_t GetLength(EdgeId edge) const override {
        return edge.getLength();
    }
};

template <>
class ForwardNeighbourIterator<ScaffoldGraph> : public NeighbourIterator<ScaffoldGraph>{
    typedef typename ScaffoldGraph::VertexId VertexId;
    typedef typename ScaffoldGraph::EdgeId EdgeId;
    typedef typename vector<ScaffoldGraph::EdgeId>::const_iterator edge_const_iterator;
    vector<EdgeId> out_edges_;
    edge_const_iterator current_;
 public:
    ForwardNeighbourIterator(const ScaffoldGraph &graph, VertexId vertex) :
        NeighbourIterator<ScaffoldGraph>(graph, vertex),
        out_edges_(graph.OutgoingEdges(vertex)), current_(out_edges_.begin()) { }

    bool HasNext() override {
        return current_ != out_edges_.end();
    }

    vertex_neighbour<ScaffoldGraph> Next() override {
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
class BackwardNeighbourIterator<ScaffoldGraph> : public NeighbourIterator<ScaffoldGraph>{
    typedef typename ScaffoldGraph::VertexId VertexId;
    typedef typename ScaffoldGraph::EdgeId EdgeId;
    typedef typename vector<EdgeId>::const_iterator edge_const_iterator;

    vector<EdgeId> in_edges_;
    edge_const_iterator current_;
 public:
    BackwardNeighbourIterator(const ScaffoldGraph &graph, VertexId vertex) :
        NeighbourIterator<ScaffoldGraph>(graph, vertex),
        in_edges_(graph.IncomingEdges(vertex)), current_(in_edges_.begin()) { }

    bool HasNext() override {
        return current_ != in_edges_.end();
    }

    vertex_neighbour<ScaffoldGraph> Next() override {
        vertex_neighbour<ScaffoldGraph> res(current_->getStart(), *current_);
        current_++;
        return res;
    }
};

template<class Graph, typename distance_t = size_t>
class ScaffoldBarcodedPathPutChecker : public VertexPutChecker<Graph, distance_t> {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    const Graph& g_;
    const VertexId first_;
    const VertexId second_;
    shared_ptr<path_extend::ScaffoldVertexPredicate> predicate_;

 public:
    ScaffoldBarcodedPathPutChecker(const Graph& g, const VertexId& first, const VertexId& second,
                                   shared_ptr<path_extend::ScaffoldVertexPredicate> predicate) :
        VertexPutChecker<Graph, distance_t>(),
        g_(g),
        first_(first),
        second_(second),
        predicate_(predicate) {
        TRACE("Construction");
        TRACE("First id: " << first_.int_id());
        TRACE("Second id: " << second_.int_id());
    }

    bool Check(VertexId vertex, EdgeId /*unused*/, distance_t distance) const override {
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
class LengthBasedProcessChecker : public VertexProcessChecker<Graph, distance_t> {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    const Graph& g_;
    const VertexId start_;
    const size_t length_bound_;
 public:
    LengthBasedProcessChecker(const Graph &g_, const VertexId start_, const size_t length_bound_)
        : g_(g_), start_(start_), length_bound_(length_bound_) {}

    bool Check(VertexId vertex, distance_t /*distance*/) override {
        return vertex == start_ or g_.length(vertex) <= length_bound_;
    }
};

template<class Graph, typename distance_t = size_t>
class TrivialScaffoldPutChecker : public VertexPutChecker<Graph, distance_t> {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

 public:
    TrivialScaffoldPutChecker() {}

    bool Check(VertexId /*unused*/, EdgeId /*unused*/, distance_t /*unused*/) const override {
        return true;
    }
    DECL_LOGGER("TrivialScaffoldPutChecker");
};

typedef ComposedDijkstraSettings<path_extend::scaffold_graph::ScaffoldGraph,
                                 SimpleScaffoldGraphLengthCalculator<path_extend::scaffold_graph::ScaffoldGraph>,
                                 LengthBasedProcessChecker<path_extend::scaffold_graph::ScaffoldGraph>,
                                 TrivialScaffoldPutChecker<path_extend::scaffold_graph::ScaffoldGraph>,
                                 ForwardNeighbourIteratorFactory<path_extend::scaffold_graph::ScaffoldGraph> >
    LengthBasedScaffoldDijkstraSettings;

typedef Dijkstra<path_extend::scaffold_graph::ScaffoldGraph, LengthBasedScaffoldDijkstraSettings> LengthBasedScaffoldDijkstra;

static LengthBasedScaffoldDijkstra CreateLengthBasedScaffoldDijkstra(
    const path_extend::scaffold_graph::ScaffoldGraph& graph,
    const path_extend::scaffold_graph::ScaffoldGraph::ScaffoldGraphVertex& vertex,
    size_t length_bound,
    size_t max_vertex_number = -1ul){
    return LengthBasedScaffoldDijkstra(graph, LengthBasedScaffoldDijkstraSettings(
        SimpleScaffoldGraphLengthCalculator<path_extend::scaffold_graph::ScaffoldGraph>(),
        LengthBasedProcessChecker<path_extend::scaffold_graph::ScaffoldGraph>(graph, vertex, length_bound),
        TrivialScaffoldPutChecker<path_extend::scaffold_graph::ScaffoldGraph>(),
        ForwardNeighbourIteratorFactory<path_extend::scaffold_graph::ScaffoldGraph>(graph)),
                                          max_vertex_number);
}

//forward scaffold dijkstra

typedef ComposedDijkstraSettings<path_extend::scaffold_graph::ScaffoldGraph,
                                 SimpleScaffoldGraphLengthCalculator<path_extend::scaffold_graph::ScaffoldGraph>,
                                 BoundProcessChecker<path_extend::scaffold_graph::ScaffoldGraph>,
                                 ScaffoldBarcodedPathPutChecker<path_extend::scaffold_graph::ScaffoldGraph>,
                                 ForwardNeighbourIteratorFactory<path_extend::scaffold_graph::ScaffoldGraph> >
    ForwardBoundedScaffoldDijkstraSettings;

typedef Dijkstra<path_extend::scaffold_graph::ScaffoldGraph, ForwardBoundedScaffoldDijkstraSettings> ForwardBoundedScaffoldDijkstra;

static ForwardBoundedScaffoldDijkstra CreateForwardBoundedScaffoldDijkstra(
        const path_extend::scaffold_graph::ScaffoldGraph& graph,
        const ScaffoldGraph::ScaffoldGraphVertex& first,
        const ScaffoldGraph::ScaffoldGraphVertex& second,
        size_t length_bound,
        shared_ptr<path_extend::ScaffoldVertexPredicate> predicate,
        size_t max_vertex_number = -1ul){
    return ForwardBoundedScaffoldDijkstra(graph, ForwardBoundedScaffoldDijkstraSettings(
        SimpleScaffoldGraphLengthCalculator<path_extend::scaffold_graph::ScaffoldGraph>(),
        BoundProcessChecker<path_extend::scaffold_graph::ScaffoldGraph>(length_bound),
        ScaffoldBarcodedPathPutChecker<path_extend::scaffold_graph::ScaffoldGraph>(graph, first, second, predicate),
        ForwardNeighbourIteratorFactory<path_extend::scaffold_graph::ScaffoldGraph>(graph)),
                           max_vertex_number);
}

//backward scaffold dijkstra

typedef ComposedDijkstraSettings<path_extend::scaffold_graph::ScaffoldGraph,
                                 SimpleScaffoldGraphLengthCalculator<path_extend::scaffold_graph::ScaffoldGraph>,
                                 BoundProcessChecker<path_extend::scaffold_graph::ScaffoldGraph>,
                                 ScaffoldBarcodedPathPutChecker<path_extend::scaffold_graph::ScaffoldGraph>,
                                 BackwardNeighbourIteratorFactory<path_extend::scaffold_graph::ScaffoldGraph> >
    BackwardBoundedScaffoldDijkstraSettings;

typedef Dijkstra<path_extend::scaffold_graph::ScaffoldGraph, BackwardBoundedScaffoldDijkstraSettings> BackwardBoundedScaffoldDijkstra;

static BackwardBoundedScaffoldDijkstra CreateBackwardBoundedScaffoldDijkstra(
        const path_extend::scaffold_graph::ScaffoldGraph&graph,
        const ScaffoldGraph::ScaffoldGraphVertex first,
        const ScaffoldGraph::ScaffoldGraphVertex second,
        size_t length_bound,
        shared_ptr<path_extend::ScaffoldVertexPredicate> predicate,
        size_t max_vertex_number = -1ul){
    return BackwardBoundedScaffoldDijkstra(graph, BackwardBoundedScaffoldDijkstraSettings(
        SimpleScaffoldGraphLengthCalculator<path_extend::scaffold_graph::ScaffoldGraph>(),
        BoundProcessChecker<path_extend::scaffold_graph::ScaffoldGraph>(length_bound),
        ScaffoldBarcodedPathPutChecker<path_extend::scaffold_graph::ScaffoldGraph>(graph, first, second, predicate),
        BackwardNeighbourIteratorFactory<path_extend::scaffold_graph::ScaffoldGraph>(graph)),
                                          max_vertex_number);
}
}
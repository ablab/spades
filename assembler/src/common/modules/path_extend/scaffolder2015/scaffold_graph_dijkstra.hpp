//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/dijkstra/dijkstra_helper.hpp"
#include "auxiliary_graphs/scaffold_graph/scaffold_graph.hpp"
#include "modules/path_extend/scaffolder2015/scaffold_vertex_predicates.hpp"

namespace omnigraph {
template<>
class ForwardNeighbourIterator<scaffold_graph::ScaffoldGraph> {
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef typename ScaffoldGraph::VertexId VertexId;
    typedef typename ScaffoldGraph::EdgeId EdgeId;
    typedef typename std::vector<ScaffoldGraph::EdgeId>::const_iterator edge_const_iterator;
    std::vector<EdgeId> out_edges_;
    edge_const_iterator current_;
  public:
    ForwardNeighbourIterator(const ScaffoldGraph &graph, VertexId vertex) :
        out_edges_(graph.OutgoingEdges(vertex)), current_(out_edges_.begin()) {}

    bool HasNext() {
        return current_ != out_edges_.end();
    }

    vertex_neighbour<ScaffoldGraph> Next() {
        TRACE("Before increment");
        TRACE(current_->getStart().int_id() << ", " << current_->getEnd().int_id());
        vertex_neighbour<ScaffoldGraph> res(current_->getEnd(), *current_);
        current_++;
        return res;
    }

    DECL_LOGGER("ScaffoldForwardNeighbourItetator");
};

template<>
class BackwardNeighbourIterator<scaffold_graph::ScaffoldGraph> {
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef typename ScaffoldGraph::VertexId VertexId;
    typedef typename ScaffoldGraph::EdgeId EdgeId;
    typedef typename std::vector<EdgeId>::const_iterator edge_const_iterator;

    std::vector<EdgeId> in_edges_;
    edge_const_iterator current_;
  public:
    BackwardNeighbourIterator(const ScaffoldGraph &graph, VertexId vertex) :
        in_edges_(graph.IncomingEdges(vertex)), current_(in_edges_.begin()) {}

    bool HasNext() {
        return current_ != in_edges_.end();
    }

    vertex_neighbour<ScaffoldGraph> Next() {
        vertex_neighbour<ScaffoldGraph> res(current_->getStart(), *current_);
        current_++;
        return res;
    }
};
}

namespace path_extend {

namespace scaffolder {

template<class Graph, typename distance_t = size_t>
class SimpleScaffoldGraphLengthCalculator {
  protected:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
  public:
    distance_t GetLength(EdgeId) const {
        return 1;
    }
};

template<class Graph, typename distance_t = size_t>
class DistanceBasedScaffoldGraphLengthCalculator {
  protected:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    const Graph &graph_;
  public:
    explicit DistanceBasedScaffoldGraphLengthCalculator(const Graph &graph) : graph_(graph) {}
    distance_t GetLength(EdgeId edge) const {
        return graph_.length(edge) + graph_.length(edge.getEnd());
    }
};

template<class Graph, typename distance_t = size_t>
class ScaffoldBarcodedPathPutChecker {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    const Graph &g_;
    const VertexId first_;
    const VertexId second_;
    std::shared_ptr<scaffolder::ScaffoldVertexPredicate> predicate_;

  public:
    ScaffoldBarcodedPathPutChecker(const Graph &g, const VertexId &first, const VertexId &second,
                                   std::shared_ptr<scaffolder::ScaffoldVertexPredicate> predicate) :
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

    const Graph &g_;
    const VertexId start_;
    const func::TypedPredicate<VertexId> &predicate_;
  public:
    StartPredicateProcessChecker(const Graph &g,
                                 const VertexId &start,
                                 const func::TypedPredicate<VertexId> &predicate)
        : g_(g), start_(start), predicate_(predicate) {}

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

typedef ComposedDijkstraSettings<scaffold_graph::ScaffoldGraph,
                                 DistanceBasedScaffoldGraphLengthCalculator<scaffold_graph::ScaffoldGraph>,
                                 StartPredicateProcessChecker<scaffold_graph::ScaffoldGraph>,
                                 TrivialScaffoldPutChecker<scaffold_graph::ScaffoldGraph>,
                                 ForwardNeighbourIteratorFactory<scaffold_graph::ScaffoldGraph> >
    PredicateBasedScaffoldDijkstraSettings;

typedef Dijkstra<scaffold_graph::ScaffoldGraph, PredicateBasedScaffoldDijkstraSettings>
    PredicateBasedScaffoldDijkstra;

//forward scaffold dijkstra

typedef ComposedDijkstraSettings<scaffold_graph::ScaffoldGraph,
                                 SimpleScaffoldGraphLengthCalculator<scaffold_graph::ScaffoldGraph>,
                                 BoundedVertexTargetedProcessChecker<scaffold_graph::ScaffoldGraph>,
                                 ScaffoldBarcodedPathPutChecker<scaffold_graph::ScaffoldGraph>,
                                 ForwardNeighbourIteratorFactory<scaffold_graph::ScaffoldGraph> >
    ForwardBoundedScaffoldDijkstraSettings;

typedef Dijkstra<scaffold_graph::ScaffoldGraph, ForwardBoundedScaffoldDijkstraSettings>
    ForwardBoundedScaffoldDijkstra;

//backward scaffold dijkstra

typedef ComposedDijkstraSettings<scaffold_graph::ScaffoldGraph,
                                 SimpleScaffoldGraphLengthCalculator<scaffold_graph::ScaffoldGraph>,
                                 BoundedVertexTargetedProcessChecker<scaffold_graph::ScaffoldGraph>,
                                 ScaffoldBarcodedPathPutChecker<scaffold_graph::ScaffoldGraph>,
                                 BackwardNeighbourIteratorFactory<scaffold_graph::ScaffoldGraph> >
    BackwardBoundedScaffoldDijkstraSettings;

typedef Dijkstra<scaffold_graph::ScaffoldGraph, BackwardBoundedScaffoldDijkstraSettings>
    BackwardBoundedScaffoldDijkstra;

class ScaffoldDijkstraHelper {
  public:
    static BackwardBoundedScaffoldDijkstra CreateBackwardBoundedScaffoldDijkstra(
        const scaffold_graph::ScaffoldGraph &graph,
        const scaffold_graph::ScaffoldVertex first,
        const scaffold_graph::ScaffoldVertex second,
        size_t length_bound,
        std::shared_ptr<scaffolder::ScaffoldVertexPredicate> predicate,
        size_t max_vertex_number = -1ul) {
        return BackwardBoundedScaffoldDijkstra(graph, BackwardBoundedScaffoldDijkstraSettings(
            SimpleScaffoldGraphLengthCalculator<scaffold_graph::ScaffoldGraph>(),
            BoundedVertexTargetedProcessChecker<scaffold_graph::ScaffoldGraph>(first, length_bound),
            ScaffoldBarcodedPathPutChecker<scaffold_graph::ScaffoldGraph>(graph, first, second, predicate),
            BackwardNeighbourIteratorFactory<scaffold_graph::ScaffoldGraph>(graph)),
                                               max_vertex_number);
    }

    static ForwardBoundedScaffoldDijkstra CreateForwardBoundedScaffoldDijkstra(
        const scaffold_graph::ScaffoldGraph &graph,
        const scaffold_graph::ScaffoldVertex &first,
        const scaffold_graph::ScaffoldVertex &second,
        size_t length_bound,
        std::shared_ptr<scaffolder::ScaffoldVertexPredicate> predicate,
        size_t max_vertex_number = -1ul) {
        return ForwardBoundedScaffoldDijkstra(graph, ForwardBoundedScaffoldDijkstraSettings(
            SimpleScaffoldGraphLengthCalculator<scaffold_graph::ScaffoldGraph>(),
            BoundedVertexTargetedProcessChecker<scaffold_graph::ScaffoldGraph>(second, length_bound),
            ScaffoldBarcodedPathPutChecker<scaffold_graph::ScaffoldGraph>(graph, first, second, predicate),
            ForwardNeighbourIteratorFactory<scaffold_graph::ScaffoldGraph>(graph)),
                                              max_vertex_number);
    }

    static PredicateBasedScaffoldDijkstra CreatePredicateBasedScaffoldDijkstra(
        const scaffold_graph::ScaffoldGraph &graph,
        const scaffold_graph::ScaffoldVertex &vertex,
        const func::TypedPredicate<scaffold_graph::ScaffoldVertex> &predicate,
        size_t max_vertex_number = -1ul) {
        return PredicateBasedScaffoldDijkstra(graph, PredicateBasedScaffoldDijkstraSettings(
            DistanceBasedScaffoldGraphLengthCalculator<scaffold_graph::ScaffoldGraph>(graph),
            StartPredicateProcessChecker<scaffold_graph::ScaffoldGraph>(graph, vertex, predicate),
            TrivialScaffoldPutChecker<scaffold_graph::ScaffoldGraph>(),
            ForwardNeighbourIteratorFactory<scaffold_graph::ScaffoldGraph>(graph)),
                                              max_vertex_number);
    }
};
}
}
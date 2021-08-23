//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once
#include "dijkstra_algorithm.hpp"

namespace omnigraph {

template<class Graph>
class DijkstraHelper {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
public:
    typedef Dijkstra<Graph,
                     ComposedDijkstraSettings<Graph,
                                              LengthCalculator<Graph>,
                                              VertexProcessChecker<Graph>,
                                              VertexPutChecker<Graph>,
                                              UnorientedNeighbourIteratorFactory<Graph> > > UnorientedDijkstra;

    //------------------------------

    typedef Dijkstra<Graph,
                     ComposedDijkstraSettings<Graph,
                                              LengthCalculator<Graph>,
                                              VertexProcessChecker<Graph>,
                                              VertexPutChecker<Graph>,
                                              BackwardNeighbourIteratorFactory<Graph> > > BackwardDijkstra;

    //------------------------------
    // bounded dijkstra
    //------------------------------
    typedef ComposedDijkstraSettings<Graph,
                                     LengthCalculator<Graph>,
                                     BoundProcessChecker<Graph>,
                                     BoundPutChecker<Graph>,
                                     ForwardNeighbourIteratorFactory<Graph> > BoundedDijkstraSettings;

    template<bool EnableTraceback = false>
    using BoundedDijkstraGeneric = Dijkstra<Graph, BoundedDijkstraSettings, EnableTraceback>;
    using BoundedDijkstra = BoundedDijkstraGeneric<false>;
    using BoundedDijkstraWithTraceback = BoundedDijkstraGeneric<true>;
    
    static auto CreateBoundedDijkstra(const Graph &graph, size_t length_bound,
                                                   size_t max_vertex_number = -1ul) {
        return BoundedDijkstra(graph,
                               BoundedDijkstraSettings(
                                   LengthCalculator<Graph>(graph),
                                   BoundProcessChecker<Graph>(length_bound),
                                   BoundPutChecker<Graph>(length_bound),
                                   ForwardNeighbourIteratorFactory<Graph>(graph)),
                               max_vertex_number);
    }

    static auto CreateBoundedDijkstraWithTraceback(const Graph &graph, size_t length_bound,
                                                   size_t max_vertex_number = -1ul) {
        return BoundedDijkstraWithTraceback(graph,
                                            BoundedDijkstraSettings(
                                                LengthCalculator<Graph>(graph),
                                                BoundProcessChecker<Graph>(length_bound),
                                                BoundPutChecker<Graph>(length_bound),
                                                ForwardNeighbourIteratorFactory<Graph>(graph)),
                                            max_vertex_number);
    }

    typedef ComposedDijkstraSettings<Graph,
            LengthCalculator<Graph>,
            BoundProcessChecker<Graph>,
            VertexPutChecker<Graph>,
            ForwardNeighbourIteratorFactory<Graph> > EdgeBoundedDijkstraSettings;

    typedef Dijkstra<Graph, EdgeBoundedDijkstraSettings> EdgeBoundedDijkstra;

    static EdgeBoundedDijkstra
    CreateEdgeBoundedDijkstra(const Graph &graph, size_t length_bound,
                              size_t max_vertex_number = -1ul){
        return EdgeBoundedDijkstra(graph,
                                   EdgeBoundedDijkstraSettings(
                                       LengthCalculator<Graph>(graph),
                                       BoundProcessChecker<Graph>(length_bound),
                                       VertexPutChecker<Graph>(),
                                       ForwardNeighbourIteratorFactory<Graph>(graph)),
                                   max_vertex_number);
    }

    typedef ComposedDijkstraSettings<Graph,
                                     LengthCalculator<Graph>,
                                     BoundProcessChecker<Graph>,
                                     CompositePutChecker<Graph>,
                                     ForwardNeighbourIteratorFactory<Graph> > LengthBoundedDijkstraSettings;

    typedef Dijkstra<Graph, LengthBoundedDijkstraSettings> LengthBoundedDijkstra;

    static LengthBoundedDijkstra CreateLengthBoundedDijkstra(const Graph& graph, size_t length_threshold,
                                                             size_t distance_bound, size_t max_vertex_number = -1ul) {
        auto length_put_checker = std::make_shared<LengthPutChecker<Graph>>(graph, length_threshold);
        auto bound_put_checker = std::make_shared<BoundPutChecker<Graph>>(distance_bound);
        std::vector<std::shared_ptr<VertexPutChecker<Graph>>> put_checkers;
        put_checkers.push_back(length_put_checker);
        put_checkers.push_back(bound_put_checker);
        return LengthBoundedDijkstra(graph, LengthBoundedDijkstraSettings(
            LengthCalculator<Graph>(graph),
            BoundProcessChecker<Graph>(distance_bound),
            CompositePutChecker<Graph>(put_checkers),
            ForwardNeighbourIteratorFactory<Graph>(graph)),
                                     max_vertex_number);
    }

    //------------------------------
    // bounded backward dijkstra
    //------------------------------

    typedef ComposedDijkstraSettings<Graph,
                                     LengthCalculator<Graph>,
                                     BoundProcessChecker<Graph>,
                                     BoundPutChecker<Graph>,
                                     BackwardNeighbourIteratorFactory<Graph> > BackwardBoundedDijkstraSettings;

    typedef Dijkstra<Graph, BackwardBoundedDijkstraSettings> BackwardBoundedDijkstra;

    static BackwardBoundedDijkstra
    CreateBackwardBoundedDijkstra(const Graph &graph,
                                  size_t bound,
                                  size_t max_vertex_number = size_t(-1)) {
        return BackwardBoundedDijkstra(graph,
                                       BackwardBoundedDijkstraSettings(
                                           LengthCalculator<Graph>(graph),
                                           BoundProcessChecker<Graph>(bound),
                                           BoundPutChecker<Graph>(bound),
                                           BackwardNeighbourIteratorFactory<Graph>(graph)),
                                       max_vertex_number);
    }

    //------------------------------
    // Undirected bounded dijkstra
    //------------------------------

    typedef ComposedDijkstraSettings<Graph,
            LengthCalculator<Graph>,
            BoundProcessChecker<Graph>,
            BoundPutChecker<Graph>,
            UnorientedNeighbourIteratorFactory<Graph> > UnorientedBoundedDijkstraSettings;

    typedef Dijkstra<Graph, UnorientedBoundedDijkstraSettings> UnorientedBoundedDijkstra;

    static UnorientedBoundedDijkstra
    CreateUnorientedBoundedDijkstra(const Graph &graph,
                                  size_t bound,
                                  size_t max_vertex_number = size_t(-1)) {
        return UnorientedBoundedDijkstra(graph,
                                       UnorientedBoundedDijkstraSettings(
                                               LengthCalculator<Graph>(graph),
                                               BoundProcessChecker<Graph>(bound),
                                               BoundPutChecker<Graph>(bound),
                                               UnorientedNeighbourIteratorFactory<Graph>(graph)),
                                       max_vertex_number);
    }


    typedef ComposedDijkstraSettings<Graph,
                                     LengthCalculator<Graph>,
                                     BoundProcessChecker<Graph>,
                                     VertexPutChecker<Graph>,
                                     BackwardNeighbourIteratorFactory<Graph> > BackwardEdgeBoundedDijkstraSettings;

    typedef Dijkstra<Graph, BackwardEdgeBoundedDijkstraSettings> BackwardEdgeBoundedDijkstra;

    static BackwardEdgeBoundedDijkstra
    CreateBackwardEdgeBoundedDijkstra(const Graph &graph,
                                      size_t bound, size_t max_vertex_number = size_t(-1)){
        return BackwardEdgeBoundedDijkstra(graph,
                                           BackwardEdgeBoundedDijkstraSettings(
                                               LengthCalculator<Graph>(graph),
                                               BoundProcessChecker<Graph>(bound),
                                               VertexPutChecker<Graph>(),
                                               BackwardNeighbourIteratorFactory<Graph>(graph)),
                                           max_vertex_number);
    }

    //------------------------------

    typedef Dijkstra<Graph,
                     ComposedDijkstraSettings<Graph,
                                              LengthCalculator<Graph>,
                                              VertexProcessChecker<Graph>,
                                              EdgeComponentPutChecker<Graph>,
                                              UnorientedNeighbourIteratorFactory<Graph> > > ComponentFinder;
    //------------------------------

    typedef Dijkstra<Graph,
                     ComposedDijkstraSettings<Graph,
                                              ComponentLenCalculator<Graph>,
                                              BoundProcessChecker<Graph>,
                                              VertexPutChecker<Graph>,
                                              UnorientedNeighbourIteratorFactory<Graph> > > NeighbourhoodFinder;
    //------------------------------

    typedef Dijkstra<Graph,
                     ComposedDijkstraSettings<Graph,
                                              LengthCalculator<Graph>,
                                              VertexProcessChecker<Graph>,
                                              SubgraphPutChecker<Graph>,
                                              UnorientedNeighbourIteratorFactory<Graph> > > SubgraphDijkstra;

    typedef ComposedDijkstraSettings<Graph,
                                     PathIgnoringLengthCalculator<Graph>,
                                     BoundProcessChecker<Graph>,
                                     BoundPutChecker<Graph>,
                                     ForwardNeighbourIteratorFactory<Graph> > PathIgnoringDijkstraSettings;


    //------------------------------
    // short edge dijkstra settings
    //------------------------------
    typedef ComposedDijkstraSettings<Graph,
                                     BoundedEdgeLenCalculator<Graph>,
                                     ZeroLengthProcessChecker<Graph>,
                                     VertexPutChecker<Graph>,
                                     UnorientedNeighbourIteratorFactory<Graph> > ShortEdgeDijkstraSettings;

    typedef Dijkstra<Graph, ShortEdgeDijkstraSettings> ShortEdgeDijkstra;

    static ShortEdgeDijkstra CreateShortEdgeDijkstra(const Graph &graph, size_t edge_length_bound,
                                                     size_t max_vertex_number = size_t(-1)) {
        return ShortEdgeDijkstra(graph,
                                 ShortEdgeDijkstraSettings(BoundedEdgeLenCalculator<Graph>(graph, edge_length_bound),
                                                           ZeroLengthProcessChecker<Graph>(),
                                                           VertexPutChecker<Graph>(),
                                                           UnorientedNeighbourIteratorFactory<Graph>(graph)),
                                 max_vertex_number);
    }

    //------------------------------
    // counting dijkstra
    //------------------------------
    typedef CountingDijkstraSettings<Graph,
                                     UnorientedNeighbourIteratorFactory<Graph> > UnorientCountingDijkstraSettings;

    typedef Dijkstra<Graph, UnorientCountingDijkstraSettings> CountingDijkstra;

    static CountingDijkstra CreateCountingDijkstra(const Graph &graph, size_t max_size,
                                                   size_t edge_length_bound,
                                                   size_t max_vertex_number = size_t(-1)) {
        return CountingDijkstra(graph,
                                UnorientCountingDijkstraSettings(graph,
                                                                 UnorientedNeighbourIteratorFactory<Graph>(graph),
                                                                 max_size, edge_length_bound),
                                max_vertex_number);
    }


    //------------------------------
    // targeted bounded dijkstra
    //------------------------------

    typedef ComposedDijkstraSettings<Graph,
                                     LengthCalculator<Graph>,
                                     BoundedVertexTargetedProcessChecker<Graph>,
                                     BoundPutChecker<Graph>,
                                     ForwardNeighbourIteratorFactory<Graph> > TargetedBoundedDijkstraSettings;

    typedef Dijkstra<Graph, TargetedBoundedDijkstraSettings> TargetedBoundedDijkstra;

    static TargetedBoundedDijkstra CreateTargetedBoundedDijkstra(const Graph &graph,
                                                                 VertexId target_vertex, size_t bound,
                                                                 size_t max_vertex_number = size_t(-1)) {
        return TargetedBoundedDijkstra(graph,
                                       TargetedBoundedDijkstraSettings(LengthCalculator<Graph>(graph),
                                                                       BoundedVertexTargetedProcessChecker<Graph>(target_vertex, bound),
                                                                       BoundPutChecker<Graph>(bound),
                                                                       ForwardNeighbourIteratorFactory<Graph>(graph)),
                                       max_vertex_number);
    }
    //------------------------------
    // coverage bounded dijkstra
    //------------------------------
    typedef ComposedDijkstraSettings<Graph,
                                     LengthCalculator<Graph>,
                                     BoundProcessChecker<Graph>,
                                     CoveragePutChecker<Graph>,
                                     ForwardNeighbourIteratorFactory<Graph> > CoverageBoundedDijkstraSettings;

    typedef Dijkstra<Graph, CoverageBoundedDijkstraSettings> CoverageBoundedDijkstra;

    static CoverageBoundedDijkstra CreateCoverageBoundedDijkstra(const Graph &graph, size_t length_bound, double min_coverage,
                                                                 size_t max_vertex_number = -1ul) {
        return CoverageBoundedDijkstra(graph,
                                       CoverageBoundedDijkstraSettings(LengthCalculator<Graph>(graph),
                                                                       BoundProcessChecker<Graph>(length_bound),
                                                                       CoveragePutChecker<Graph>(min_coverage, graph, length_bound),
                                                                       ForwardNeighbourIteratorFactory<Graph>(graph)),
                                       max_vertex_number);
    }

    //------------------------------
    // Length, distance, and component dijkstra
    //------------------------------
    typedef std::tuple<ForbiddenEdgesPutChecker<Graph>, LengthPutChecker<Graph>, BoundPutChecker<Graph>> DilationTuple;

    typedef ComposedDijkstraSettings<Graph,
                                     LengthCalculator<Graph>,
                                     BoundProcessChecker<Graph>,
                                     AndPutChecker<DilationTuple, Graph>,
                                     ForwardNeighbourIteratorFactory<Graph> > DilationDijkstraSettings;

    typedef Dijkstra<Graph, DilationDijkstraSettings> DilationDijkstra;

    static DilationDijkstra CreateDilationDijkstra(const Graph &graph,
                                                   const std::unordered_set<EdgeId> &forbidden_edges,
                                                   size_t length_threshold,
                                                   size_t distance_bound,
                                                   size_t max_vertex_number = -1ul) {
        ForbiddenEdgesPutChecker<Graph> forbidden_put_checker(forbidden_edges);
        LengthPutChecker<Graph> length_put_checker(graph, length_threshold);
        BoundPutChecker<Graph> bound_put_checker(distance_bound);
        auto put_checkers = std::make_tuple(forbidden_put_checker, length_put_checker, bound_put_checker);
        return DilationDijkstra(graph, DilationDijkstraSettings(
            LengthCalculator<Graph>(graph),
            BoundProcessChecker<Graph>(distance_bound),
            AndPutChecker<DilationTuple, Graph>(put_checkers),
            ForwardNeighbourIteratorFactory<Graph>(graph)),
                                     max_vertex_number);
    }
};

template<class Graph>
typename DijkstraHelper<Graph>::BoundedDijkstra CreateBoundedDijkstra(const Graph &graph, size_t length_bound,
                                                                      size_t max_vertex_number = -1ul) {
    return DijkstraHelper<Graph>::CreateBoundedDijkstra(graph, length_bound, max_vertex_number);
}

template<class Graph>
typename DijkstraHelper<Graph>::BackwardBoundedDijkstra CreateBackwardBoundedDijkstra(const Graph &graph, size_t length_bound,
                                                                                      size_t max_vertex_number = -1ul) {
    return DijkstraHelper<Graph>::CreateBackwardBoundedDijkstra(graph, length_bound, max_vertex_number);
}

template<class Graph>
typename DijkstraHelper<Graph>::EdgeBoundedDijkstra CreateEdgeBoundedDijkstra(const Graph &graph, size_t length_bound,
                                                                              size_t max_vertex_number = -1ul) {
    return DijkstraHelper<Graph>::CreateEdgeBoundedDijkstra(graph, length_bound, max_vertex_number);
}

template<class Graph>
typename DijkstraHelper<Graph>::BackwardEdgeBoundedDijkstra CreateBackwardEdgeBoundedDijkstra(const Graph &graph, size_t length_bound,
                                                                                              size_t max_vertex_number = -1ul) {
    return DijkstraHelper<Graph>::CreateBackwardEdgeBoundedDijkstra(graph, length_bound, max_vertex_number);
}

}

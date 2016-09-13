//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
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
    typedef Dijkstra<Graph, ComposedDijkstraSettings<Graph,
            LengthCalculator<Graph>,
                VertexProcessChecker<Graph>,
                VertexPutChecker<Graph>,
                UnorientedNeighbourIteratorFactory<Graph> > > UnorientedDijkstra;

    //------------------------------

    typedef Dijkstra<Graph, ComposedDijkstraSettings<Graph,
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

    typedef Dijkstra<Graph, BoundedDijkstraSettings> BoundedDijkstra;

    static BoundedDijkstra CreateBoundedDijkstra(const Graph &graph, size_t length_bound,
            size_t max_vertex_number = -1ul){
        return BoundedDijkstra(graph, BoundedDijkstraSettings(
                        LengthCalculator<Graph>(graph),
                        BoundProcessChecker<Graph>(length_bound),
                        BoundPutChecker<Graph>(length_bound),
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

    static BackwardBoundedDijkstra CreateBackwardBoundedDijkstra(const Graph &graph,
            size_t bound, size_t max_vertex_number = size_t(-1)){
        return BackwardBoundedDijkstra(graph, BackwardBoundedDijkstraSettings(
                LengthCalculator<Graph>(graph),
                BoundProcessChecker<Graph>(bound),
                BoundPutChecker<Graph>(bound),
                BackwardNeighbourIteratorFactory<Graph>(graph)), max_vertex_number);
    }

    //------------------------------

    typedef Dijkstra<Graph, ComposedDijkstraSettings<Graph,
            LengthCalculator<Graph>,
            VertexProcessChecker<Graph>,
            EdgeComponentPutChecker<Graph>,
            UnorientedNeighbourIteratorFactory<Graph> > > ComponentFinder;
    //------------------------------

    typedef Dijkstra<Graph, ComposedDijkstraSettings<Graph,
            ComponentLenCalculator<Graph>,
            BoundProcessChecker<Graph>,
            VertexPutChecker<Graph>,
            UnorientedNeighbourIteratorFactory<Graph> > > NeighbourhoodFinder;
    //------------------------------

    typedef Dijkstra<Graph, ComposedDijkstraSettings<Graph,
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
            size_t max_vertex_number = size_t(-1)){
        return ShortEdgeDijkstra(graph, ShortEdgeDijkstraSettings(
                        BoundedEdgeLenCalculator<Graph>(graph, edge_length_bound),
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
            size_t edge_length_bound, size_t max_vertex_number = size_t(-1)){
        return CountingDijkstra(graph, UnorientCountingDijkstraSettings(graph,
                        UnorientedNeighbourIteratorFactory<Graph>(graph),
                        max_size, edge_length_bound), max_vertex_number);
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
            VertexId target_vertex, size_t bound, size_t max_vertex_number = size_t(-1)){
        return TargetedBoundedDijkstra(graph,
                TargetedBoundedDijkstraSettings(LengthCalculator<Graph>(graph),
                        BoundedVertexTargetedProcessChecker<Graph>(target_vertex, bound),
                        BoundPutChecker<Graph>(bound),
                        ForwardNeighbourIteratorFactory<Graph>(graph)),
                max_vertex_number);
    }
};

}

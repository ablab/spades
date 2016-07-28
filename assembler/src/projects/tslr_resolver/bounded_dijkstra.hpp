#pragma once
#include "barcode_mapper.hpp"
#include "../../modules/algorithms/dijkstra/vertex_put_checker.hpp"

namespace omnigraph {

    template<class Graph, typename distance_t = size_t>
    class LengthPutChecker : public VertexPutChecker<Graph, distance_t> {
        typedef typename Graph::VertexId VertexId;
        typedef typename Graph::EdgeId EdgeId;
        const Graph& g_;
        const distance_t bound_;
    public:
        LengthPutChecker(distance_t bound, const Graph& g) : VertexPutChecker<Graph, distance_t>(),
                                                             g_(g), bound_(bound) { }
        bool Check(VertexId, EdgeId edge, distance_t) const {
            return g_.length(edge) < bound_;
        }
    };

    template <class Graph>
    class LengthDijkstra {
        using LengthBoundedDijkstraSettings = ComposedDijkstraSettings<Graph,
                LengthCalculator<Graph>,
                BoundProcessChecker<Graph>,
                LengthPutChecker<Graph>,
                ForwardNeighbourIteratorFactory<Graph> >;

        using LengthBoundedDijkstra = Dijkstra<Graph, LengthBoundedDijkstraSettings>;

    public:
        static LengthBoundedDijkstra CreateLengthBoundedDijkstra(const Graph &graph, size_t distance_bound,
                                    size_t length_bound, size_t max_vertex_number = -1ul){
            return LengthBoundedDijkstra(graph, LengthBoundedDijkstraSettings(
                    LengthCalculator<Graph>(graph),
                    BoundProcessChecker<Graph>(distance_bound),
                    LengthPutChecker<Graph>(length_bound, graph),
                    ForwardNeighbourIteratorFactory<Graph>(graph)),
                                         max_vertex_number);
        }
    };
}

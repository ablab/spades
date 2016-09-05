#pragma once
#include "barcode_mapper.hpp"
#include "tslr_extension_chooser.hpp"
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

    template<class Graph, typename distance_t = size_t>
    class BarcodePutChecker : public VertexPutChecker<Graph, distance_t> {
        typedef typename Graph::VertexId VertexId;
        typedef typename Graph::EdgeId EdgeId;
        typedef shared_ptr<tslr_resolver::BarcodeMapper> Bmapper;
        const Graph& g_;
        const distance_t length_bound_;
        const double barcode_threshold_;
        Bmapper mapper_;
        EdgeId decisive_edge_;
        ScaffoldingUniqueEdgeStorage unique_storage_;
        vector <EdgeId>& candidates_;
    public:
        BarcodePutChecker(const Graph& g, 
            const distance_t& length_bound, 
            const double& barcode_threshold, 
            const Bmapper& mapper, 
            const EdgeId& decisive_edge,
            const ScaffoldingUniqueEdgeStorage& unique_storage,
            vector<EdgeId>& candidates) : VertexPutChecker<Graph, distance_t> (),
                                                             g_(g), 
                                                             length_bound_(length_bound), 
                                                             barcode_threshold_(barcode_threshold),
                                                             mapper_(mapper), 
                                                             decisive_edge_(decisive_edge),
                                                             unique_storage_(unique_storage),
                                                             candidates_(candidates) { }
        bool Check(VertexId, EdgeId edge, distance_t) const {
            DEBUG("Checking edge " << edge.int_id());
            DEBUG("Length " << g_.length(edge)) 
            DEBUG("decisive_edge " << decisive_edge_.int_id())
            DEBUG("intersection " << mapper_->IntersectionSize(decisive_edge_, edge))
            DEBUG("Barcodes " << mapper_->GetSizeHeads(edge))
            DEBUG("Normalized intersection " << mapper_->IntersectionSizeNormalizedByFirst(decisive_edge_, edge))
            DEBUG("Is unique " << unique_storage_.IsUnique(edge))
            if (g_.length(edge) < length_bound_) {
                DEBUG("Short edge, passed" << endl)  //todo use short edges to reduce number of candidates
                return true;
            }
            if (! unique_storage_.IsUnique(edge)) {
                DEBUG("Long non-unique nearby edge, passed" << endl)
                return true;
            }
            else {
                candidates_.push_back(edge);
                DEBUG("Long unique nearby edge, put to candidates list and stop" << endl)
                return false;
            }
        }

        DECL_LOGGER("BarcodePutChecker")
    };

    template <class Graph>
    class BarcodeDijkstra {
        using BarcodeBoundedDijkstraSettings = ComposedDijkstraSettings<Graph,
                LengthCalculator<Graph>,
                BoundProcessChecker<Graph>,
                BarcodePutChecker<Graph>,
                ForwardNeighbourIteratorFactory<Graph> >;

        using BarcodeBoundedDijkstra = Dijkstra<Graph, BarcodeBoundedDijkstraSettings>;

    public:
        static BarcodeBoundedDijkstra CreateBarcodeBoundedDijkstra(const Graph &graph, 
            size_t distance_bound, 
            BarcodePutChecker<Graph>& put_checker, 
            size_t max_vertex_number = -1ul){
            return BarcodeBoundedDijkstra(graph, BarcodeBoundedDijkstraSettings(
                    LengthCalculator<Graph>(graph),
                    BoundProcessChecker<Graph>(distance_bound),
                    put_checker,
                    ForwardNeighbourIteratorFactory<Graph>(graph)),
                                         max_vertex_number);
        }
    }; 
}

#pragma once

#include "common/modules/path_extend/scaffolder2015/connection_condition2015.hpp"
#include "common/assembly_graph/graph_support/scaff_supplementary.hpp"

namespace omnigraph {
//------------------------------
// unique based dijkstra with paired connections
//------------------------------
template<class Graph>
class PairedConnectionIterator : public NeighbourIterator<Graph> {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef path_extend::PairedLibConnectionCondition paired_condition_t;
    const paired_condition_t paired_connection_condition_;
    map<EdgeId, double> connection_map_;
    typename map<EdgeId, double>::const_iterator iterator_;

 public:
    PairedConnectionIterator(const Graph& graph, VertexId vertex, EdgeId edge,
                             paired_condition_t paired_connection_condition) :
        NeighbourIterator<Graph>(graph, vertex),
        paired_connection_condition_(paired_connection_condition),
        connection_map_(paired_connection_condition.ConnectedWith(edge)),
        iterator_(connection_map_.begin()) {}

    bool HasNext() {
        return iterator_ != connection_map_.end();
    }

    vertex_neighbour<Graph> Next() {
        vertex_neighbour<Graph> result(this->graph_.EdgeEnd(iterator_->first), iterator_->first);
        iterator_++;
        return result;
    }

};

template<class Graph>
class NeighbourIteratorPointerWrapper : public NeighbourIterator<Graph> {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    shared_ptr<NeighbourIterator<Graph>> iterator_ptr_;

 public:
    NeighbourIteratorPointerWrapper(const Graph& graph,
                                    VertexId vertex,
                                    const shared_ptr<NeighbourIterator<Graph>> actual_iterator) :
        NeighbourIterator<Graph>(graph, vertex),
        iterator_ptr_(actual_iterator) {}

    bool HasNext() {
        return iterator_ptr_->HasNext();
    }

    vertex_neighbour<Graph> Next() {
        return iterator_ptr_->Next();
    }
};

template<class Graph>
class PairedConnectionNeighbourIteratorFactory {
    typedef path_extend::PairedLibConnectionCondition paired_condition_t;
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    const Graph& graph_;
    paired_condition_t paired_connection_condition_;
 public:
    typedef NeighbourIteratorPointerWrapper<Graph> NeighbourIterator;
    PairedConnectionNeighbourIteratorFactory(const Graph& graph, const paired_condition_t& paired_connection_condition)
        :
        graph_(graph), paired_connection_condition_(paired_connection_condition) {}

    NeighbourIterator CreateIterator(VertexId vertex) {
        if (IsEndOfTip(vertex)) {
            auto iterator_ptr = make_shared<ForwardNeighbourIterator<Graph>>(graph_, vertex);
            return NeighbourIteratorPointerWrapper<Graph>(graph_, vertex, iterator_ptr);
        } else {
            const EdgeId tip = *graph_.IncomingEdges(vertex).begin();
            auto iterator_ptr =
                make_shared<PairedConnectionIterator<Graph>>(graph_, vertex, tip, paired_connection_condition_);
            return NeighbourIteratorPointerWrapper<Graph>(graph_, vertex, iterator_ptr);
        }
    }

 private:
    bool IsEndOfTip(const VertexId& vertex) const {
        return (graph_.OutgoingEdgeCount(vertex) == 0 and graph_.IncomingEdgeCount(vertex) == 1);
    }
};

//    typedef ComposedDijkstraSettings<Graph,
//            LengthCalculator<Graph>,
//            BoundProcessChecker<Graph>,
//            UniquePutChecker<Graph>,
//            PairedConnectionNeighbourIteratorFactory<Graph> > PairedDijkstraSettings;
//
//    typedef Dijkstra<Graph, PairedDijkstraSettings> PairedDijkstra;
//
//    static PairedDijkstra CreatePairedDijkstra(const Graph &graph, size_t length_bound,
//                                               const path_extend::ScaffoldingUniqueEdgeStorage unique_storage,
//                                               const path_extend::PairedLibConnectionCondition& paired_connection_condition,
//                                               size_t max_vertex_number = -1ul) {
//        return PairedDijkstra(graph,
//                PairedDijkstraSettings(
//                    LengthCalculator<Graph>(graph),
//                    BoundProcessChecker<Graph>(length_bound),
//                    UniquePutChecker<Graph>(graph, unique_storage),
//                    PairedConnectionNeighbourIteratorFactory<Graph>(graph, paired_connection_condition)),
//                max_vertex_number);
//    }

//------------------------------
// unique based dijkstra
//------------------------------

template<class Graph, typename distance_t = size_t>
class UniquePutChecker : public VertexPutChecker<Graph, distance_t> {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    const Graph& g_;
    const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage_;

 public:
    UniquePutChecker(const Graph& g,
                     const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage) :
        VertexPutChecker<Graph, distance_t>(),
        g_(g),
        unique_storage_(unique_storage) {}

    //This method performs simple checks on nearby edges to compose candidate set for further filtering
    bool Check(VertexId, EdgeId edge, distance_t dist) const {
        DEBUG("Checking edge " << edge.int_id());
        DEBUG("Length " << g_.length(edge))
        size_t gap = dist - g_.length(edge);
        DEBUG("Gap " << gap)
        DEBUG("Is unique " << unique_storage_.IsUnique(edge));

        if (g_.length(edge) < unique_storage_.min_length()) {
            DEBUG("Short edge, passed" << endl)  //todo use short edges to reduce number of candidates
            return true;
        }
        if (!unique_storage_.IsUnique(edge)) {
            DEBUG("Long non-unique nearby edge, passed" << endl)
            return true;
        } else {
            DEBUG("Long unique nearby edge, put to candidates list and stop" << endl)
            return false;
        }
    }

    DECL_LOGGER("UniquePutChecker")
};

typedef ComposedDijkstraSettings<Graph,
                                 LengthCalculator<Graph>,
                                 BoundProcessChecker<Graph>,
                                 UniquePutChecker<Graph>,
                                 ForwardNeighbourIteratorFactory<Graph> > UniqueDijkstraSettings;

typedef Dijkstra<Graph, UniqueDijkstraSettings> UniqueDijkstra;

static UniqueDijkstra CreateUniqueDijkstra(const Graph& graph, size_t length_bound,
                                           const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage,
                                           size_t max_vertex_number = -1ul) {
    return UniqueDijkstra(graph,
                          UniqueDijkstraSettings(
                              LengthCalculator<Graph>(graph),
                              BoundProcessChecker<Graph>(length_bound),
                              UniquePutChecker<Graph>(graph, unique_storage),
                              ForwardNeighbourIteratorFactory<Graph>(graph)),
                          max_vertex_number);
}

//------------------------------
// barcode based dijkstra
//------------------------------

template<class Graph, typename distance_t = size_t>
class BarcodedPathPutChecker : public VertexPutChecker<Graph, distance_t> {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    const Graph& g_;
    const barcode_index::FrameBarcodeIndexInfoExtractor& barcode_extractor_;
    const vector<barcode_index::BarcodeId>& barcode_intersection_;
    const size_t barcode_threshold_;
    const size_t count_threshold_;
    const size_t tail_threshold_;
    const size_t edge_length_threshold_;

 public:
    BarcodedPathPutChecker(const Graph& g, const barcode_index::FrameBarcodeIndexInfoExtractor& barcode_info_extractor,
                           const vector<barcode_index::BarcodeId>& barcodes,
                           size_t barcode_threshold, size_t count_threshold,
                           size_t tail_threshold, size_t edge_length_threshold) :
        VertexPutChecker<Graph, distance_t>(),
        g_(g),
        barcode_extractor_(barcode_info_extractor),
        barcode_intersection_(barcodes),
        barcode_threshold_(barcode_threshold),
        count_threshold_(count_threshold),
        tail_threshold_(tail_threshold),
        edge_length_threshold_(edge_length_threshold) {}

    bool Check(VertexId, EdgeId edge, distance_t ) const override {
        DEBUG("Checking edge " << edge.int_id());
        DEBUG("Length " << g_.length(edge))
        if (g_.length(edge) <= edge_length_threshold_) {
            return true;
        }
        vector<barcode_index::BarcodeId> middle_intersection;
        vector<barcode_index::BarcodeId> middle_barcodes = barcode_extractor_.GetBarcodes(edge);
        std::set_intersection(barcode_intersection_.begin(), barcode_intersection_.end(),
                              middle_barcodes.begin(), middle_barcodes.end(),
                              std::back_inserter(middle_intersection));
        return middle_intersection.size() >= barcode_threshold_;
    }

 private:
    bool CheckBarcode(const EdgeId& first, const EdgeId& second, barcode_index::BarcodeId barcode,
                      size_t count_threshold, size_t tail_threshold) const {
        return barcode_extractor_.GetMaxPos(first, barcode) + tail_threshold > g_.length(first) and
            barcode_extractor_.GetMinPos(second, barcode) < tail_threshold and
            barcode_extractor_.GetNumberOfReads(first, barcode) >= count_threshold and
            barcode_extractor_.GetNumberOfReads(second, barcode) >= count_threshold;
    }

    DECL_LOGGER("BarcodePutChecker")
};

static CompositePutChecker<Graph>
CreateLongGapCloserPutChecker(const Graph& g, const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage,
                              const barcode_index::FrameBarcodeIndexInfoExtractor& barcode_extractor,
                              const vector<barcode_index::BarcodeId>& barcode_intersection,
                              size_t barcode_threshold, size_t count_threshold, size_t tail_threshold, size_t len_threshold) {
    auto barcode_put_checker = make_shared<BarcodedPathPutChecker<Graph>>(g, barcode_extractor, barcode_intersection,
                                                                          barcode_threshold, count_threshold,
                                                                          tail_threshold, len_threshold);
    auto unique_put_checker = make_shared<UniquePutChecker<Graph>>(g, unique_storage);
    vector<shared_ptr<VertexPutChecker<Graph>>> put_checkers({unique_put_checker, barcode_put_checker});
    return CompositePutChecker<Graph>(put_checkers);
}

typedef ComposedDijkstraSettings<Graph,
                                 LengthCalculator<Graph>,
                                 BoundProcessChecker<Graph>,
                                 CompositePutChecker<Graph>,
                                 ForwardNeighbourIteratorFactory<Graph> > LongGapCloserDijkstraSettings;

typedef Dijkstra<Graph, LongGapCloserDijkstraSettings> LongGapCloserDijkstra;

static LongGapCloserDijkstra CreateLongGapCloserDijkstra(const Graph& graph, size_t length_bound,
                                                         const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage,
                                                         const barcode_index::FrameBarcodeIndexInfoExtractor& barcode_extractor,
                                                         const vector<barcode_index::BarcodeId>& barcode_intersection,
                                                         size_t barcode_threshold, size_t count_threshold, size_t tail_threshold,
                                                         size_t len_threshold, size_t max_vertex_number = -1ul) {
    return LongGapCloserDijkstra(graph,
                          LongGapCloserDijkstraSettings(
                              LengthCalculator<Graph>(graph),
                              BoundProcessChecker<Graph>(length_bound),
                              CreateLongGapCloserPutChecker(graph, unique_storage, barcode_extractor, barcode_intersection,
                                                            barcode_threshold, count_threshold, tail_threshold,
                                                            len_threshold),
                              ForwardNeighbourIteratorFactory<Graph>(graph)),
                          max_vertex_number);
}

}
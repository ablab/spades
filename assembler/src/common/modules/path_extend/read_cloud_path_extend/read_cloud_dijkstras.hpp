#pragma once

#include "common/modules/path_extend/scaffolder2015/connection_condition2015.hpp"
#include "common/assembly_graph/graph_support/scaff_supplementary.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/intermediate_scaffolding/scaffold_vertex_predicates.hpp"

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


//------------------------------
// barcode based dijkstra
//------------------------------

template<class Graph, typename distance_t = size_t>
class BarcodedPathPutChecker : public VertexPutChecker<Graph, distance_t> {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    const Graph& g_;
    const path_extend::LongEdgePairGapCloserPredicate gap_closer_predicate_;

 public:
    BarcodedPathPutChecker(const Graph& g,
                           shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> barcode_extractor,
                           const shared_ptr<path_extend::PairEntryProcessor> pair_entry_processor,
                           const path_extend::scaffold_graph::ScaffoldVertex& start,
                           const path_extend::scaffold_graph::ScaffoldVertex& end,
                           const path_extend::LongEdgePairGapCloserParams& params) :
        VertexPutChecker<Graph, distance_t>(), g_(g),
        gap_closer_predicate_(g, barcode_extractor, params, start, end, pair_entry_processor) {}

    bool Check(VertexId, EdgeId edge, distance_t ) const override {
        DEBUG("Checking edge " << edge.int_id());
        DEBUG("Length " << g_.length(edge))
        return gap_closer_predicate_.Check(edge);
    }
    DECL_LOGGER("BarcodePutChecker")
};

template<class Graph, typename distance_t = size_t>
class SimpleBarcodePutChecker : public VertexPutChecker<Graph, distance_t> {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    const Graph& g_;
    shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> barcode_extractor_;
    barcode_index::SimpleVertexEntry target_barcodes_;
    size_t threshold_;

 public:
    SimpleBarcodePutChecker(const Graph& g,
                            shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> barcode_extractor,
                            const barcode_index::SimpleVertexEntry& target_barcodes, size_t threshold) :
        VertexPutChecker<Graph, distance_t>(),
        g_(g),
        barcode_extractor_(barcode_extractor),
        target_barcodes_(target_barcodes),
        threshold_(threshold) {}

    bool Check(VertexId, EdgeId edge, distance_t ) const override {
        DEBUG("Checking edge " << edge.int_id());
        DEBUG("Length " << g_.length(edge))
        size_t intersection_size = barcode_extractor_->GetIntersectionSize(edge, target_barcodes_);
        return intersection_size >= threshold_;
    }
    DECL_LOGGER("SimpleBarcodePutChecker");
};

class ReadCloudDijkstraHelper {

 public:
    typedef ComposedDijkstraSettings<Graph,
                                     LengthCalculator<Graph>,
                                     BoundProcessChecker<Graph>,
                                     UniquePutChecker<Graph>,
                                     ForwardNeighbourIteratorFactory<Graph> > UniqueDijkstraSettings;

    typedef Dijkstra<Graph, UniqueDijkstraSettings> UniqueDijkstra;

    static UniqueDijkstra CreateUniqueDijkstra(const Graph& graph, size_t distance_bound,
                                               const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage,
                                               size_t max_vertex_number = -1ul) {
        return UniqueDijkstra(graph,
                              UniqueDijkstraSettings(
                                  LengthCalculator<Graph>(graph),
                                  BoundProcessChecker<Graph>(distance_bound),
                                  UniquePutChecker<Graph>(graph, unique_storage),
                                  ForwardNeighbourIteratorFactory<Graph>(graph)),
                              max_vertex_number);
    }

    static CompositePutChecker<Graph>
    CreateLongGapCloserPutChecker(const Graph& g, const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage,
                                  shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> short_edge_extractor,
                                  shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> long_edge_extractor,
                                  const path_extend::scaffold_graph::ScaffoldVertex& start,
                                  const path_extend::scaffold_graph::ScaffoldVertex& end,
                                  const path_extend::LongEdgePairGapCloserParams& params) {
        auto score_function = make_shared<path_extend::RepetitiveVertexEntryScoreFunction>(short_edge_extractor);
        auto pair_entry_processor = make_shared<path_extend::TwoSetsBasedPairEntryProcessor>(
            long_edge_extractor->GetTailEntry(start),
            long_edge_extractor->GetHeadEntry(end),
            score_function);
//    auto pair_entry_processor = make_shared<path_extend::IntersectionBasedPairEntryProcessor>(barcode_intersection,
//                                                                                              short_edge_extractor);
        auto barcode_put_checker = make_shared<BarcodedPathPutChecker<Graph>>(g, short_edge_extractor, pair_entry_processor,
                                                                              start, end, params);
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

    static LongGapCloserDijkstra CreateLongGapCloserDijkstra(
            const Graph& graph, size_t length_bound,
            const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage,
            shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> short_edge_extractor,
            shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> long_edge_extractor,
            const path_extend::scaffold_graph::ScaffoldVertex& start,
            const path_extend::scaffold_graph::ScaffoldVertex& end,
            const path_extend::LongEdgePairGapCloserParams& vertex_predicate_params,
            size_t max_vertex_number = -1ul) {
        return LongGapCloserDijkstra(graph,
                                     LongGapCloserDijkstraSettings(
                                         LengthCalculator<Graph>(graph),
                                         BoundProcessChecker<Graph>(length_bound),
                                         CreateLongGapCloserPutChecker(graph, unique_storage, short_edge_extractor, long_edge_extractor,
                                                                       start, end, vertex_predicate_params),
                                         ForwardNeighbourIteratorFactory<Graph>(graph)),
                                     max_vertex_number);
    }

    typedef ComposedDijkstraSettings<Graph,
                                     LengthCalculator<Graph>,
                                     BoundProcessChecker<Graph>,
                                     CompositePutChecker<Graph>,
                                     ForwardNeighbourIteratorFactory<Graph> > LengthBoundedDijkstraSettings;

    typedef Dijkstra<Graph, LengthBoundedDijkstraSettings> LengthBoundedDijkstra;

    static LengthBoundedDijkstra CreateSimpleCloudBoundedDijkstra(
            const Graph& graph,
            shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> barcode_extractor,
            const barcode_index::SimpleVertexEntry target_barcodes,
            size_t length_threshold,
            size_t distance_bound,
            size_t barcode_threshold,
            size_t max_vertex_number = -1ul) {
        auto length_put_checker = std::make_shared<LengthPutChecker<Graph>>(graph, length_threshold);
        auto bound_put_checker = std::make_shared<BoundPutChecker<Graph>>(distance_bound);
        auto simple_barcode_checker = std::make_shared<SimpleBarcodePutChecker<Graph>>(graph,
                                                                                       barcode_extractor,
                                                                                       target_barcodes,
                                                                                       barcode_threshold);
        vector<shared_ptr<VertexPutChecker<Graph>>> put_checkers;
        put_checkers.push_back(length_put_checker);
        put_checkers.push_back(bound_put_checker);
        put_checkers.push_back(simple_barcode_checker);
        return LengthBoundedDijkstra(graph, LengthBoundedDijkstraSettings(
            LengthCalculator<Graph>(graph),
            BoundProcessChecker<Graph>(distance_bound),
            CompositePutChecker<Graph>(put_checkers),
            ForwardNeighbourIteratorFactory<Graph>(graph)),
                                     max_vertex_number);
    }
};
}
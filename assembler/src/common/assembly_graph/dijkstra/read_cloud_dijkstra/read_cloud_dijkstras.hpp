#pragma once

#include "common/modules/path_extend/scaffolder2015/connection_condition2015.hpp"
#include "common/assembly_graph/graph_support/scaff_supplementary.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/intermediate_scaffolding/scaffold_vertex_predicates.hpp"

namespace omnigraph {
//------------------------------
// unique based dijkstra
//------------------------------

template<class Graph, typename distance_t = size_t>
class UniquePutChecker {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    const Graph &g_;
    const path_extend::ScaffoldingUniqueEdgeStorage &unique_storage_;

  public:
    UniquePutChecker(const Graph &g,
                     const path_extend::ScaffoldingUniqueEdgeStorage &unique_storage) :
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
class BarcodedPathPutChecker {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    const Graph &g_;
    const path_extend::read_cloud::LongEdgePairGapCloserPredicate gap_closer_predicate_;

  public:
    BarcodedPathPutChecker(const Graph &g,
                           shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> barcode_extractor,
                           const shared_ptr<path_extend::read_cloud::PairEntryProcessor> pair_entry_processor,
                           const path_extend::scaffold_graph::ScaffoldVertex &start,
                           const path_extend::scaffold_graph::ScaffoldVertex &end,
                           const path_extend::read_cloud::LongEdgePairGapCloserParams &params) :
        g_(g), gap_closer_predicate_(g, barcode_extractor, params, start, end, pair_entry_processor) {}

    bool Check(VertexId, EdgeId edge, distance_t) const {
        DEBUG("Checking edge " << edge.int_id());
        DEBUG("Length " << g_.length(edge))
        return gap_closer_predicate_.Check(edge);
    }
    DECL_LOGGER("BarcodePutChecker")
};

template<class Graph, typename distance_t = size_t>
class SimpleBarcodePutChecker {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    const Graph &g_;
    shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> barcode_extractor_;
    barcode_index::SimpleVertexEntry target_barcodes_;
    size_t threshold_;

  public:
    SimpleBarcodePutChecker(const Graph &g,
                            shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> barcode_extractor,
                            const barcode_index::SimpleVertexEntry &target_barcodes, size_t threshold) :
        g_(g),
        barcode_extractor_(barcode_extractor),
        target_barcodes_(target_barcodes),
        threshold_(threshold) {}

    bool Check(VertexId, EdgeId edge, distance_t) const {
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
    typedef barcode_index::SimpleIntersectingScaffoldVertexExtractor BarcodeExtractor;

    static UniqueDijkstra CreateUniqueDijkstra(const Graph &graph, size_t distance_bound,
                                               const path_extend::ScaffoldingUniqueEdgeStorage &unique_storage,
                                               size_t max_vertex_number = -1ul) {
        return UniqueDijkstra(graph,
                              UniqueDijkstraSettings(
                                  LengthCalculator<Graph>(graph),
                                  BoundProcessChecker<Graph>(distance_bound),
                                  UniquePutChecker<Graph>(graph, unique_storage),
                                  ForwardNeighbourIteratorFactory<Graph>(graph)),
                              max_vertex_number);
    }

    typedef std::tuple<BarcodedPathPutChecker<Graph>, UniquePutChecker<Graph>> GapCloserCheckerTuple;

    static AndPutChecker<GapCloserCheckerTuple, Graph>
    CreateLongGapCloserPutChecker(const Graph &g, const path_extend::ScaffoldingUniqueEdgeStorage &unique_storage,
                                  shared_ptr<BarcodeExtractor> short_edge_extractor,
                                  shared_ptr<BarcodeExtractor> long_edge_extractor,
                                  const path_extend::scaffold_graph::ScaffoldVertex &start,
                                  const path_extend::scaffold_graph::ScaffoldVertex &end,
                                  const path_extend::read_cloud::LongEdgePairGapCloserParams &params) {
        auto score_function =
            make_shared<path_extend::read_cloud::RepetitiveVertexEntryScoreFunction>(short_edge_extractor);
        auto pair_entry_processor = make_shared<path_extend::read_cloud::TwoSetsBasedPairEntryProcessor>(
            long_edge_extractor->GetTailEntry(start),
            long_edge_extractor->GetHeadEntry(end),
            score_function);
//    auto pair_entry_processor = make_shared<path_extend::IntersectionBasedPairEntryProcessor>(barcode_intersection,
//                                                                                              short_edge_extractor);
        BarcodedPathPutChecker<Graph> barcode_put_checker(g, short_edge_extractor, pair_entry_processor,
                                                          start, end, params);
        UniquePutChecker<Graph> unique_put_checker(g, unique_storage);

        auto put_checkers = std::make_tuple(barcode_put_checker, unique_put_checker);
        return AndPutChecker<GapCloserCheckerTuple, Graph>(put_checkers);
    }

    typedef ComposedDijkstraSettings<Graph,
                                     LengthCalculator<Graph>,
                                     BoundProcessChecker<Graph>,
                                     AndPutChecker<GapCloserCheckerTuple, Graph>,
                                     ForwardNeighbourIteratorFactory<Graph> > LongGapCloserDijkstraSettings;

    typedef Dijkstra<Graph, LongGapCloserDijkstraSettings> LongGapCloserDijkstra;

    static LongGapCloserDijkstra CreateLongGapCloserDijkstra(
        const Graph &graph, size_t length_bound,
        const path_extend::ScaffoldingUniqueEdgeStorage &unique_storage,
        shared_ptr<BarcodeExtractor> short_edge_extractor,
        shared_ptr<BarcodeExtractor> long_edge_extractor,
        const path_extend::scaffold_graph::ScaffoldVertex &start,
        const path_extend::scaffold_graph::ScaffoldVertex &end,
        const path_extend::read_cloud::LongEdgePairGapCloserParams &vertex_predicate_params,
        size_t max_vertex_number = -1ul) {
        return LongGapCloserDijkstra(graph,
                                     LongGapCloserDijkstraSettings(
                                         LengthCalculator<Graph>(graph),
                                         BoundProcessChecker<Graph>(length_bound),
                                         CreateLongGapCloserPutChecker(graph, unique_storage, short_edge_extractor,
                                                                       long_edge_extractor,
                                                                       start, end, vertex_predicate_params),
                                         ForwardNeighbourIteratorFactory<Graph>(graph)),
                                     max_vertex_number);
    }

    typedef std::tuple<LengthPutChecker<Graph>, BoundPutChecker<Graph>, SimpleBarcodePutChecker<Graph>>
        SimpleCloudTuple;
    typedef ComposedDijkstraSettings<Graph,
                                     LengthCalculator<Graph>,
                                     BoundProcessChecker<Graph>,
                                     AndPutChecker<SimpleCloudTuple, Graph>,
                                     ForwardNeighbourIteratorFactory<Graph> > LengthBoundedDijkstraSettings;

    typedef Dijkstra<Graph, LengthBoundedDijkstraSettings> LengthBoundedDijkstra;

    static LengthBoundedDijkstra CreateSimpleCloudBoundedDijkstra(
        const Graph &graph,
        shared_ptr<BarcodeExtractor> barcode_extractor,
        const barcode_index::SimpleVertexEntry target_barcodes,
        size_t length_threshold,
        size_t distance_bound,
        size_t barcode_threshold,
        size_t max_vertex_number = -1ul) {
        LengthPutChecker<Graph> length_put_checker(graph, length_threshold);
        BoundPutChecker<Graph> bound_put_checker(distance_bound);
        SimpleBarcodePutChecker<Graph> simple_barcode_checker(graph, barcode_extractor, target_barcodes,
                                                              barcode_threshold);
        auto put_checkers = std::make_tuple(length_put_checker, bound_put_checker, simple_barcode_checker);
        return LengthBoundedDijkstra(graph, LengthBoundedDijkstraSettings(
            LengthCalculator<Graph>(graph),
            BoundProcessChecker<Graph>(distance_bound),
            AndPutChecker<SimpleCloudTuple, Graph>(put_checkers),
            ForwardNeighbourIteratorFactory<Graph>(graph)),
                                     max_vertex_number);
    }
};
}
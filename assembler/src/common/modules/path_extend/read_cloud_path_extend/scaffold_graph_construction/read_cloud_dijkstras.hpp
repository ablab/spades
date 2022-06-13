//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/dijkstra/dijkstra_settings.hpp"
#include "assembly_graph/dijkstra/dijkstra_algorithm.hpp"
#include "modules/path_extend//scaff_supplementary.hpp"
#include "modules/path_extend/read_cloud_path_extend/intermediate_scaffolding/pair_entry_processors.hpp"

namespace path_extend {

namespace read_cloud {

using namespace omnigraph;

//------------------------------
// unique based dijkstra
//------------------------------

template<class Graph, typename distance_t = size_t>
class UniquePutChecker {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    const Graph &g_;
    const ScaffoldingUniqueEdgeStorage &unique_storage_;

  public:
    UniquePutChecker(const Graph &g,
                     const ScaffoldingUniqueEdgeStorage &unique_storage) :
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
            DEBUG("Short edge, passed" << std::endl)  //todo use short edges to reduce number of candidates
            return true;
        }
        if (!unique_storage_.IsUnique(edge)) {
            DEBUG("Long non-unique nearby edge, passed" << std::endl)
            return true;
        } else {
            DEBUG("Long unique nearby edge, put to candidates list and stop" << std::endl)
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
    const LongEdgePairGapCloserPredicate gap_closer_predicate_;

  public:
    BarcodedPathPutChecker(const Graph &g,
                           std::shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> barcode_extractor,
                           std::shared_ptr<PairEntryProcessor> pair_entry_processor,
                           const scaffold_graph::ScaffoldVertex &start,
                           const scaffold_graph::ScaffoldVertex &end,
                           const LongEdgePairGapCloserParams &params) :
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
    std::shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> barcode_extractor_;
    barcode_index::SimpleVertexEntry target_barcodes_;
    size_t threshold_;

  public:
    SimpleBarcodePutChecker(const Graph &g,
                            std::shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> barcode_extractor,
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
    typedef debruijn_graph::Graph Graph;

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
    CreateLongGapCloserPutChecker(const Graph &g, const ScaffoldingUniqueEdgeStorage &unique_storage,
                                  std::shared_ptr<BarcodeExtractor> short_edge_extractor,
                                  std::shared_ptr<BarcodeExtractor> long_edge_extractor,
                                  const scaffold_graph::ScaffoldVertex &start,
                                  const scaffold_graph::ScaffoldVertex &end,
                                  const LongEdgePairGapCloserParams &params) {
        auto score_function =
            std::make_shared<RepetitiveVertexEntryScoreFunction>(short_edge_extractor);
        auto pair_entry_processor = std::make_shared<TwoSetsBasedPairEntryProcessor>(
            long_edge_extractor->GetTailEntry(start),
            long_edge_extractor->GetHeadEntry(end),
            score_function);
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
        const ScaffoldingUniqueEdgeStorage &unique_storage,
        std::shared_ptr<BarcodeExtractor> short_edge_extractor,
        std::shared_ptr<BarcodeExtractor> long_edge_extractor,
        const scaffold_graph::ScaffoldVertex &start,
        const scaffold_graph::ScaffoldVertex &end,
        const LongEdgePairGapCloserParams &vertex_predicate_params,
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
        std::shared_ptr<BarcodeExtractor> barcode_extractor,
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
}
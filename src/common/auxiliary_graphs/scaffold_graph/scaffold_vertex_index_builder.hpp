//***************************************************************************
//* Copyright (c) 2023-2025 SPAdes team
//* Copyright (c) 2019-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "scaffold_vertex_index.hpp"
#include "barcode_info_extractor.hpp"

namespace barcode_index {

    template <class EdgeEntryT>
    class AbstractScaffoldVertexEntryExtractor {
     public:
        virtual EdgeEntryT ExtractEntry(const scaffold_graph::ScaffoldVertex &vertex) const = 0;
    };

    class TailThresholdGetter {
     public:
        virtual size_t GetTailThreshold(const scaffold_graph::ScaffoldVertex &vertex) const = 0;
    };

    class ConstTailThresholdGetter: public TailThresholdGetter {
     public:
        explicit ConstTailThresholdGetter(const size_t tail_threshold_) : tail_threshold_(tail_threshold_) {}
        size_t GetTailThreshold(const scaffold_graph::ScaffoldVertex &/*vertex*/) const override {
            return tail_threshold_;
        }
     private:
        const size_t tail_threshold_;
    };

    class FractionTailThresholdGetter: public TailThresholdGetter {
     public:
        typedef debruijn_graph::Graph Graph;

        FractionTailThresholdGetter(const Graph &g_, const double edge_length_fraction_)
            : g_(g_), edge_length_fraction_(edge_length_fraction_) {}

        size_t GetTailThreshold(const scaffold_graph::ScaffoldVertex &vertex) const override {
            return static_cast<size_t>(static_cast<double>(vertex.GetLengthFromGraph(g_)) * edge_length_fraction_);
        }
     private:
        const Graph& g_;
        const double edge_length_fraction_;
    };

    class ScaffoldVertexSimpleEntryExtractor: public AbstractScaffoldVertexEntryExtractor<SimpleVertexEntry> {
     public:
        typedef debruijn_graph::Graph Graph;
        typedef debruijn_graph::EdgeId EdgeId;
        typedef typename scaffold_graph::EdgeIdVertex EdgeIdVertex;
        typedef typename scaffold_graph::PathVertex PathVertex;

        ScaffoldVertexSimpleEntryExtractor(const Graph &g_,
                                           const FrameBarcodeIndexInfoExtractor &barcode_extractor_,
                                           std::shared_ptr<TailThresholdGetter> tail_threshold_getter,
                                           const size_t count_threshold_,
                                           const size_t length_threshold_)
            : g_(g_),
              barcode_extractor_(barcode_extractor_),
              tail_threshold_getter_(tail_threshold_getter),
              count_threshold_(count_threshold_),
              length_threshold_(length_threshold_) {}

        SimpleVertexEntry ExtractEntry(const scaffold_graph::ScaffoldVertex &vertex) const override {
            auto inner_vertex = vertex.GetInnerVertex();

            SimpleVertexEntry empty;
            auto type = vertex.GetType();
            switch (type) {
                case scaffold_graph::ScaffoldVertexT::Edge: {
                    auto edge_vertex = std::static_pointer_cast<EdgeIdVertex>(inner_vertex);
                    return ExtractEntryInner(edge_vertex);
                }
                case scaffold_graph::ScaffoldVertexT::Path: {
                    auto path_vertex = std::static_pointer_cast<PathVertex>(inner_vertex);
                    return ExtractEntryInner(path_vertex);
                }
            }
            WARN("ScaffoldVertex of unknown type");
            return empty;
        }

     private:
        SimpleVertexEntry ExtractEntryInner(std::shared_ptr<EdgeIdVertex> simple_edge_vertex) const {
            SimpleVertexEntry result;
            TRACE("Extracting entry from edge");
            auto edge = simple_edge_vertex->get();
            size_t tail_threshold = tail_threshold_getter_->GetTailThreshold(edge);
            TRACE("Tail threshold: " << tail_threshold);
            auto entry = barcode_extractor_.GetBarcodesFromHead(edge, count_threshold_, tail_threshold);
            std::copy(entry.begin(), entry.end(), std::inserter(result, result.end()));
            TRACE("Entry size: " << entry.size());
            return result;
        }

        //fixme optimize later
        SimpleVertexEntry ExtractEntryInner(std::shared_ptr<PathVertex> path_vertex) const {
            TRACE("Extracting entry from path");
            size_t current_prefix = 0;
            path_extend::BidirectionalPath* path = path_vertex->get();
            size_t path_size = path->Size();
            SimpleVertexEntry result;
            const size_t global_count_threshold = 5;
            std::unordered_map<BarcodeId, size_t> barcode_to_count;
            size_t tail_threshold = tail_threshold_getter_->GetTailThreshold(path_vertex->get());
            TRACE("Tail threshold: " << tail_threshold);
            for (size_t i = 0; i < path_size and current_prefix <= tail_threshold; ++i) {
                EdgeId current_edge = path->At(i);
                if (g_.length(current_edge) < length_threshold_) {
                    current_prefix += g_.length(current_edge);
                    continue;
                }
                size_t current_tail = tail_threshold - current_prefix;
                TRACE("Current tail: " << current_tail);
                const auto &current_entry = barcode_extractor_.GetBarcodesAndCountsFromHead(current_edge,
                                                                                            count_threshold_,
                                                                                            current_tail);
                for (const auto& barcode_and_reads: current_entry) {
                    barcode_to_count[barcode_and_reads.first] += barcode_and_reads.second;
                }
                TRACE("Current entry size: " << barcode_to_count.size());
                current_prefix += g_.length(current_edge);
            }
            for (const auto& barcode_and_count: barcode_to_count) {
                if (barcode_and_count.second >= global_count_threshold) {
                    result.insert(barcode_and_count.first);
                }
            }
            TRACE("Result size: " << result.size());
            return result;
        }

        const debruijn_graph::Graph &g_;
        const FrameBarcodeIndexInfoExtractor& barcode_extractor_;
        std::shared_ptr<TailThresholdGetter> tail_threshold_getter_;
        const size_t count_threshold_;
        const size_t length_threshold_;

        DECL_LOGGER("ScaffoldVertexSimpleEntryExtractor");
    };

    template <class EdgeEntryT>
    class ScaffoldVertexIndexBuilder {
     public:
        typedef debruijn_graph::Graph Graph;
        typedef debruijn_graph::EdgeId EdgeId;
        typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
        typedef std::shared_ptr<AbstractScaffoldVertexEntryExtractor<EdgeEntryT>> EntryExtractorPtr;

        ScaffoldVertexIndexBuilder(const Graph &g, EntryExtractorPtr vertex_entry_extractor, size_t max_threads)
            : g_(g), vertex_entry_extractor_(vertex_entry_extractor),
              index_(std::make_shared<ScaffoldVertexIndex<EdgeEntryT>>(g_)), max_threads_(max_threads) {}

        template <class ContainerT>
        std::shared_ptr<ScaffoldVertexIndex<EdgeEntryT>> GetConstructedIndex(const ContainerT& vertex_container) {

            //todo make parallel using iterator chunks
            DEBUG("Constructing long edge index in " << max_threads_ << " threads");
//            size_t counter = 0;
//            size_t block_size = vertex_container.size() / 10;
            for (const auto& vertex: vertex_container)
            {
                auto entry = vertex_entry_extractor_->ExtractEntry(vertex);
                TRACE("Entry size: " << entry.size());
                {
                    index_->InsertEntry(vertex, std::move(entry));
//                    ++counter;
                }
//                if (counter % block_size == 0) {
//                    INFO("Processed " << counter << " edges out of " << vertex_container.size());
//                }
            }
            DEBUG("Constructed long edge index");
            return index_;
        }

     private:
        const Graph& g_;
        EntryExtractorPtr vertex_entry_extractor_;
        std::shared_ptr<ScaffoldVertexIndex<EdgeEntryT>> index_;
        size_t max_threads_;
    };

    class SimpleScaffoldVertexIndexBuilderHelper {
     public:
        typedef debruijn_graph::Graph Graph;
        typedef std::shared_ptr<SimpleScaffoldVertexIndex> ScaffoldIndexPtr;
        typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;

        template <class ContainerT>
        ScaffoldIndexPtr ConstructScaffoldVertexIndex(const Graph& g_, const FrameBarcodeIndexInfoExtractor& extractor,
                                                      std::shared_ptr<TailThresholdGetter> tail_threshold_getter,
                                                      size_t count_threshold, size_t length_threshold,
                                                      size_t max_threads, const ContainerT& vertex_container) {
            DEBUG("Building simple long edge barcode index with parameters");
            DEBUG("Count threshold: " << count_threshold);
            DEBUG("Length threshold: " << length_threshold);
            auto entry_extractor = std::make_shared<ScaffoldVertexSimpleEntryExtractor>(g_, extractor,
                                                                                        tail_threshold_getter,
                                                                                        count_threshold,
                                                                                        length_threshold);
            ScaffoldVertexIndexBuilder<SimpleVertexEntry> builder(g_, entry_extractor, max_threads);
            return builder.GetConstructedIndex(vertex_container);
        }

        template <class ContainerT>
        ScaffoldIndexPtr HalfEdgeScaffoldVertexIndex(const Graph& g, const FrameBarcodeIndexInfoExtractor& extractor,
                                                     const ContainerT& vertex_container, size_t count_threshold,
                                                     size_t max_threads) {
            const size_t length_threshold = 1000;
            const size_t linkage_distance = 10;
            const double EDGE_LENGTH_FRACTION = 0.5;
            auto threshold_getter = std::make_shared<barcode_index::FractionTailThresholdGetter>(g, EDGE_LENGTH_FRACTION);
            auto split_scaffold_vertex_index = ConstructScaffoldVertexIndex(g, extractor,
                                                                            threshold_getter,
                                                                            count_threshold, length_threshold,
                                                                            max_threads, vertex_container);
            return split_scaffold_vertex_index;
        }

        template <class ContainerT>
        ScaffoldIndexPtr TailEdgeScaffoldVertexIndex(const Graph& g, const FrameBarcodeIndexInfoExtractor& extractor,
                                                     const ContainerT& vertex_container, size_t count_threshold,
                                                     size_t tail_threshold, size_t max_threads) {
            const size_t length_threshold = 1000;
            auto tail_threshold_getter = std::make_shared<barcode_index::ConstTailThresholdGetter>(tail_threshold);
            auto scaffold_vertex_index = ConstructScaffoldVertexIndex(g, extractor, tail_threshold_getter,
                                                                      count_threshold, length_threshold,
                                                                      max_threads, vertex_container);
            return scaffold_vertex_index;
        }
    };

}
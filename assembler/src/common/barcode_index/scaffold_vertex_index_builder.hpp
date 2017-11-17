#pragma once
#include "scaffold_vertex_index.hpp"
#include "barcode_info_extractor.hpp"

namespace barcode_index {

    template <class EdgeEntryT>
    class AbstractScaffoldVertexEntryExtractor {
     public:
        virtual EdgeEntryT ExtractEntry(const path_extend::scaffold_graph::ScaffoldVertex &vertex) const = 0;
    };

    class SimpleScaffoldVertexEntryExtractor: public AbstractScaffoldVertexEntryExtractor<SimpleVertexEntry> {
     public:
        typedef typename path_extend::scaffold_graph::EdgeIdVertex EdgeIdVertex;
        typedef typename path_extend::scaffold_graph::PathVertex PathVertex;
     private:
        const debruijn_graph::Graph &g_;
        const FrameBarcodeIndexInfoExtractor& barcode_extractor_;
        const size_t tail_threshold_;
        const size_t count_threshold_;
        const size_t length_threshold_;

     public:
        SimpleScaffoldVertexEntryExtractor(const Graph &g_,
                                           const FrameBarcodeIndexInfoExtractor &barcode_extractor_,
                                           const size_t tail_threshold_,
                                           const size_t count_threshold_,
                                           const size_t length_threshold_)
            : g_(g_),
              barcode_extractor_(barcode_extractor_),
              tail_threshold_(tail_threshold_),
              count_threshold_(count_threshold_),
              length_threshold_(length_threshold_) {}

        SimpleVertexEntry ExtractEntry(const path_extend::scaffold_graph::ScaffoldVertex &vertex) const override {
            auto inner_vertex = vertex.getInnerVertex();

            SimpleVertexEntry empty;
            auto type = vertex.getType();
            switch (type) {
                case path_extend::scaffold_graph::ScaffoldVertexT::Edge: {
                    auto edge_vertex = std::static_pointer_cast<EdgeIdVertex>(inner_vertex);
                    return ExtractEntryInner(edge_vertex);
                }
                case path_extend::scaffold_graph::ScaffoldVertexT::Path: {
                    auto path_vertex = std::static_pointer_cast<PathVertex>(inner_vertex);
                    return ExtractEntryInner(path_vertex);
                }
            }
            WARN("ScaffoldVertex of unknown type");
            return empty;
        }

     private:
        SimpleVertexEntry ExtractEntryInner(shared_ptr<EdgeIdVertex> simple_edge_vertex) const {
            SimpleVertexEntry result;
            const auto& entry = barcode_extractor_.GetBarcodesFromRange(simple_edge_vertex->get(), count_threshold_, 0, tail_threshold_);
            std::copy(entry.begin(), entry.end(), std::inserter(result, result.end()));
            return result;
        }

        //fixme optimize later
        SimpleVertexEntry ExtractEntryInner(shared_ptr<PathVertex> path_vertex) const {
            size_t current_prefix = 0;
            path_extend::BidirectionalPath* path = path_vertex->get();
            size_t path_size = path->Size();
            SimpleVertexEntry old_result;
            SimpleVertexEntry new_result;
            for (size_t i = 0; i < path_size and current_prefix <= tail_threshold_; ++i) {
                EdgeId current_edge = path->At(i);
                if (g_.length(current_edge) < length_threshold_) {
                    continue;
                }
                size_t current_tail = tail_threshold_ - current_prefix;
                const auto &current_entry = barcode_extractor_.GetBarcodesFromRange(current_edge, count_threshold_, 0, current_tail);
                std::set_union(old_result.begin(), old_result.end(), current_entry.begin(),
                               current_entry.end(), std::inserter(new_result, new_result.end()));
                old_result.clear();
                std::swap(old_result, new_result);
            }
            return new_result;
        }
    };

    template <class EdgeEntryT>
    class ScaffoldVertexIndexBuilder {
        typedef path_extend::scaffold_graph::ScaffoldVertex ScaffoldVertex;

        const Graph& g_;
        shared_ptr<AbstractScaffoldVertexEntryExtractor<EdgeEntryT>> vertex_entry_extractor_;
        shared_ptr<ScaffoldVertexIndex<EdgeEntryT>> index_;
        size_t max_threads_;

     public:
        ScaffoldVertexIndexBuilder(const Graph &g_,
                                   shared_ptr<AbstractScaffoldVertexEntryExtractor<EdgeEntryT>> vertex_entry_extractor_,
                                   size_t max_threads)
            : g_(g_), vertex_entry_extractor_(vertex_entry_extractor_),
              index_(std::make_shared<ScaffoldVertexIndex<EdgeEntryT>>(g_)), max_threads_(max_threads) {}

        template <class ContainerT>
        shared_ptr<ScaffoldVertexIndex<EdgeEntryT>> GetConstructedIndex(const ContainerT& vertex_container) {
            vector<ScaffoldVertex> vertices_copy;
            std::copy(vertex_container.begin(), vertex_container.end(), std::back_inserter(vertices_copy));
            INFO("Constructing long edge index in " << max_threads_ << " threads");
            size_t counter = 0;
            const size_t block_size = vertices_copy.size() / 10;
#pragma omp parallel for num_threads(max_threads_) 
            for (size_t i = 0; i < vertices_copy.size(); ++i)
            {
                const ScaffoldVertex& vertex = vertices_copy[i];
                auto entry = vertex_entry_extractor_->ExtractEntry(vertex);
#pragma omp critical
                {
                    index_->InsertEntry(vertex, std::move(entry));
                    ++counter;
                }
                if (counter % block_size == 0) {
                    INFO("Processed " << counter << " edges out of " << vertices_copy.size());
                }
            }
            INFO("Constructed long edge index");
            return index_;
        }
    };

    class SimpleScaffoldVertexIndexBuilderHelper {
     public:
        typedef path_extend::scaffold_graph::ScaffoldVertex ScaffoldVertex;

        template <class ContainerT>
        shared_ptr<SimpleScaffoldVertexIndex> ConstructScaffoldVertexIndex(const Graph& g_,
                                                                           const FrameBarcodeIndexInfoExtractor& extractor,
                                                                           size_t tail_threshold,
                                                                           size_t count_threshold,
                                                                           size_t length_threshold,
                                                                           size_t max_threads,
                                                                           const ContainerT& vertex_container) {
            auto entry_extractor = make_shared<SimpleScaffoldVertexEntryExtractor>(g_, extractor, tail_threshold,
                                                                                       count_threshold, length_threshold);
            ScaffoldVertexIndexBuilder<SimpleVertexEntry> builder(g_, entry_extractor, max_threads);
            return builder.GetConstructedIndex(vertex_container);
        }
    };

}
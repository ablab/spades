#pragma once

#include "modules/path_extend/scaffolder2015/scaffold_vertex.hpp"
#include "common/assembly_graph/core/graph.hpp"
#include "adt/iterator_range.hpp"
#include "barcode_info_extractor.hpp"

namespace barcode_index {

    template <class VertexEntryT>
    class ScaffoldVertexIndexBuilder;

    template <class VertexEntryT>
    class ScaffoldVertexIndex {
        friend class ScaffoldVertexIndexBuilder<VertexEntryT>;
     public:
        typedef path_extend::scaffold_graph::ScaffoldVertex ScaffoldVertex;
        typedef typename VertexEntryT::const_iterator const_iterator;

     private:
        const Graph& g_;
        std::unordered_map<ScaffoldVertex, VertexEntryT> vertex_to_entry_;

     public:
        ScaffoldVertexIndex(const Graph &g_): g_(g_), vertex_to_entry_() {}

        const VertexEntryT& GetHeadEntry(const ScaffoldVertex& vertex) const {
            return vertex_to_entry_.at(vertex);
        }

        const VertexEntryT& GetTailEntry(const ScaffoldVertex& vertex) const {
            return vertex_to_entry_.at(vertex.getConjugateFromGraph(g_));
        }

        const_iterator GetHeadBegin(const ScaffoldVertex& vertex) const {
            return vertex_to_entry_.at(vertex).begin();
        }

        const_iterator GetHeadEnd(const ScaffoldVertex& vertex) const {
            return vertex_to_entry_.at(vertex).end();
        }

        adt::iterator_range<const_iterator> GetHeadRange(const ScaffoldVertex& vertex) const {
            return adt::make_range(GetHeadBegin(vertex), GetHeadEnd(vertex));
        }

        const_iterator GetTailBegin(const ScaffoldVertex& vertex) const {
            return vertex_to_entry_.at(vertex.getConjugateFromGraph(g_)).begin();
        }

        const_iterator GetTailEnd(const ScaffoldVertex& vertex) const {
            return vertex_to_entry_.at(vertex.getConjugateFromGraph(g_)).end();
        }

        adt::iterator_range<const_iterator> GetTailRange(const ScaffoldVertex& vertex) const {
            return adt::make_range(GetTailBegin(vertex.getConjugateFromGraph(g_)),
                                   GetTailEnd(vertex.getConjugateFromGraph(g_)));
        }

     private:
        void InsertEntry(const ScaffoldVertex& vertex, VertexEntryT&& entry) {
            vertex_to_entry_.insert({vertex, entry});
        }
    };

    typedef std::set<BarcodeId> SimpleVertexEntry;

    typedef ScaffoldVertexIndex<SimpleVertexEntry> SimpleScaffoldVertexIndex;

    class ScaffoldVertexIndexInfoExtractor {
     public:
        typedef typename path_extend::scaffold_graph::ScaffoldVertex ScaffoldVertex;
     public:
        virtual size_t GetHeadSize(const ScaffoldVertex &vertex) const = 0;
        virtual size_t GetTailSize(const ScaffoldVertex &vertex) const = 0;

        virtual size_t GetIntersectionSize(const ScaffoldVertex &first, const ScaffoldVertex &second) const = 0;

        /**
         * @note second is supposed to be between first and third
         */
        virtual size_t GetIntersectionSize(const ScaffoldVertex &first, const ScaffoldVertex &second,
                                           const ScaffoldVertex &third) const = 0;
    };

    template <class VertexEntryT>
    class IntersectingScaffoldVertexIndexInfoExtractor: public ScaffoldVertexIndexInfoExtractor {
     public:
        using ScaffoldVertexIndexInfoExtractor::ScaffoldVertex;

     public:
        virtual SimpleVertexEntry GetIntersection(const VertexEntryT &first, const VertexEntryT &second) const = 0;

        virtual SimpleVertexEntry GetIntersection(const ScaffoldVertex &first, const ScaffoldVertex &second) const = 0;
        /**
         * @note second is supposed to be between first and third
         */
        virtual size_t GetIntersectionSize(const ScaffoldVertex &middle, const VertexEntryT &entry) const = 0;

        size_t GetIntersectionSize(const VertexEntryT &first, const VertexEntryT &second) {
            return GetIntersection(first, second).size();
        }

        virtual SimpleVertexEntry GetHeadEntry(const ScaffoldVertex &vertex) = 0;

        virtual SimpleVertexEntry GetTailEntry(const ScaffoldVertex &vertex) = 0;
    };

class BarcodeIndexInfoExtractorWrapper: public IntersectingScaffoldVertexIndexInfoExtractor<SimpleVertexEntry> {
 private:
    const Graph& g_;
    shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_;
 public:
    BarcodeIndexInfoExtractorWrapper(const Graph &g,shared_ptr<FrameBarcodeIndexInfoExtractor> barcode_index_)
        : g_(g), barcode_extractor_(barcode_index_) {}

    size_t GetHeadSize(const ScaffoldVertex &vertex) const override {
        return barcode_extractor_->GetNumberOfBarcodes(GetEdge(vertex));
    }
    size_t GetTailSize(const ScaffoldVertex &vertex) const override {
        return barcode_extractor_->GetNumberOfBarcodes(g_.conjugate(GetEdge(vertex)));
    }
    size_t GetIntersectionSize(const ScaffoldVertex &first, const ScaffoldVertex &second) const override {
        return barcode_extractor_->GetNumberOfSharedBarcodes(GetEdge(first), GetEdge(second));
    }
    size_t GetIntersectionSize(const ScaffoldVertex &first,
                               const ScaffoldVertex &second,
                               const ScaffoldVertex &third) const override {
        return GetIntersectionSize(third, GetIntersection(first, second));
    }

    SimpleVertexEntry GetIntersection(const SimpleVertexEntry &first,
                                      const SimpleVertexEntry &second) const override {
        SimpleVertexEntry result;
        std::set_intersection(first.begin(), first.end(), second.begin(), second.end(), std::inserter(result, result.end()));
        return result;
    }

    SimpleVertexEntry GetIntersection(const ScaffoldVertex &first, const ScaffoldVertex &second) const override {
        auto intersection = barcode_extractor_->GetSharedBarcodes(GetEdge(first), GetEdge(second));
        std::set<BarcodeId> result;
        std::copy(intersection.begin(), intersection.end(), std::inserter(result, result.begin()));
        return result;
    }
    size_t GetIntersectionSize(const ScaffoldVertex &middle, const SimpleVertexEntry &entry) const override {
        auto barcodes = barcode_extractor_->GetBarcodes(GetEdge(middle));
        SimpleVertexEntry intersection;
        std::set_intersection(barcodes.begin(), barcodes.end(), entry.begin(), entry.end(),
                              std::inserter(intersection, intersection.begin()));
        return intersection.size();
    }

    //fixme Can not collect from part of the edge. Slow.
    SimpleVertexEntry GetHeadEntry(const ScaffoldVertex &vertex) override {
        VERIFY_MSG(false, "Head entry extractor from BarcodeIndexInfoExtractorWrapper is currently not supported");
        SimpleVertexEntry result;
        auto edge = GetEdge(vertex);
        auto barcodes = barcode_extractor_->GetBarcodes(edge);
        std::copy(barcodes.begin(), barcodes.end(), std::inserter(result, result.begin()));
        return result;
    }

    SimpleVertexEntry GetTailEntry(const ScaffoldVertex &vertex) override {
        return GetHeadEntry(vertex);
    }

 private:
        EdgeId GetEdge(const ScaffoldVertex& vertex) const {
            path_extend::scaffold_graph::EdgeGetter edge_getter;
            return edge_getter.GetEdgeFromScaffoldVertex(vertex);
        }
    };

    class SimpleScaffoldVertexIndexInfoExtractor: public IntersectingScaffoldVertexIndexInfoExtractor<SimpleVertexEntry> {
        typedef path_extend::scaffold_graph::ScaffoldVertex ScaffoldVertex;
     private:
        shared_ptr<ScaffoldVertexIndex<SimpleVertexEntry>> index_;

     public:
        explicit SimpleScaffoldVertexIndexInfoExtractor(shared_ptr<ScaffoldVertexIndex<SimpleVertexEntry>> index_)
            : index_(index_) {}

        SimpleVertexEntry GetIntersection(const ScaffoldVertex &first, const ScaffoldVertex &second) const override {
            SimpleVertexEntry result;
            auto first_begin = index_->GetTailBegin(first);
            auto first_end = index_->GetTailEnd(first);
            auto second_begin = index_->GetHeadBegin(second);
            auto second_end = index_->GetHeadEnd(second);
            std::set_intersection(first_begin, first_end, second_begin, second_end, std::inserter(result, result.end()));
            return result;
        }

        size_t GetIntersectionSize(const ScaffoldVertex &first, const ScaffoldVertex &second) const override {
            return GetIntersection(first, second).size();
        }

        size_t GetIntersectionSize(const ScaffoldVertex &middle, const SimpleVertexEntry &entry) const override {
            auto middle_begin = index_->GetHeadBegin(middle);
            auto middle_end = index_->GetHeadEnd(middle);
            SimpleVertexEntry intersection;
            std::set_intersection(entry.begin(), entry.end(), middle_begin, middle_end, std::inserter(intersection, intersection.end()));
            return intersection.size();
        }

        size_t GetIntersectionSize(const ScaffoldVertex &first,
                                   const ScaffoldVertex &second,
                                   const ScaffoldVertex &third) const override {
            const auto& entry = GetIntersection(first, third);
            return GetIntersectionSize(second, entry);
        }

        size_t GetHeadSize(const ScaffoldVertex &vertex) const override {
            return (index_->GetHeadEntry(vertex)).size();
        }
        size_t GetTailSize(const ScaffoldVertex &vertex) const override {
            return (index_->GetTailEntry(vertex)).size();
        }

        SimpleVertexEntry GetIntersection(const SimpleVertexEntry &first,
                                          const SimpleVertexEntry &second) const override {
            SimpleVertexEntry result;
            std::set_intersection(first.begin(), first.end(), second.begin(), second.end(), std::inserter(result, result.end()));
            return result;
        }

        SimpleVertexEntry GetHeadEntry(const ScaffoldVertex &vertex) override {
            return index_->GetHeadEntry(vertex);
        }

        SimpleVertexEntry GetTailEntry(const ScaffoldVertex &vertex) override {
            return index_->GetTailEntry(vertex);
        }
    };

    typedef IntersectingScaffoldVertexIndexInfoExtractor<SimpleVertexEntry> SimpleIntersectingScaffoldVertexExtractor;
}
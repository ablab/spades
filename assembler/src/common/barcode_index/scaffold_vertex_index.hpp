//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "barcode_info_extractor.hpp"
#include "auxiliary_graphs/scaffold_graph/scaffold_vertex.hpp"
#include "assembly_graph/core/graph.hpp"
#include "adt/iterator_range.hpp"

namespace barcode_index {

template <class VertexEntryT>
class ScaffoldVertexIndexBuilder;

template <class VertexEntryT>
class ScaffoldVertexIndex {
    friend class ScaffoldVertexIndexBuilder<VertexEntryT>;
 public:
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef typename VertexEntryT::const_iterator const_iterator;
    typedef debruijn_graph::Graph Graph;

    ScaffoldVertexIndex(const Graph &g): g_(g), vertex_to_entry_() {}

    const VertexEntryT& GetHeadEntry(const ScaffoldVertex& vertex) const {
        return vertex_to_entry_.at(vertex);
    }
    const VertexEntryT& GetTailEntry(const ScaffoldVertex& vertex) const {
        return vertex_to_entry_.at(vertex.GetConjugateFromGraph(g_));
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
        return vertex_to_entry_.at(vertex.GetConjugateFromGraph(g_)).begin();
    }
    const_iterator GetTailEnd(const ScaffoldVertex& vertex) const {
        return vertex_to_entry_.at(vertex.GetConjugateFromGraph(g_)).end();
    }
    adt::iterator_range<const_iterator> GetTailRange(const ScaffoldVertex& vertex) const {
        return adt::make_range(GetTailBegin(vertex.GetConjugateFromGraph(g_)),
                               GetTailEnd(vertex.GetConjugateFromGraph(g_)));
    }

    bool Contains(const ScaffoldVertex &vertex) const {
        return vertex_to_entry_.find(vertex) != vertex_to_entry_.end();
    }
 private:
    void InsertEntry(const ScaffoldVertex& vertex, VertexEntryT&& entry) {
        vertex_to_entry_.insert({vertex, entry});
    }

    const Graph& g_;
    std::unordered_map<ScaffoldVertex, VertexEntryT> vertex_to_entry_;
};

typedef std::set<BarcodeId> SimpleVertexEntry;

typedef ScaffoldVertexIndex<SimpleVertexEntry> SimpleScaffoldVertexIndex;

class ScaffoldVertexIndexInfoExtractor {
 public:
    typedef typename scaffold_graph::ScaffoldVertex ScaffoldVertex;
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
  public:
    using Graph = debruijn_graph::Graph;

    BarcodeIndexInfoExtractorWrapper(const Graph &g, std::shared_ptr<FrameBarcodeIndexInfoExtractor> barcode_index_)
        : g_(g), barcode_extractor_(barcode_index_) {}

    size_t GetHeadSize(const ScaffoldVertex &vertex) const override {
        return barcode_extractor_->GetNumberOfBarcodes(vertex.GetFirstEdge());
    }
    size_t GetTailSize(const ScaffoldVertex &vertex) const override {
        return barcode_extractor_->GetNumberOfBarcodes(vertex.GetConjugateFromGraph(g_).GetFirstEdge());
    }
    size_t GetIntersectionSize(const ScaffoldVertex &first, const ScaffoldVertex &second) const override {
        return barcode_extractor_->GetNumberOfSharedBarcodes(first.GetLastEdge(), second.GetFirstEdge());
    }
    size_t GetIntersectionSize(const ScaffoldVertex &first,
                               const ScaffoldVertex &second,
                               const ScaffoldVertex &third) const override {
        return GetIntersectionSize(third, GetIntersection(first, second));
    }

    SimpleVertexEntry GetIntersection(const SimpleVertexEntry &first,
                                      const SimpleVertexEntry &second) const override {
        SimpleVertexEntry result;
        std::set_intersection(first.begin(), first.end(), second.begin(), second.end(),
                              std::inserter(result, result.end()));
        return result;
    }
    SimpleVertexEntry GetIntersection(const ScaffoldVertex &first, const ScaffoldVertex &second) const override {
        auto intersection = barcode_extractor_->GetSharedBarcodes(first.GetLastEdge(), second.GetFirstEdge());
        std::set<BarcodeId> result;
        std::copy(intersection.begin(), intersection.end(), std::inserter(result, result.begin()));
        return result;
    }
    size_t GetIntersectionSize(const ScaffoldVertex &middle, const SimpleVertexEntry &entry) const override {
        auto barcodes = barcode_extractor_->GetBarcodes(middle.GetFirstEdge());
        SimpleVertexEntry intersection;
        std::set_intersection(barcodes.begin(), barcodes.end(), entry.begin(), entry.end(),
                              std::inserter(intersection, intersection.begin()));
        return intersection.size();
    }
    SimpleVertexEntry GetHeadEntry(const ScaffoldVertex &/*vertex*/) override {
        VERIFY_MSG(false, "Head entry extractor from BarcodeIndexInfoExtractorWrapper is currently not supported");
        SimpleVertexEntry result;
        return result;
    }

    SimpleVertexEntry GetTailEntry(const ScaffoldVertex &vertex) override {
        return GetHeadEntry(vertex);
    }

  private:
    const Graph &g_;
    std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_;
};

class SimpleScaffoldVertexIndexInfoExtractor: public IntersectingScaffoldVertexIndexInfoExtractor<SimpleVertexEntry> {
  public:
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;

    explicit SimpleScaffoldVertexIndexInfoExtractor(std::shared_ptr<ScaffoldVertexIndex<SimpleVertexEntry>> index_)
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

  private:
    std::shared_ptr<ScaffoldVertexIndex<SimpleVertexEntry>> index_;
};

typedef IntersectingScaffoldVertexIndexInfoExtractor<SimpleVertexEntry> SimpleIntersectingScaffoldVertexExtractor;
}
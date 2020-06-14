//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <limits>
#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/core/action_handlers.hpp"
#include "assembly_graph/index/edge_info_updater.hpp"
#include "edge_index_refiller.hpp"


namespace io { namespace binary {
template<class Graph>
class EdgeIndexIO;
} }

namespace debruijn_graph {

/**
 * EdgeIndex is a structure to store info about location of certain k-mers in graph. It delegates all
 * container procedures to inner_index_ and all handling procedures to updater_.
 */
template<class Graph>
class EdgeIndex: public omnigraph::GraphActionHandler<Graph> {
    using InnerIndex32 = KmerFreeEdgeIndex<Graph, uint32_t>;
    using InnerIndex64 = KmerFreeEdgeIndex<Graph, uint64_t>;

    typedef typename Graph::EdgeId EdgeId;
public:
    typedef RtSeq KMer;
    static constexpr size_t NOT_FOUND = size_t(-1);

private:
    bool large_index_;
    void *inner_index_;

    EdgeInfoUpdater<Graph> updater_;
    EdgeIndexRefiller refiller_;

    template<class Index>
    std::pair<EdgeId, size_t> get(const Index *index, const KMer& kmer) const {
        auto kwh = index->ConstructKWH(kmer);
        if (index->contains(kwh)) {
            auto entry = index->get_value(kwh);
            return { entry.edge(), (size_t)entry.offset() };
        }
        
        return { EdgeId(), NOT_FOUND };
    }

    template<class Index>
    bool contains(const Index *index, const KMer& kmer) const {
        return index->contains(index->ConstructKWH(kmer));
    }

    template<class Index>
    void UpdateKmers(Index *index, EdgeId e) {
        updater_.UpdateKmers(this->g(), e, *index);
    }

    template<class Index>
    void DeleteKmers(Index *index, EdgeId e) {
        updater_.DeleteKmers(this->g(), e, *index);
    }

    template<class Index>
    void clear(Index *index) {
        if (!inner_index_)
            return;

        index->clear();
        delete index;
        // inner_index_ is always type-erased index here.
        inner_index_ = nullptr;
    }

    template<class Index>
    void Refill(Index *) {
        auto index = new Index(this->g());
        refiller_.Refill(*index, this->g());
        inner_index_ = index;
    }

    template<class Index>
    void Refill(Index *, const std::vector<EdgeId> &edges) {
        auto index = new Index(this->g());
        refiller_.Refill(*index, this->g(), edges);
        inner_index_ = index;
    }

    template<class Writer, class Index>
    void BinWrite(const Index *index, Writer &writer) const {
        index->BinWrite(writer);
    }

    template<class Reader, class Index>
    void BinRead(Index *, Reader &reader) {
        auto index = new Index(this->g());
        index->BinRead(reader);
        inner_index_ = index;
    }

public:
    EdgeIndex(const Graph& g, const std::string &workdir)
            : omnigraph::GraphActionHandler<Graph>(g, "EdgeIndex"),
              large_index_(true), inner_index_(nullptr),
              refiller_(workdir) {
        INFO("Size of edge index entries: "
             << sizeof(typename InnerIndex64::KmerPos) << "/"
             << sizeof(typename InnerIndex32::KmerPos));
    }

    virtual ~EdgeIndex() {
        clear();
        TRACE("~EdgeIndex OK")
    }

    size_t k() const {
        return this->g().k() + 1;
    }

#define DISPATCH_TO(method, ...)                                        \
    do {                                                                \
        if (large_index_) {                                             \
            return method(static_cast<InnerIndex64*>(inner_index_),##__VA_ARGS__); \
        } else {                                                        \
            return method(static_cast<InnerIndex32*>(inner_index_),##__VA_ARGS__); \
        }                                                               \
    } while(0)

    void HandleAdd(EdgeId e) override {
        DISPATCH_TO(UpdateKmers, e);
    }

    void HandleDelete(EdgeId e) override {
        DISPATCH_TO(DeleteKmers, e);
    }

    bool contains(const KMer& kmer) const {
        DISPATCH_TO(contains, kmer);
    }

    std::pair<EdgeId, size_t> get(const KMer& kmer) const {
        DISPATCH_TO(get, kmer);
    }

    void Refill() {
        clear();
        uint64_t max_id = this->g().max_eid();
        large_index_ = (max_id > std::numeric_limits<uint32_t>::max());
        INFO("Using " << (large_index_ ? "large" : "small") << " index (max_id = " << max_id << ")");
        DISPATCH_TO(Refill);

        INFO("Index refilled");
    }

    void Refill(const std::vector<EdgeId> &edges) {
        clear();

        uint64_t max_id = this->g().max_eid();
        large_index_ = (max_id > std::numeric_limits<uint32_t>::max());
        INFO("Using " << (large_index_ ? "large" : "small") << " index (max_id = " << max_id << ")");
        DISPATCH_TO(Refill, edges);

        INFO("Index refilled");
    }

    void clear() {
        DISPATCH_TO(clear);
    }

    static bool IsInvertable() {
        static_assert(InnerIndex32::storing_type::IsInvertable() == InnerIndex64::storing_type::IsInvertable(),
                      "Indices must be compatible");
        return InnerIndex32::storing_type::IsInvertable();
    }

    template<class Writer>
    void BinWrite(Writer &writer) const {
        writer << large_index_;
        DISPATCH_TO(BinWrite, writer);
    }

    template<class Reader>
    void BinRead(Reader &reader) {
        VERIFY(inner_index_ == nullptr);
        reader >> large_index_;
        DISPATCH_TO(BinRead, reader);
    }

};

#undef DISPATCH_TO

template<class Graph>
constexpr size_t EdgeIndex<Graph>::NOT_FOUND;

}

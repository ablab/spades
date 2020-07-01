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
 * container procedures to inner_index_ and all handling procedures to
 * renewer_ which is DataHashRenewer.
 */
template<class Graph>
class EdgeIndex: public omnigraph::GraphActionHandler<Graph> {
    using InnerIndex = KmerFreeEdgeIndex<Graph>;

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
    const std::pair<EdgeId, size_t> get(const KMer& kmer, const Index &index) const {
        auto kwh = index.ConstructKWH(kmer);
        if (!index.contains(kwh)) {
            return { EdgeId(), NOT_FOUND };
        } else {
            auto entry = index.get_value(kwh);
            return { entry.edge(), (size_t)entry.offset() };
        }
    }

    template<class Index>
    bool contains(const KMer& kmer, const Index &index) const {
        return index.contains(index.ConstructKWH(kmer));
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
        if (large_index_)
            delete static_cast<InnerIndex64*>(inner_index_);
        else
            delete static_cast<InnerIndex32*>(inner_index_);
        TRACE("~EdgeIndex OK")
    }

    size_t k() const {
        return this->g().k() + 1;
    }

    void HandleAdd(EdgeId e) override {
        if (large_index_)
            updater_.UpdateKmers(this->g(), e,
                                 *static_cast<InnerIndex64*>(inner_index_));
        else
            updater_.UpdateKmers(this->g(), e,
                                 *static_cast<InnerIndex32*>(inner_index_));
    }

    void HandleDelete(EdgeId e) override {
        if (large_index_)
            updater_.DeleteKmers(this->g(), e,
                                 *static_cast<InnerIndex64*>(inner_index_));
        else
            updater_.DeleteKmers(this->g(), e,
                                 *static_cast<InnerIndex32*>(inner_index_));
    }

    bool contains(const KMer& kmer) const {
        if (large_index_)
            return contains(kmer,
                            *static_cast<const InnerIndex64*>(inner_index_));
        else
            return contains(kmer,
                            *static_cast<const InnerIndex32*>(inner_index_));
    }

    const std::pair<EdgeId, size_t> get(const KMer& kmer) const {
        if (large_index_)
            return get(kmer,
                       *static_cast<const InnerIndex64*>(inner_index_));
        else
            return get(kmer,
                       *static_cast<const InnerIndex32*>(inner_index_));
    }

    void Refill() {
        clear();
        uint64_t max_id = this->g().max_eid();
        large_index_ = (max_id > std::numeric_limits<uint32_t>::max());
        if (large_index_) {
            INFO("Using large index (max_id = " << max_id << ")");
            auto index = new InnerIndex64(this->g());
            refiller_.Refill(*index, this->g());
            inner_index_ = index;
        } else {
            INFO("Using small index (max_id = " << max_id << ")");
            auto index = new InnerIndex32(this->g());
            refiller_.Refill(*index, this->g());
            inner_index_ = index;
        }

        INFO("Index refilled");
    }

    void Update() {
        updater_.UpdateAll();
    }

    void clear() {
        if (!inner_index_)
            return;
        
        if (large_index_) {
            auto index = static_cast<InnerIndex64*>(inner_index_);
            index->clear();
            delete index;
        } else {
            auto index = static_cast<InnerIndex32*>(inner_index_);
            index->clear();
            delete index;
        }
        inner_index_ = nullptr;
    }

    static bool IsInvertable() {
        return InnerIndex::storing_type::IsInvertable();
    }

    template<class Writer>
    void BinWrite(Writer &writer) const {
        writer << large_index_;
        if (large_index_)
            static_cast<const InnerIndex64*>(inner_index_)->BinWrite(writer);
        else
            static_cast<const InnerIndex32*>(inner_index_)->BinWrite(writer);
    }

    template<class Reader>
    void BinRead(Reader &reader) {
        VERIFY(inner_index_ == nullptr);
        reader >> large_index_;
        if (large_index_) {
            inner_index_ = new InnerIndex64(this->g());
            static_cast<InnerIndex64*>(inner_index_)->BinRead(reader);
        } else {
            inner_index_ = new InnerIndex32(this->g());
            static_cast<InnerIndex32*>(inner_index_)->BinRead(reader);
        }
    }

};

template<class Graph>
constexpr size_t EdgeIndex<Graph>::NOT_FOUND;

}

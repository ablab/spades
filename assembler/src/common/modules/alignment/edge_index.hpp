//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

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

    typedef typename Graph::EdgeId EdgeId;
public:
    typedef RtSeq KMer;
    static constexpr size_t NOT_FOUND = size_t(-1);

private:
    InnerIndex inner_index_;
    EdgeInfoUpdater<Graph> updater_;
    EdgeIndexRefiller refiller_;
    
public:
    EdgeIndex(const Graph& g, const std::string &workdir)
            : omnigraph::GraphActionHandler<Graph>(g, "EdgeIndex"),
              inner_index_(g),
              refiller_(workdir) {
        INFO("Size of edge index entry: " << sizeof(EdgeInfo<EdgeId>));
    }

    virtual ~EdgeIndex() {
        TRACE("~EdgeIndex OK")
    }

    size_t k() const {
        return inner_index_.k();
    }

    void HandleAdd(EdgeId e) override {
        updater_.UpdateKmers(this->g(), e, inner_index_);
    }

    void HandleDelete(EdgeId e) override {
        updater_.DeleteKmers(this->g(), e, inner_index_);
    }

    bool contains(const KMer& kmer) const {
        VERIFY(this->IsAttached());
        return inner_index_.contains(inner_index_.ConstructKWH(kmer));
    }

    const std::pair<EdgeId, size_t> get(const KMer& kmer) const {
        VERIFY(this->IsAttached());
        auto kwh = inner_index_.ConstructKWH(kmer);
        if (!inner_index_.contains(kwh)) {
            return { EdgeId(), NOT_FOUND };
        } else {
            EdgeInfo<EdgeId> entry = inner_index_.get_value(kwh);
            return { entry.edge(), (size_t)entry.offset() };
        }
    }

    void Refill() {
        clear();
        refiller_.Refill(inner_index_, this->g());
        INFO("Index refilled");
    }

    void Update() {
        updater_.UpdateAll();
    }

    void clear() {
        inner_index_.clear();
    }

    static bool IsInvertable() {
        return InnerIndex::storing_type::IsInvertable();
    }

    template<class Writer>
    void BinWrite(Writer &writer) const {
        inner_index_.BinWrite(writer);
    }

    template<class Reader>
    void BinRead(Reader &reader) {
        inner_index_.BinRead(reader);
    }

};

template<class Graph>
constexpr size_t EdgeIndex<Graph>::NOT_FOUND;

}

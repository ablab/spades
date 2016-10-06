//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "common/assembly_graph/core/graph.hpp"
#include "common/assembly_graph/core/action_handlers.hpp"
#include "utils/indices/edge_info_updater.hpp"
#include "edge_index_refiller.hpp"
    
namespace debruijn_graph {

/**
 * EdgeIndex is a structure to store info about location of certain k-mers in graph. It delegates all
 * container procedures to inner_index_ and all handling procedures to
 * renewer_ which is DataHashRenewer.
 */
template<class Graph>
class EdgeIndex: public omnigraph::GraphActionHandler<Graph> {

public:
    typedef typename Graph::EdgeId EdgeId;
    using InnerIndex = KmerFreeEdgeIndex<Graph, DefaultStoring>;
    typedef Graph GraphT;
    typedef typename InnerIndex::KMer KMer;
    typedef typename InnerIndex::KMerIdx KMerIdx;
    typedef typename InnerIndex::KmerPos Value;

private:
    InnerIndex inner_index_;
    EdgeInfoUpdater<InnerIndex, Graph> updater_;
    EdgeIndexRefiller refiller_;
    bool delete_index_;

public:
    EdgeIndex(const Graph& g, const std::string &workdir)
            : omnigraph::GraphActionHandler<Graph>(g, "EdgeIndex"),
              inner_index_(g, workdir),
              updater_(g, inner_index_),
              delete_index_(true) {
    }

    virtual ~EdgeIndex() {
        TRACE("~EdgeIndex OK")
    }

    InnerIndex &inner_index() {
        return inner_index_;
    }

    size_t k() const {
        return inner_index_.k();
    }

    const InnerIndex &inner_index() const {
        VERIFY(this->IsAttached());
        return inner_index_;
    }

    void HandleAdd(EdgeId e) override {
        updater_.UpdateKmers(e);
    }

    void HandleDelete(EdgeId e) override {
        updater_.DeleteKmers(e);
    }

    bool contains(const KMer& kmer) const {
        VERIFY(this->IsAttached());
        return inner_index_.contains(inner_index_.ConstructKWH(kmer));
    }

    const pair<EdgeId, size_t> get(const KMer& kmer) const {
        VERIFY(this->IsAttached());
        auto kwh = inner_index_.ConstructKWH(kmer);
        if (!inner_index_.contains(kwh)) {
            return make_pair(EdgeId(0), -1u);
        } else {
            EdgeInfo<EdgeId> entry = inner_index_.get_value(kwh);
            return std::make_pair(entry.edge_id, (size_t)entry.offset);
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

};
}

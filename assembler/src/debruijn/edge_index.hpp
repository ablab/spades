//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "openmp_wrapper.h"

#include "omni/omni_utils.hpp"
#include "adt/kmer_map.hpp"

#include "debruijn_graph.hpp"
#include "standard.hpp"
#include "indices/edge_index_builders.hpp"

namespace debruijn_graph {

/**
 * EdgeIndex is a structure to store info about location of certain k-mers in graph. It delegates all
 * container procedures to inner_index_ and all handling procedures to
 * renewer_ which is DataHashRenewer.
 * @see DeBruijnKMerIndex
 * @see DataHashRenewer
 */
//fixme template params
template<class Graph, class Seq /*= runtime_k::RtSeq*/,
        class Index /*= KmerFreeEdgeIndex<Graph, Seq>*/>
class EdgeIndex: public GraphActionHandler<Graph> {

public:
    typedef typename Graph::EdgeId EdgeId;
    typedef Index InnerIndexT;
    typedef Graph GraphT;
    typedef typename Index::KMer KMer;
    typedef typename Index::KMerIdx KMerIdx;
    typedef typename Index::Value Value;

private:
    Index inner_index_;
    EdgeInfoUpdater<Index, Graph> updater_;
    bool delete_index_;

public:

    EdgeIndex(const Graph& g, const std::string &workdir)
            : GraphActionHandler<Graph>(g, "EdgeIndex"),
              inner_index_(g, workdir),
              updater_(g, inner_index_),
              delete_index_(true) {
    }

    virtual ~EdgeIndex() {
        TRACE("~EdgeIndex OK")
    }

    Index &inner_index() {
        return inner_index_;
    }

    size_t k() const {
        return inner_index_.k();
    }

    const Index &inner_index() const {
        VERIFY(this->IsAttached());
        return inner_index_;
    }

    virtual void HandleAdd(EdgeId e) {
        updater_.UpdateKmers(e);
    }

    virtual void HandleDelete(EdgeId e) {
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
        typedef typename EdgeIndexHelper<InnerIndexT>::GraphPositionFillingIndexBuilderT IndexBuilder;
        //also makes an update!
        //todo pass appropriate 3-rd arg
        IndexBuilder().BuildIndexFromGraph(inner_index_, this->g());
    }

    void Update() {
        updater_.UpdateAll();
    }

    void clear() {
        inner_index_.clear();
    }

};
}

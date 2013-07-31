//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
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
 * container procedures to inner_index_ which is DeBruijnKMerIndex and all handling procedures to
 * renewer_ which is DataHashRenewer.
 * @see DeBruijnKMerIndex
 * @see DataHashRenewer
 */
//fixme template params
template<class Graph, class Seq /*= runtime_k::RtSeq*/,
        class Index /*= DeBruijnEdgeIndex<KmerFreeDeBruijnEdgeIndex<Graph, Seq>>*/>
class EdgeIndex: public GraphActionHandler<Graph> {

public:
    typedef typename Graph::EdgeId EdgeId;
    typedef Index InnerIndexT;
    typedef Graph GraphT;
    typedef typename Index::KMer KMer;
    typedef typename Index::KMerIdx KMerIdx;

private:
    Index inner_index_;
    EdgeInfoUpdater<Index, Graph> updater_;
    bool delete_index_;

public:

    EdgeIndex(const Graph& g, size_t k, const std::string &workdir)
            : GraphActionHandler<Graph>(g, "EdgeIndex"),
              inner_index_((unsigned) k, g, workdir),
              updater_(g, inner_index_),
              delete_index_(true) {
    }

    virtual ~EdgeIndex() {
        TRACE("~EdgeIndex OK")
    }

    Index &inner_index() {
        return inner_index_;
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
        return inner_index_.contains(kmer);
    }

    const pair<EdgeId, size_t> get(const KMer& kmer) const {
        VERIFY(this->IsAttached());
        KMerIdx idx = inner_index_.seq_idx(kmer);
        if (!inner_index_.contains(idx, kmer)) {
            return make_pair(EdgeId(0), -1u);
        } else {
            return inner_index_.get(idx, kmer);
        }
    }

//  KMerIdx seq_idx(const Kmer& kmer) const {
//    VERIFY(this->IsAttached());
//    return inner_index_.seq_idx(kmer);
//  }
//  bool contains(const KMerIdx idx) const {
//    VERIFY(this->IsAttached());
//    return inner_index_.contains(idx);
//  }
//  bool contains(const KMerIdx idx, const Kmer& kmer) const {
//  return inner_index_.contains(idx, kmer);
//  }
//  const pair<EdgeId, size_t> get(KMerIdx idx) const {
//    VERIFY(this->IsAttached());
//    return inner_index_.get(idx);
//  }

    void Refill() {
        clear();
        typedef typename EdgeIndexHelper<InnerIndexT>::GraphPositionFillingIndexBuilderT IndexBuilder;
        IndexBuilder().BuildIndexFromGraph(inner_index_, this->g());
        //Update();
    }

    void Update() {
        updater_.UpdateAll();
    }

    void clear() {
        inner_index_.clear();
    }

};}

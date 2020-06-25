//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "edge_position_index.hpp"

#include "assembly_graph/core/graph_iterators.hpp"
#include "sequence/sequence.hpp"

#include "utils/parallel/openmp_wrapper.h"

namespace debruijn_graph {

template<typename Index, typename Graph>
class EdgeInfoUpdater {
    typedef typename Index::KMer Kmer;
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Index::KeyWithHash KeyWithHash;

    const Graph &g_;
    Index &index_;

    bool DeleteIfEqual(const KeyWithHash& kwh, EdgeId e) {
        if (!index_.contains(kwh))
            return false;

        if (index_.get_value(kwh).edge() == e) {
            index_.get_raw_value_reference(kwh).clear();
            return true;
        }
        return false;
    }

    void UpdateKMers(const Sequence &nucls, EdgeId e) {
        VERIFY(nucls.size() >= index_.k());
        KeyWithHash kwh = index_.ConstructKWH(Kmer(index_.k(), nucls));
        if (kwh.is_minimal())
            index_.PutInIndex(kwh, e, 0);
        for (size_t i = index_.k(), n = nucls.size(); i < n; ++i) {
            kwh <<= nucls[i];
            if (kwh.is_minimal())
                index_.PutInIndex(kwh, e, i - index_.k() + 1);
        }
    }

    void DeleteKMers(const Sequence &nucls, EdgeId e) {
        VERIFY(nucls.size() >= index_.k());
        KeyWithHash kwh = index_.ConstructKWH(Kmer(index_.k(), nucls));
        DeleteIfEqual(kwh, e);
        for (size_t i = index_.k(), n = nucls.size(); i < n; ++i) {
            kwh <<= nucls[i];
            DeleteIfEqual(kwh, e);
        }
    }

 public:
    /**
     * Creates EdgeInfoUpdater for specified graph and index
     * @param g graph to be indexed
     * @param index index to be synchronized with graph
     */
    EdgeInfoUpdater(const Graph& g, Index& index)
            : g_(g), index_(index) { }

    void UpdateKmers(EdgeId e) {
        UpdateKMers(g_.EdgeNucls(e), e);
    }

    void DeleteKmers(EdgeId e) {
        DeleteKMers(g_.EdgeNucls(e), e);
    }

    void UpdateAll() {
        unsigned nthreads = omp_get_max_threads();

        omnigraph::IterationHelper<Graph, EdgeId> edges(g_);
        auto iters = edges.Chunks(16 * nthreads);

        #pragma omp parallel for schedule(guided)
        for (size_t i = 0; i < iters.size() - 1; ++i) {
            TRACE("Processing chunk #" << i);
            for (auto it = iters[i]; it != iters[i + 1]; ++it) {
                UpdateKmers(*it);
            }
        }
    }

 private:
    DECL_LOGGER("EdgeInfoUpdater")
};

}

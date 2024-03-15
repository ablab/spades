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

template<typename Graph>
class EdgeInfoUpdater {
    typedef typename Graph::EdgeId EdgeId;

    template<class Index>
    bool DeleteIfEqual(const typename Index::KeyWithHash& kwh, EdgeId e, Index &index) {
        if (!index.contains(kwh))
            return false;

        if (index.get_value(kwh).edge() == e) {
            index.get_raw_value_reference(kwh).clear();
            return true;
        }
        return false;
    }

    template<class Index>
    void UpdateKMers(const Sequence &nucls, EdgeId e, Index &index) {
        VERIFY(nucls.size() >= index.k());
        typename Index::KeyWithHash kwh = index.ConstructKWH(typename Index::KMer(index.k(), nucls));
        if (kwh.is_minimal())
            index.PutInIndex(kwh, e, 0);
        for (size_t i = index.k(), n = nucls.size(); i < n; ++i) {
            kwh <<= nucls[i];
            if (kwh.is_minimal())
                index.PutInIndex(kwh, e, i - index.k() + 1);
        }
    }

    template<class Index>
    void DeleteKMers(const Sequence &nucls, EdgeId e, Index &index) {
        VERIFY(nucls.size() >= index.k());
        typename Index::KeyWithHash kwh = index.ConstructKWH(typename Index::KMer(index.k(), nucls));
        DeleteIfEqual(kwh, e, index);
        for (size_t i = index.k(), n = nucls.size(); i < n; ++i) {
            kwh <<= nucls[i];
            DeleteIfEqual(kwh, e, index);
        }
    }

 public:
    template<class Index>
    void UpdateKmers(const Graph &g, EdgeId e, Index &index) {
        UpdateKMers(g.EdgeNucls(e), e, index);
    }

    template<class Index>
    void DeleteKmers(const Graph &g, EdgeId e, Index &index) {
        DeleteKMers(g.EdgeNucls(e), e, index);
    }

    template<class Index>
    void UpdateAll(const Graph &g, Index &index) {
        unsigned nthreads = omp_get_max_threads();

        omnigraph::IterationHelper<Graph, EdgeId> edges(g);
        auto iters = edges.Chunks(16 * nthreads);

        #pragma omp parallel for schedule(guided)
        for (size_t i = 0; i < iters.size() - 1; ++i) {
            TRACE("Processing chunk #" << i);
            for (auto it = iters[i]; it != iters[i + 1]; ++it) {
                UpdateKmers(g, *it, index);
            }
        }
    }

    template<class Index>
    void Update(const Graph &g, Index &index, const std::vector<EdgeId> &edges) {
#pragma omp parallel for schedule(guided)
        for (size_t i = 0; i < edges.size(); ++i) {
            UpdateKmers(g, edges[i], index);
        }
    }

 private:
    DECL_LOGGER("EdgeInfoUpdater")
};

}

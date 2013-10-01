#pragma once

#include "standard.hpp"

namespace debruijn_graph {

template<typename Index, typename Graph>
class EdgeInfoUpdater {
    typedef typename Index::KMer Kmer;
//    typedef typename Index::GraphT Graph;
    typedef typename Graph::EdgeId EdgeId;

    const Graph &g_;
    Index &index_;

    void UpdateKMers(const Sequence &nucls, EdgeId e) {
        VERIFY(nucls.size() >= index_.k());
        Kmer kmer(index_.k(), nucls);

        index_.PutInIndex(kmer, e, 0);
        for (size_t i = index_.k(), n = nucls.size(); i < n; ++i) {
            kmer <<= nucls[i];
            index_.PutInIndex(kmer, e, i - index_.k() + 1);
        }
    }

    void DeleteKMers(const Sequence &nucls, EdgeId e) {
        VERIFY(nucls.size() >= index_.k());
        Kmer kmer(index_.k(), nucls);
        index_.DeleteIfEqual(kmer, e);
        for (size_t i = index_.k(), n = nucls.size(); i < n; ++i) {
            kmer <<= nucls[i];
            index_.DeleteIfEqual(kmer, e);
        }
    }

 public:
    /**
     * Creates DataHashRenewer for specified graph and index
     * @param g graph to be indexed
     * @param index index to be synchronized with graph
     */
    EdgeInfoUpdater(const Graph& g, Index& index)
            : g_(g),
              index_(index) {
    }

    void UpdateKmers(EdgeId e) {
        Sequence nucls = g_.EdgeNucls(e);
        UpdateKMers(nucls, e);
    }

    void DeleteKmers(EdgeId e) {
        Sequence nucls = g_.EdgeNucls(e);
        DeleteKMers(nucls, e);
    }

    void UpdateAll() {
        for (auto it = g_.ConstEdgeBegin(); !it.IsEnd(); ++it) {
          UpdateKmers(*it);
        }
    }

 private:
    DECL_LOGGER("EdgeInfoUpdater")
};

}

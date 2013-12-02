#pragma once

#include "standard.hpp"

namespace debruijn_graph {

template<typename Index, typename Graph>
class EdgeInfoUpdater {
    typedef typename Index::KMer Kmer;
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Index::KeyWithHash KeyWithHash;
    typedef typename Index::Value EdgeInfo;

    const Graph &g_;
    Index &index_;


    void PutInIndex(const KeyWithHash &kwh, EdgeId id, size_t offset) {
        if (index_.valid(kwh)) {
            auto &entry = index_.get_raw_value_reference(kwh);
            if (!entry.valid() || index_.contains(kwh)) {
                index_.put_value(kwh, EdgeInfo(id, (unsigned)offset, entry.count));
            }
        }
    }

  	//todo why do we need to check equality???!!!
  	bool DeleteIfEqual(const KeyWithHash &kwh, EdgeId e) {
  		if (!index_.contains(kwh))
  			return false;
  		if (index_.get_value(kwh).edge_id == e) {
  		    index_.get_raw_value_reference(kwh).invalidate();
  			return true;
  		}
  		return false;
  	}

    void UpdateKMers(const Sequence &nucls, EdgeId e) {
        VERIFY(nucls.size() >= index_.k());
        KeyWithHash kwh = index_.ConstructKWH(Kmer(index_.k(), nucls));
        index_.PutInIndex(kwh, e, 0);
        for (size_t i = index_.k(), n = nucls.size(); i < n; ++i) {
        	kwh <<= nucls[i];
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

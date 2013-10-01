#pragma once

#include "indices/debruijn_kmer_index.hpp"
#include "omni/coverage.hpp"
#include "graph_pack.hpp"
#include "verify.hpp"
#include <vector>
#include <map>
#include <set>
#include <string>
#include <iostream>
#include <fstream>

namespace debruijn_graph {

template<class Graph, class Index>
class FlankingCoverage : public GraphActionHandler<Graph>, public omnigraph::AbstractFlankingCoverage<Graph> {
  typedef GraphActionHandler<Graph> base;
  typedef typename Graph::EdgeId EdgeId;
  typedef typename Graph::VertexId VertexId;
  /*
   Iterates over kmer index saves values of coverage on the ends of edges
   */

  std::map<EdgeId, double> in_coverage_;
  std::map<EdgeId, double> out_coverage_;
  const Index& kmer_index_;
  const size_t averaging_range_;

  double CountAvgCoverage(EdgeId e, size_t offset) const {
    size_t k = this->g().k();
    VERIFY(offset == 0 || offset + averaging_range_ == this->g().length(e));
    unsigned size_bound = (unsigned) std::min(averaging_range_, this->g().length(e));
    const Sequence& seq = this->g().EdgeNucls(e);

    size_t edge_coverage_in = 0;

    runtime_k::RtSeq kpomer(k + 1, seq, offset);
    kpomer >>= 0;
    for (size_t i = 0; i < size_bound; ++i) {
      kpomer <<= seq[offset + i + k];
//      VERIFY(kmer_index_.contains(kpomer));
      if(kmer_index_.contains(kpomer))
    	  edge_coverage_in += kmer_index_[kpomer].count;
      else {
//    	  cout << "warn! kmer not in index" << endl;
      }
    }

    return double(edge_coverage_in) / size_bound;
  }

  double CountInCoverage(EdgeId e) const {
    return CountAvgCoverage(e, 0);
  }

  double CountOutCoverage(EdgeId e) const {
    return CountAvgCoverage(e, std::max(this->g().length(e), averaging_range_) - averaging_range_);
  }

 public:

  FlankingCoverage(const Graph& g, const Index& kmer_index,
                   unsigned averaging_range)
      : base(g, "FlankingCoverage"),
        kmer_index_(kmer_index),
        averaging_range_(averaging_range) {
    for (auto it = g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
      EdgeId e = *it;
      in_coverage_.insert(std::make_pair(*it, CountInCoverage(e)));
      out_coverage_.insert(std::make_pair(*it, CountOutCoverage(e)));
    }
  }

  //todo rename
  double GetInCov(EdgeId edge) const {
    return get(in_coverage_, edge);
  }

  //todo rename
  double GetOutCov(EdgeId edge) const {
    return get(out_coverage_, edge);
  }

  /*virtual */
  void HandleAdd(EdgeId e) {
    in_coverage_.insert(std::make_pair(e, CountInCoverage(e)));
    out_coverage_.insert(std::make_pair(e, CountOutCoverage(e)));
  }

  /*virtual*/
  void HandleDelete(EdgeId e) {
    in_coverage_.erase(e);
    out_coverage_.erase(e);
  }

  double LocalCoverage(EdgeId e, VertexId v) const {
    if (this->g().EdgeStart(e) == v) {
      return GetInCov(e);
    } else if (this->g().EdgeEnd(e) == v) {
      return GetOutCov(e);
    } else {
      VERIFY(false);
      return 0.0;
    }
  }

 private:
  DECL_LOGGER("FlankingCoverage");

};

}

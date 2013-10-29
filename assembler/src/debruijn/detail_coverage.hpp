#pragma once

#include "verify.hpp"
#include "omni/coverage.hpp"
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

template<class Graph>
class NewFlankingCoverage : public GraphActionHandler<Graph>,
        public omnigraph::AbstractFlankingCoverage<Graph> {
    typedef GraphActionHandler<Graph> base;
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef pair<EdgeId, unsigned> Pos;

    Graph& g_;
    const size_t averaging_range_;

    template<class CoverageIndex>
    unsigned Count(const runtime_k::RtSeq& kpomer,
                   const CoverageIndex& index) const {
        //todo optimize
        if (index.contains(kpomer)) {
            return index[kpomer].count;
        } else {
            return 0;
        }
    }

    template<class CoverageIndex>
    pair<Pos, unsigned> Next(const Pos& curr, const CoverageIndex& index) const {
        const Graph& g = this->g();
        if (curr.second == g.length(curr.first) - 1) {
            EdgeId best(0);
            unsigned best_cnt = 0;
            FOREACH(EdgeId e, g.OutgoingEdges(g.EdgeEnd(curr.first))) {
                unsigned cnt = Count(runtime_k::RtSeq(g.k() + 1, g.EdgeNucls(e)), index);
                if (cnt > best_cnt) {
                    best_cnt = cnt;
                    best = e;
                }
            }
            if (best_cnt > 0) {
                return make_pair(make_pair(best, 0), best_cnt);
            } else {
                return make_pair(make_pair(EdgeId(0), -1u), 0);
            }
        } else {
            return make_pair(
                    make_pair(curr.first, curr.second + 1),
                    Count(runtime_k::RtSeq(g.k() + 1, g.EdgeNucls(curr.first), curr.second + 1),
                          index));
        }
    }

    double CountInCoverage(EdgeId e) const {
        return CountAvgCoverage(e, 0);
    }

    void SetRawCoverage(EdgeId e, unsigned cov) {
        g_.data(e).set_flanking_coverage(cov);
    }

    unsigned RawCoverage(EdgeId e) const {
        return g_.data(e).flanking_coverage();
    }

    unsigned InterpolateCoverage(EdgeId e, size_t l) const {
        VERIFY(l <= averaging_range_);
        VERIFY(l < g_.length(e));
        return unsigned(math::round(double(RawCoverage(e)) / double(averaging_range_) * double(l)));
    }

public:

    //todo think about interactions with gap closer
    NewFlankingCoverage(Graph& g, size_t averaging_range)
            : base(g, "NewFlankingCoverage"), g_(g),
              averaging_range_(averaging_range) {
    }

    size_t averaging_range() const {
        return averaging_range_;
    }

    //todo currently left for saves compatibility! remove later!
    template<class CoverageIndex>
    void Fill(const CoverageIndex& count_index) {
        TRACE("Filling flanking coverage from index");

        for (auto I = count_index.value_cbegin(), E = count_index.value_cend();
                I != E; ++I) {
            const auto& edge_info = *I;
            EdgeId e = edge_info.edge_id;
            unsigned offset = edge_info.offset;
            unsigned count = edge_info.count;
            VERIFY(offset != -1u);
            VERIFY(e.get() != NULL);
            if (offset < averaging_range_) {
                IncRawCoverage(e, count);
            }
        }
    }

    void IncRawCoverage(EdgeId e, unsigned count) {
        g_.data(e).inc_flanking_coverage(count);
    }

    double CoverageOfStart(EdgeId e) const {
        size_t averaging = std::min(this->g().length(e), averaging_range_);
        return double(RawCoverage(e)) / double(averaging);
    }

//    double operator[](EdgeId e) const {
//        return CoverageOfStart(e);
//    }

    double CoverageOfEnd(EdgeId e) const {
        return CoverageOfStart(this->g().conjugate(e));
    }

    virtual void HandleAdd(EdgeId /*e*/) {
    }

    virtual void HandleMerge(const vector<EdgeId>& old_edges, EdgeId new_edge) {
//        SetRawCoverage(new_edge, RawCoverage(old_edges.front()));
        size_t kpomers_left = averaging_range_;
        unsigned acc = 0;
        FOREACH(EdgeId e, old_edges) {
            if (kpomers_left >= g_.length(e)) {
                acc += RawCoverage(e);
                kpomers_left -= g_.length(e);
            } else {
                if (kpomers_left != 0)
                    acc += InterpolateCoverage(e, kpomers_left);
                break;
            }
        }
        SetRawCoverage(new_edge, acc);
    }

    virtual void HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) {
        SetRawCoverage(new_edge, RawCoverage(edge1) + RawCoverage(edge2));
    }

    virtual void HandleSplit(EdgeId old_edge, EdgeId new_edge_1,
                             EdgeId /*new_edge_2*/) {
        //hard to think of any other solution
        SetRawCoverage(new_edge_1, RawCoverage(old_edge));
    }

    virtual void HandleDelete(EdgeId e) {
        SetRawCoverage(e, 0);
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

    //left for compatibility
    //todo rename
    double GetInCov(EdgeId e) const {
        return CoverageOfStart(e);
    }

    //todo rename
    double GetOutCov(EdgeId e) const {
        return CoverageOfEnd(e);
    }

    //////////////////////////

    void Save(EdgeId e, ostream& out) const {
        out << RawCoverage(e);
    }

    void Load(EdgeId e, istream& in) {
        unsigned cov;
        in >> cov;
        SetRawCoverage(e, cov);
    }

private:
    DECL_LOGGER("NewFlankingCoverage")
    ;
};


template<class Graph, class CountIndex>
class SimultaneousCoverageFiller {
    const Graph& g_;
    const CountIndex& count_index_;
    NewFlankingCoverage<Graph>& flanking_coverage_;
    CoverageIndex<Graph>& coverage_index_;

public:
    SimultaneousCoverageFiller(const Graph& g, const CountIndex& count_index,
                               NewFlankingCoverage<Graph>& flanking_coverage,
                               CoverageIndex<Graph>& coverage_index) :
                                   g_(g),
                                   count_index_(count_index),
                                   flanking_coverage_(flanking_coverage),
                                   coverage_index_(coverage_index) {
    }

    void Fill() {
        for (auto I = count_index_.value_cbegin(), E = count_index_.value_cend();
                I != E; ++I) {
            const auto& edge_info = *I;
            EdgeId e = edge_info.edge_id;
            unsigned offset = edge_info.offset;
            unsigned count = edge_info.count;
            VERIFY(offset != -1u);
            VERIFY(e.get() != NULL);
            coverage_index_.IncRawCoverage(e, count);
            if (offset < flanking_coverage_.averaging_range()) {
                flanking_coverage_.IncRawCoverage(e, count);
            }
        }
    }

};

template<class Graph, class CountIndex>
void FillCoverageAndFlanking(const CountIndex& count_index, Graph& g,
                             NewFlankingCoverage<Graph>& flanking_coverage) {
    SimultaneousCoverageFiller<Graph, CountIndex> filler(g, count_index, flanking_coverage, g.coverage_index());
    filler.Fill();
}

}

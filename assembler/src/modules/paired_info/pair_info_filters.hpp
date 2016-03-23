//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef PAIR_INFO_FILTERS_HPP_
#define PAIR_INFO_FILTERS_HPP_

#include "paired_info_helpers.hpp"

namespace omnigraph {

namespace de {

template<class Graph>
class AbstractPairInfoChecker{
private:
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    typedef PairInfo<EdgeId> PairInfoT;

protected:
    const Graph& graph_;

public:
    AbstractPairInfoChecker(const Graph &graph) : graph_(graph) { }

    virtual bool Check(const PairInfoT&) {
        return true;
    }

    virtual bool Check(EdgeId, EdgeId) {
        return true;
    }

    virtual ~AbstractPairInfoChecker() {    }
};

template<class Graph>
class PairInfoWeightChecker : public AbstractPairInfoChecker<Graph>{
 private:
  typedef typename Graph::EdgeId EdgeId;
  typedef PairInfo<EdgeId> PairInfoT;
  double weight_threshold_;

 public:
  PairInfoWeightChecker(const Graph& graph, double weight_threshold) :
    AbstractPairInfoChecker<Graph>(graph), weight_threshold_(weight_threshold) {
  }

  bool Check(const PairInfoT& info) {
    return math::ge(info.weight(), weight_threshold_);
  }
};

template<class Graph>
class PairInfoWeightCheckerWithCoverage: public AbstractPairInfoChecker<Graph> {
 private:
  typedef typename Graph::EdgeId EdgeId;
  typedef PairInfo<EdgeId> PairInfoT;
  double weight_threshold_;

 public:
  PairInfoWeightCheckerWithCoverage(const Graph& graph, double weight_threshold) :
    AbstractPairInfoChecker<Graph>(graph), weight_threshold_(weight_threshold){
  }

  bool Check(const PairInfoT& info) {
    double info_weight = info.weight();
    return math::ge(info_weight, weight_threshold_)
       || (math::ge(info_weight, 0.1 * this->graph_.coverage(info.first)))
       || (math::ge(info_weight, 0.1 * this->graph_.coverage(info.second)));
  }
};

template <class Graph>
class AmbiguousPairInfoChecker : public AbstractPairInfoChecker<Graph> {

  typedef typename Graph::EdgeId EdgeId;
  typedef typename Graph::VertexId VertexId;
  typedef PairInfo<EdgeId> PairInfoT;
  typedef boost::optional<EdgeId> OptEdgeId;

  AbstractPairInfoChecker<Graph> &standard_filter_;
  const PairedInfoIndexT<Graph>& index_;

  double haplom_threshold_;
  double relative_length_threshold_;
  double relative_seq_threshold_;

  bool IsEdgeOneHaplome(EdgeId edge){
      return this->graph_.coverage(edge) < 1.5 * haplom_threshold_;
  }

  bool IsPairInfoGood(EdgeId edge1, EdgeId edge2){
      return index_.Get(edge1, edge2).size() <= 1;
  }

  bool EdgesAreFromSimpleBulgeWithAmbPI(const PairInfoT& info){
      EdgeId edge1 = info.first;
      EdgeId edge2 = info.second;
      // edge is auto reverse complementary
      TRACE("Check for auto reverse complementary");
      if(this->graph_.conjugate(edge1) == info.second)
          return false;
      TRACE("Done");

      TRACE("Check for coverage 1x haplome for edge from pair info");
      if(!IsEdgeOneHaplome(edge1) || !IsEdgeOneHaplome(edge2))
          return false;
      TRACE("Done");

      // first edge is not side of simple bulge
      TRACE("Check for bulge side for the 1st edge");
      OptEdgeId edge1_alt = GetOtherSideOfSimpleBulge(edge1);
      if(!edge1_alt.is_initialized())
          return false;
      TRACE("Done");

      // second edge is not side of simple bulge
      TRACE("Check for bulge side for the 2nd edge");
      OptEdgeId edge2_alt = GetOtherSideOfSimpleBulge(edge2);
      if(!edge2_alt.is_initialized())
          return false;
      TRACE("Done");

      TRACE("Check for coverage 1x haplome for edge from alternative bulge sides");
      if(!IsEdgeOneHaplome(edge1_alt.get()) || !IsEdgeOneHaplome(edge2_alt.get()))
          return false;
      TRACE("Done");

      TRACE("Check for multiplicity of pair info");
      if(!(IsPairInfoGood(edge1, edge2_alt.get()) &&
              IsPairInfoGood(edge1_alt.get(), edge2) &&
              IsPairInfoGood(edge1_alt.get(), edge2_alt.get())))
          return false;
      TRACE("Done");

      return true;
  }

  DEWeight GetPairInfoWeight(EdgeId edge1, EdgeId edge2){
      auto hist = index_.Get(edge1, edge2);
      return (hist.size() == 1) ? hist.begin()->weight : DEWeight(0);
  }

  bool InnerCheck(const PairInfoT& info){

      EdgeId edge1 = info.first;
      EdgeId edge2 = info.second;

      // get second edges of simple bulge
      OptEdgeId opt_edge1_alt = GetOtherSideOfSimpleBulge(edge1);
      VERIFY(opt_edge1_alt.is_initialized());
      EdgeId edge1_alt = opt_edge1_alt.get();

      OptEdgeId opt_edge2_alt = GetOtherSideOfSimpleBulge(edge2);
      VERIFY(opt_edge2_alt.is_initialized());
      EdgeId edge2_alt = opt_edge2_alt.get();

      double direct_weight = GetPairInfoWeight(edge1, edge2) +
              GetPairInfoWeight(edge1_alt, edge2_alt);

      double reverse_weight = GetPairInfoWeight(edge1, edge2_alt) +
              GetPairInfoWeight(edge1_alt, edge2);

      TRACE("Direct_weight " << direct_weight << ", reverse_weight " << reverse_weight);
      return direct_weight > reverse_weight;
  }

public:
  AmbiguousPairInfoChecker(const Graph& graph, const PairedInfoIndexT<Graph>& index,
    AbstractPairInfoChecker<Graph> &standard_filter, double haplom_threshold,
    double relative_length_threshold, double relative_seq_threshold) :
        AbstractPairInfoChecker<Graph>(graph),
        standard_filter_(standard_filter),
        index_(index),
        haplom_threshold_(haplom_threshold),
        relative_length_threshold_(relative_length_threshold),
        relative_seq_threshold_(relative_seq_threshold) { }

  bool Check(const PairInfoT& info) {
      TRACE(this->graph_.int_id(info.first) << " " << this->graph_.int_id(info.second));
      if(EdgesAreFromSimpleBulgeWithAmbPI(info)){
        TRACE("Forward directed edges form a simple bulge");
        return InnerCheck(info);
      }

      if(EdgesAreFromSimpleBulgeWithAmbPI(BackwardInfo(info))){
          TRACE("Backward directed edges form a simple bulge");
          return InnerCheck(BackwardInfo(info));
      }

      TRACE("Edges do not form a bulge. Applying default checker");
      return standard_filter_.Check(info);
  }

private:
  OptEdgeId GetOtherSideOfSimpleBulge(EdgeId edge){
      auto edges = this->graph_.GetEdgesBetween(this->graph_.EdgeStart(edge),
              this->graph_.EdgeEnd(edge));
      TRACE("Number alternative edges -  " << edges.size());
      if(edges.size() == 1)
          return OptEdgeId();

      size_t edge_length = this->graph_.length(edge);
      Sequence edge_seq = this->graph_.EdgeNucls(edge);
      for(auto it_edge = edges.begin(); it_edge != edges.end(); it_edge++)
          if(*it_edge != edge){
              size_t it_edge_length = this->graph_.length(*it_edge);
              Sequence it_edge_seq = this->graph_.EdgeNucls(*it_edge);
              double length_ratio = double(min<size_t>(edge_length, it_edge_length)) /
                      double(max<size_t>(edge_length, it_edge_length));
              if(length_ratio >= relative_length_threshold_){
    //              size_t edit_dist = EditDistance(edge_seq, it_edge_seq);
    //              double seq_ratio = edit_dist / min<size_t> (edge_seq.size(), it_edge_seq.size());
                  return *it_edge;
              }
          }
      return OptEdgeId();
  }
};

template<class Graph>
class PairInfoFilter{
private:
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    typedef PairInfo<EdgeId> PairInfoT;

protected:
    AbstractPairInfoChecker<Graph> &pair_info_checker_;

public:
    PairInfoFilter(AbstractPairInfoChecker<Graph> &pair_info_checker) :
            pair_info_checker_(pair_info_checker)
    {}

    void Filter(PairedInfoIndexT<Graph>& index) {
        INFO("Start filtering; index size: " << index.size());
        //We can't filter while traversing, because Remove may invalidate iterators
        //So let's save edge pairs first
        using EdgePair = std::pair<EdgeId, EdgeId>;
        std::vector<EdgePair> pairs;
        for (auto i = pair_begin(index); i != pair_end(index); ++i)
            if (pair_info_checker_.Check(i.first(), i.second()))
                pairs.push_back({i.first(), i.second()});

        //TODO: implement fast removing of the whole set of points
        for (const auto& pair : pairs) {
            //Same thing with invalidation
            HistogramWithWeight hist;
            for (auto point : index[pair])
                if (!pair_info_checker_.Check(PairInfoT(pair.first, pair.second, point)))
                    hist.insert(point);
            //index.RemoveMany(pair_hist.first.first, pair_hist.first.second, pair_hist.second);
            for (const auto& point : hist)
                index.Remove(pair.first, pair.second, point);
        }

        INFO("Done filtering");
    }
};

}

}

#endif /* PAIR_INFO_FILTERS_HPP_ */

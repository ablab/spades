//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef PAIR_INFO_FILTERS_HPP_
#define PAIR_INFO_FILTERS_HPP_

namespace omnigraph {

template<class Graph>
class AbstractPairInfoFilter {

 private:
  typedef typename Graph::VertexId VertexId;
  typedef typename Graph::EdgeId EdgeId;
  typedef PairInfo<EdgeId> PairInfoT;
  typedef set<Point> Histogram;

 protected:
  virtual bool Check(const PairInfoT&) const {
    return true;
  }

  virtual bool Check(EdgeId, EdgeId) const {
    return true;
  }

  const Graph& graph_;

 public:
  AbstractPairInfoFilter(const Graph& graph) :
      graph_(graph)
  {
  }

  void Filter(const PairedInfoIndexT<Graph>& index,
                    PairedInfoIndexT<Graph>& new_index) const
  {
    TRACE("index size: " << index.size());
    for (auto it = index.begin(); it != index.end(); ++it) {
      const Histogram& infos = *it;
      const EdgeId& e1 = it.first();
      const EdgeId& e2 = it.second();
      if (Check(e1, e2)) {
        for (auto p_iter = infos.begin(); p_iter != infos.end(); ++p_iter) {
          const Point& point = *p_iter;
          if (Check(PairInfoT(e1, e2, point)))
            new_index.AddPairInfo(e1, e2, point, false);
        }
      }
    }
  }

  virtual ~AbstractPairInfoFilter()
  {
  }
};

template<class Graph>
class JumpingPairInfoChecker: public AbstractPairInfoFilter<Graph> {

 private:
  typedef typename Graph::VertexId VertexId;
  typedef typename Graph::EdgeId EdgeId;
  typedef PairInfo<EdgeId> PairInfoT;
  const GraphDistanceFinder<Graph> finder_;
  mutable size_t filtered_;
  mutable size_t passed_;

 public:
  JumpingPairInfoChecker(const Graph& graph, size_t short_is,
      size_t read_length, size_t delta) :
      AbstractPairInfoFilter<Graph>(graph), finder_(graph,
          short_is, read_length, delta),
          filtered_(0), passed_(0) {
  }

  virtual ~JumpingPairInfoChecker() {
    TRACE("In destructor of JumpingPairInfoChecker");
    TRACE("Filtered edge pairs " << filtered_);
    TRACE("Passed edge pairs " << passed_);
  }

 protected:
  virtual bool Check(EdgeId e1, EdgeId e2) const {
    vector<size_t> result1 = finder_.GetGraphDistancesLengths(e1, e2);
    vector<size_t> result2 = finder_.GetGraphDistancesLengths(e2, e1);
    bool result = result1.empty() && result2.empty();
    if (result)
      ++passed_;
    else
      ++filtered_;
    return result;
  }

 private:
  DECL_LOGGER("JumpingPairInfoChecker");
};

template<class Graph>
class PairInfoWeightFilter: public AbstractPairInfoFilter<Graph> {

 private:
  typedef typename Graph::EdgeId EdgeId;
  typedef PairInfo<EdgeId> PairInfoT;
  double weight_threshold_;

 public:
  PairInfoWeightFilter(const Graph& graph, double weight_threshold) :
    AbstractPairInfoFilter<Graph>(graph), weight_threshold_(weight_threshold) {
  }

 protected:
  virtual bool Check(const PairInfoT& info) const {
    return math::ge(info.weight(), weight_threshold_);
  }
};


template<class Graph>
class PairInfoWeightFilterWithCoverage: public AbstractPairInfoFilter<Graph> {

 private:
  typedef typename Graph::EdgeId EdgeId;
  typedef PairInfo<EdgeId> PairInfoT;
  double weight_threshold_;

 public:
  PairInfoWeightFilterWithCoverage(const Graph& graph, double weight_threshold) :
    AbstractPairInfoFilter<Graph>(graph), weight_threshold_(weight_threshold)
  {
  }

 protected:
  virtual bool Check(const PairInfoT& info) const {
    double info_weight = info.weight();   
    return math::ge(info_weight, weight_threshold_) 
       || (math::ge(info_weight, 0.1 * this->graph_.coverage(info.first))) 
       || (math::ge(info_weight, 0.1 * this->graph_.coverage(info.second)));
  }
};

}

#endif /* PAIR_INFO_FILTERS_HPP_ */

//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef PAIR_INFO_FILTERS_HPP_
#define PAIR_INFO_FILTERS_HPP_

namespace omnigraph {

namespace de {

template<class Graph>
class AbstractPairInfoFilter {

 private:
  typedef typename Graph::VertexId VertexId;
  typedef typename Graph::EdgeId EdgeId;
  typedef PairInfo<EdgeId> PairInfoT;

 protected:
  virtual bool Check(const PairInfoT&) const {
    return true;
  }

  virtual bool Check(EdgeId, EdgeId, const Point&) const {
    return true;
  }

  const Graph& graph_;

 public:
  AbstractPairInfoFilter(const Graph& graph) :
          graph_(graph) {}

  void Filter(PairedInfoIndexT<Graph>& index) const {
    TRACE("index size: " << index.size());
    for (auto it = index.begin(); it != index.end(); ++it) {
      const Histogram& infos = *it;
      const EdgeId& e1 = it.first();
      const EdgeId& e2 = it.second();

      for (auto p_iter = infos.begin(); p_iter != infos.end(); ++p_iter) {
        const Point& point = *p_iter;
        if (!Check(e1, e2, point))
          index.DeletePairInfo(e1, e2, point);
      }
    }
  }

  virtual ~AbstractPairInfoFilter()
  {
  }
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
  virtual bool Check(EdgeId, EdgeId, const Point& p) const {
    return math::ge(p.weight, weight_threshold_);
  }

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
    {}

 protected:
  virtual bool Check(EdgeId e1, EdgeId e2, const Point& p) const {
    double info_weight = p.weight;
    return math::ge(info_weight, weight_threshold_) ||
           math::ge(info_weight, 0.1 * this->graph_.coverage(e1)) ||
           math::ge(info_weight, 0.1 * this->graph_.coverage(e2));
  }

  virtual bool Check(const PairInfoT& info) const {
    return Check(info.first, info.second, info.point);
  }
};

}

}

#endif /* PAIR_INFO_FILTERS_HPP_ */

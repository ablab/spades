//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef SMOOTHING_DISTANCE_ESTIMATION_HPP_
#define SMOOTHING_DISTANCE_ESTIMATION_HPP_

#include "paired_info.hpp"
#include "omni/omni_utils.hpp"
#include "data_divider.hpp"
#include "peak_finder.hpp"
#include "extensive_distance_estimation.hpp"

namespace omnigraph {

namespace de {

template<class Graph>
class SmoothingDistanceEstimator: public ExtensiveDistanceEstimator<Graph> {
  typedef ExtensiveDistanceEstimator<Graph> base;

 public:
  SmoothingDistanceEstimator(const Graph& graph,
                             const PairedInfoIndexT<Graph>& histogram,
                             const GraphDistanceFinder<Graph>& dist_finder,
                             boost::function<double(int)> weight_f,
                             size_t linkage_distance,
                             size_t max_distance,
                             size_t threshold,
                             double range_coeff,
                             double delta_coeff,
                             size_t cutoff,
                             size_t min_peak_points,
                             double inv_density,
                             double percentage,
                             double derivative_threshold,
                             bool only_scaffolding = false) :
        base( graph, histogram, dist_finder, weight_f, linkage_distance, max_distance),
        threshold_(threshold),
        range_coeff_(range_coeff),
        delta_coeff_(delta_coeff),
        cutoff_((int) cutoff),
        min_peak_points_(min_peak_points),
        inv_density_(inv_density),
        percentage_(percentage),
        deriv_thr(derivative_threshold),
        only_scaffolding_(only_scaffolding),
        gap_distances(0)
    {
    }

    virtual ~SmoothingDistanceEstimator()
    {
    }

protected:
  typedef typename Graph::EdgeId EdgeId;
  typedef pair<EdgeId, EdgeId> EdgePair;
  typedef vector<pair<int, double> > EstimHist;
  typedef vector<PairInfo<EdgeId> > PairInfos;
  typedef vector<size_t> GraphLengths;

  virtual EstimHist EstimateEdgePairDistances(EdgePair /*ep*/,
                                const Histogram& /*raw_data*/,
                                const vector<size_t>& /*forward*/) const {
    VERIFY_MSG(false, "Sorry, the SMOOOOTHING estimator is not available anymore." <<
               "SPAdes is going to terminate");

    return EstimHist();
  }

private:
  typedef pair<size_t, size_t> Interval;

  size_t threshold_;
  double range_coeff_;
  double delta_coeff_;
  int    cutoff_;
  size_t min_peak_points_;
  double inv_density_;
  double percentage_;
  double deriv_thr;
  bool   only_scaffolding_;
  mutable size_t gap_distances;

  EstimHist FindEdgePairDistances(EdgePair ep,
                                  const Histogram& raw_hist) const
  {
    size_t first_len = this->graph().length(ep.first);
    size_t second_len = this->graph().length(ep.second);
    TRACE("Lengths are " << first_len << " " << second_len);
    Histogram data;
    for (auto I = raw_hist.begin(), E = raw_hist.end(); I != E; ++I)
    {
      const Point& p = *I;
      if (math::ge(2 * (long) rounded_d(p) + (long) second_len, (long) first_len))
        if ((long) rounded_d(p) >= (long) first_len)
          data.insert(p);
    }
    EstimHist result;
    double picture_weight = 0.;
    for (auto I = data.begin(), E = data.end(); I != E; ++I)
      picture_weight += I->weight;
    if (math::ls(picture_weight, 3.))
      return result;

    DataDivider<EdgeId> data_divider(threshold_,
                                     vector<Point>(data.begin(), data.end()));

    PairInfos infos;
    infos.reserve(data.size());
    const vector<Interval>& clusters =
      data_divider.DivideAndSmoothData(ep, infos, this->weight_f_);
    DEBUG("Seeking for distances");
    TRACE("size " << infos.size());

    for (size_t i = 0; i < clusters.size(); ++i) {
      size_t begin = clusters[i].first;
      size_t end = clusters[i].second;
      TRACE("begin " << begin << " at " << rounded_d(infos[begin])
        <<  ", " << " end " << end << " at " << rounded_d(infos[end - 1]));
      size_t data_length = rounded_d(infos[end - 1]) - rounded_d(infos[begin]) + 1;
      TRACE("data length " << data_length);
      if (end - begin > min_peak_points_) {
        size_t range = (size_t) math::round((double) data_length * range_coeff_);
        size_t delta = (size_t) math::round((double) data_length * delta_coeff_);
        PeakFinder<EdgeId> peakfinder(infos, begin, end, range, delta, percentage_, deriv_thr);
        DEBUG("Processing window : " << rounded_d(infos[begin])
                              << " " << rounded_d(infos[end - 1]));
        peakfinder.FFTSmoothing(cutoff_);
        TRACE("Listing peaks");
        const EstimHist& peaks = peakfinder.ListPeaks();
        //for (auto iter = peaks.begin(); iter != peaks.end(); ++iter) {
        //TRACE("PEAKS " << iter->first << " " << iter->second);
        //}
        if (peaks.size() == 0)
          continue;
        size_t index_of_max_weight = 0;
        for (size_t i = 0; i < peaks.size(); ++i)
          if (math::ls(peaks[index_of_max_weight].second, peaks[i].second))
            index_of_max_weight = i;
        result.push_back(peaks[index_of_max_weight]);
      }
    }

    if (result.size() == 0)
      return result;
    size_t index_of_max_weight = 0;
    for (size_t i = 0; i < result.size(); ++i)
      if (math::ls(result[index_of_max_weight].second, result[i].second))
        index_of_max_weight = i;

    EstimHist new_result;
    for (size_t i = 0; i < result.size(); ++i)
      if (result[i].second > .5 * result[index_of_max_weight].second)
        new_result.push_back(result[i]);
    return new_result;
    }

  virtual void ProcessEdge(EdgeId e1,
                           const typename PairedInfoIndexT<Graph>::InnerMap& inner_map,
                           PairedInfoIndexT<Graph>& result,
                           perf_counter& pc) const
  {
    pc.reset();
    set<EdgeId> second_edges;
    for (auto I = inner_map.begin(), E = inner_map.end(); I != E; ++I)
      second_edges.insert(I->first);

    vector<GraphLengths> lens_array = this->GetGraphDistancesLengths(e1, second_edges);

    size_t i = 0;
    for (auto I = inner_map.begin(), E = inner_map.end(); I != E; ++I) {
      EdgeId e2 = I->first;
      EdgePair ep(e1, e2);
      if (ep <= this->ConjugatePair(ep)) {
        const GraphLengths& forward = lens_array[i++];

        TRACE("Processing edge pair " << this->graph().int_id(e1)
                               << " " << this->graph().int_id(e2));
        Histogram hist = I->second;
        EstimHist estimated;
        //DEBUG("Extending paired information");
        //DEBUG("Extend left");
        //this->base::ExtendInfoLeft(e1, e2, hist, 1000);
        DEBUG("Extend right");
        this->ExtendInfoRight(e1, e2, hist, 1000);
        if (forward.size() == 0) {
          estimated = FindEdgePairDistances(ep, hist);
          ++gap_distances;
        }
        else if (forward.size() > 0 && (!only_scaffolding_)) {
          estimated = this->base::EstimateEdgePairDistances(ep, hist, forward);
        }
        DEBUG(gap_distances << " distances between gap edge pairs have been found");
        Histogram res = this->ClusterResult(ep, estimated);
        this->AddToResult(res, ep, result);
        this->AddToResult(this->ConjugateInfos(ep, res), this->ConjugatePair(ep), result);
      }
    }
  }

  bool IsTipTip(EdgeId e1, EdgeId e2) const {
    return (this->graph().OutgoingEdgeCount(this->graph().EdgeEnd(e1)) == 0 &&
            this->graph().IncomingEdgeCount(this->graph().EdgeEnd(e1)) == 1 &&
            this->graph().IncomingEdgeCount(this->graph().EdgeStart(e2)) == 0 &&
            this->graph().OutgoingEdgeCount(this->graph().EdgeStart(e2)) == 1);
  }

  virtual const string Name() const {
    static const string my_name = "SMOOTHING";
    return my_name;
  }

  DECL_LOGGER("SmoothingDistanceEstimator")
};

}

}

#endif /* SMOOTHING_DISTANCE_ESTIMATION_HPP_ */

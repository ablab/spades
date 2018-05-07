//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef SMOOTHING_DISTANCE_ESTIMATION_HPP_
#define SMOOTHING_DISTANCE_ESTIMATION_HPP_

#include "weighted_distance_estimation.hpp"
#include "data_divider.hpp"
#include "peak_finder.hpp"

namespace omnigraph {

namespace de {

class SmoothingDistanceEstimator : public WeightedDistanceEstimator {
    //FIXME configure
    static const size_t OVERLAP_TOLERANCE = 1000;
protected:
    typedef WeightedDistanceEstimator base;
    typedef typename base::InPairedIndex InPairedIndex;
    typedef typename base::OutPairedIndex OutPairedIndex;
    typedef typename base::InHistogram InHistogram;
    typedef typename base::OutHistogram OutHistogram;
    typedef typename InPairedIndex::Histogram TempHistogram;

public:
    SmoothingDistanceEstimator(const debruijn_graph::Graph &graph,
                               const InPairedIndex &histogram,
                               const GraphDistanceFinder &dist_finder,
                               std::function<double(int)> weight_f,
                               size_t linkage_distance, size_t max_distance, size_t threshold,
                               double range_coeff, double delta_coeff,
                               size_t cutoff,
                               size_t min_peak_points,
                               double inv_density,
                               double percentage,
                               double derivative_threshold) :
            base(graph, histogram, dist_finder, weight_f, linkage_distance, max_distance),
            threshold_(threshold),
            range_coeff_(range_coeff),
            delta_coeff_(delta_coeff),
            cutoff_((int) cutoff),
            min_peak_points_(min_peak_points),
            inv_density_(inv_density),
            percentage_(percentage),
            deriv_thr(derivative_threshold),
            gap_distances(0) { }

    virtual ~SmoothingDistanceEstimator() { }

protected:
    typedef pair<debruijn_graph::EdgeId, debruijn_graph::EdgeId> EdgePair;
    typedef vector<pair<int, double> > EstimHist;
    typedef vector<PairInfo<debruijn_graph::EdgeId> > PairInfos;
    typedef vector<size_t> GraphLengths;

    EstimHist EstimateEdgePairDistances(EdgePair /*ep*/,
                                        const InHistogram & /*raw_data*/,
                                        const vector<size_t> & /*forward*/) const override {
        VERIFY_MSG(false, "Sorry, the SMOOOOTHING estimator is not available anymore." <<
                          "SPAdes is going to terminate");

        return EstimHist();
    }

private:
    typedef pair<size_t, size_t> Interval;

    size_t threshold_;
    double range_coeff_;
    double delta_coeff_;
    int cutoff_;
    size_t min_peak_points_;
    double inv_density_;
    double percentage_;
    double deriv_thr;
    mutable size_t gap_distances;

    EstimHist FindEdgePairDistances(EdgePair ep,
                                    const TempHistogram &raw_hist) const;

    void ProcessEdge(debruijn_graph::EdgeId e1,
                     const InPairedIndex &pi,
                     PairedInfoBuffer<debruijn_graph::Graph> &result) const override;

    bool IsTipTip(debruijn_graph::EdgeId e1, debruijn_graph::EdgeId e2) const;

    void ExtendInfoRight(debruijn_graph::EdgeId e1, debruijn_graph::EdgeId e2, TempHistogram &data,
                         size_t max_shift) const {
        ExtendRightDFS(e1, e2, data, 0, max_shift);
    }

    void MergeInto(const InHistogram &what, TempHistogram &where, int shift) const;

    void ExtendRightDFS(const debruijn_graph::EdgeId &first, debruijn_graph::EdgeId current, TempHistogram &data,
                        int shift, size_t max_shift) const;

    const string Name() const override {
        static const string my_name = "SMOOTHING";
        return my_name;
    }

    DECL_LOGGER("SmoothingDistanceEstimator")
};

}

}

#endif /* SMOOTHING_DISTANCE_ESTIMATION_HPP_ */

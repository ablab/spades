//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef SMOOTHING_DISTANCE_ESTIMATION_HPP_
#define SMOOTHING_DISTANCE_ESTIMATION_HPP_

#include "paired_info.hpp"
#include "omni_utils.hpp"
#include "data_divider.hpp"
#include "peak_finder.hpp"
#include "distance_estimation.hpp"

namespace omnigraph {

template<class Graph>
class SmoothingDistanceEstimator: public DistanceEstimator<Graph> {
	typedef DistanceEstimator<Graph> base;
	typedef typename Graph::EdgeId EdgeId;
	typedef pair<int, int> Interval;


    size_t threshold_;
    double range_coeff_;
    double delta_coeff_;
	int    cutoff_;
    size_t min_peak_points_;    
    double inv_density_;
    double percentage_;
    double derivative_threshold_;

    int round(double x) const { 
        int res = (int) (abs(x) + 0.5 + 1e-9);
        if (x < 0)
            res = -res;
        return res;
    }

protected:
    vector<pair<size_t, double> > EstimateEdgePairDistances(EdgeId first, EdgeId second, const vector<PairInfo<EdgeId> >& data, const vector<size_t>& forward) const {
        vector<pair<size_t, double> > result;
        if (data.size() <= 1) 
            return result;

        DataDivider data_divider(threshold_);
		const vector<Interval>& clusters = data_divider.DivideData<EdgeId>(data);
		size_t cur = 0;
        stringstream ss;
        for (size_t i = 0; i < forward.size(); i++){
            ss << forward[i] << " ";
        }
        DEBUG("Possible distances : " << ss.str());

		for (size_t i = 0; i < clusters.size(); i++) {
            int begin = clusters[i].first;
            int end = clusters[i].second;
            if (end - begin >= (int) min_peak_points_) {
                size_t data_length = rounded_d(data[end - 1]) - rounded_d(data[begin]) + 1;
                while (cur < forward.size() && (int) forward[cur] < rounded_d(data[begin]))
                    cur++;
                if (cur == forward.size()) {
                    DEBUG("BREAKING " << rounded_d(data[begin]));
                    break;
                }
                if ((int) forward[cur] > rounded_d(data[end - 1])) 
                    continue;
                else {
                    PeakFinder<EdgeId> peakfinder(data, begin, end, round(data_length * range_coeff_), round(data_length * delta_coeff_), percentage_, derivative_threshold_);
                    DEBUG("Processing window : " << rounded_d(data[begin]) << " " << rounded_d(data[end - 1]));
                    peakfinder.FFTSmoothing(cutoff_);
                    if ( (cur + 1) == forward.size() || (int) forward[cur + 1] > rounded_d(data[end - 1] )) {
                        if (round(inv_density_ * (end - begin)) > (int) data_length) {
                            result.push_back(make_pair(forward[cur], peakfinder.GetNormalizedWeight()));       // default weight is one
                            DEBUG("Pair made " << forward[cur] << " " << peakfinder.GetNormalizedWeight());
                        }
                        cur++;
                    } 
                    else {
                        while (cur < forward.size() && (int) forward[cur] <= rounded_d(data[end - 1] )) {
                            if (peakfinder.IsPeak(forward[cur])) { 
                                result.push_back(make_pair(forward[cur], peakfinder.GetNormalizedWeight()));
                                DEBUG("Pair made " << forward[cur] << " " << peakfinder.GetNormalizedWeight());
                            }   
                            cur++;
                        }
                    }
                }
			}
		}
		return result;
	}

	virtual void ProcessEdgePair(EdgeId first, EdgeId second, const vector<PairInfo<EdgeId>>& data, PairedInfoIndex<Graph> &result) const {
		if (make_pair(first, second) <= ConjugatePair(first, second)) {
			vector<size_t> forward = this->GetGraphDistances(first, second);
		    vector<pair<size_t, double> > estimated;
            TRACE("Processing edge pair " << first << " " << second);
            if (forward.size() > 0) 
                estimated = EstimateEdgePairDistances(first, second, data, forward);
            else 
                estimated = FindEdgePairDistances(data);

			vector<PairInfo<EdgeId>> res = this->ClusterResult(first, second, estimated);
			this->AddToResult(result, res);
			this->AddToResult(result, ConjugateInfos(res));
		}
	}

public:
	SmoothingDistanceEstimator(const Graph& graph, const PairedInfoIndex<Graph>& histogram, const GraphDistanceFinder<Graph>& dist_finder, size_t linkage_distance, size_t threshold, double range_coeff, double delta_coeff, size_t cutoff, size_t min_peak_points, double inv_density, double percentage, double derivative_threshold) : 
        base(graph, histogram, dist_finder, linkage_distance, 0),
        threshold_(threshold),
        range_coeff_(range_coeff), 
        delta_coeff_(delta_coeff), 
        cutoff_(cutoff), 
        min_peak_points_(min_peak_points), 
        inv_density_(inv_density), 
        percentage_(percentage), 
        derivative_threshold_(derivative_threshold) 
    {
    }

    vector<pair<size_t, double> > FindEdgePairDistances(const vector<PairInfo<EdgeId> >& data) const {
        vector<pair<size_t, double> > result;
        if (data.size() <= 1)
            return result;
		
        DataDivider data_divider(threshold_);
		const vector<Interval>& clusters = data_divider.DivideData<EdgeId>(data);
        DEBUG("Seeking for distances");
		
        for (size_t i = 0; i < clusters.size(); i++) {
            int begin = clusters[i].first;
            int end = clusters[i].second;
            size_t data_length = rounded_d(data[end - 1]) - rounded_d(data[begin]) + 1;
            if (end - (int) begin > (int) min_peak_points_) {
                PeakFinder<EdgeId> peakfinder(data, begin, end, round(data_length * range_coeff_), round(data_length * delta_coeff_), percentage_, derivative_threshold_);
                DEBUG("Processing window : " << rounded_d(data[begin]) << " " << rounded_d(data[end - 1]));

                TRACE("Smoothing via FFT");
                
                peakfinder.FFTSmoothing(cutoff_);
                //peakfinder.PrintStats();

                TRACE("Listing peaks");

                vector<pair<size_t, double> > peaks = peakfinder.ListPeaks();

                for (auto iter = peaks.begin(); iter != peaks.end(); iter++) {
                    result.push_back(*iter);
                    TRACE("PEAKS " << iter->first << " " << iter->second);
                }
            }
        }
        return result;
    }

    private:
	    DECL_LOGGER("SmoothingDistanceEstimator")
};

}

#endif /* SMOOTHING_DISTANCE_ESTIMATION_HPP_ */

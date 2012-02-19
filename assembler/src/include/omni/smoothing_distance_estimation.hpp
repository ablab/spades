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
#include "extensive_distance_estimation.hpp"

namespace omnigraph {

template<class Graph>
class SmoothingDistanceEstimator: public ExtensiveDistanceEstimator<Graph> {
	typedef ExtensiveDistanceEstimator<Graph> base;
	typedef typename Graph::EdgeId EdgeId;
	typedef pair<size_t, size_t> Interval;

    size_t threshold_;
    double range_coeff_;
    double delta_coeff_;
	int    cutoff_;
    size_t min_peak_points_;    
    double inv_density_;
    double percentage_;
    double derivative_threshold_;
    mutable size_t gap_distances;

    int round(double x) const { 
        int res = (int) (abs(x) + 0.5 + 1e-9);
        if (x < 0)
            res = -res;
        return res;
    }

protected: 
vector<pair<int, double> > EstimateEdgePairDistances(EdgeId first, EdgeId second, const vector<PairInfo<EdgeId> >& raw_data, const vector<size_t>& forward) const {
        VERIFY(false);
        vector<PairInfo<EdgeId> > data;
        size_t first_len = this->graph().length(raw_data[0].first);
        size_t second_len = this->graph().length(raw_data[0].second);
        for (size_t i = 0 ; i < raw_data.size(); ++i) {
            if (math::ge(2. * rounded_d(raw_data[i]) + second_len, (double) first_len))
                data.push_back(raw_data[i]);
        }
        
        vector<pair<int, double> > result;
        double picture_weight = 0.;
        for (size_t i = 0; i < data.size(); ++i) 
            picture_weight += data[i].weight;
        if (math::ls(picture_weight, 40000000.))
            return result;
        // dirty hack against false negative
        //if (data.size() > 1 && forward.size() == 1) {
            //double total_weight = 0.;
            //for (size_t i = 0; i < data.size(); ++i) 
                //total_weight += data[i].weight;
            //result.push_back(make_pair(forward[0], total_weight));
            //return result;
        //}

        vector<PairInfo<EdgeId> > new_data;

        DataDivider data_divider(threshold_);
		const vector<Interval>& clusters = data_divider.DivideAndSmoothData<EdgeId>(data, new_data, this->weight_f_);
		size_t cur = 0;
        stringstream ss;
        for (size_t i = 0; i < forward.size(); i++){
            ss << forward[i] << " ";
        }
        DEBUG("Possible distances : " << ss.str());

		for (size_t i = 0; i < clusters.size(); i++) {
            size_t begin = clusters[i].first;
            size_t end = clusters[i].second;
            TRACE("begin " << begin << " at " << rounded_d(new_data[begin]) <<  ", " << "end " << end << " at " << rounded_d(new_data[end - 1]));
            size_t data_length = rounded_d(new_data[end - 1]) - rounded_d(new_data[begin]) + 1;
            TRACE("data length " << data_length);
            if ((end - begin) >= min_peak_points_) {
                while (cur < forward.size() && (int) forward[cur] < rounded_d(new_data[begin]))
                    cur++;
                if (cur == forward.size()) {
                    DEBUG("Breaking " << rounded_d(new_data[begin]));
                    break;
                }
                if ((int) forward[cur] > rounded_d(new_data[end - 1])) 
                    continue;
                else {
                    PeakFinder<EdgeId> peakfinder(new_data, begin, end, round(data_length * range_coeff_), round(data_length * delta_coeff_), percentage_, derivative_threshold_);
                    DEBUG("Processing window : " << rounded_d(new_data[begin]) << " " << rounded_d(new_data[end - 1]));
                    peakfinder.FFTSmoothing(cutoff_);
                    if ( (cur + 1) == forward.size() || (int) forward[cur + 1] > rounded_d(new_data[end - 1] )) {
                        if (round(inv_density_ * (end - begin)) > (int) data_length) {
                            result.push_back(make_pair(forward[cur], peakfinder.GetNormalizedWeight()));       // default weight is one
                            DEBUG("Pair made " << forward[cur] << " " << peakfinder.GetNormalizedWeight());
                        }
                        cur++;
                    } 
                    else {
                        while (cur < forward.size() && (int) forward[cur] <= rounded_d(new_data[end - 1] )) {
                            if (peakfinder.IsPeak(forward[cur])) { 
                                result.push_back(make_pair(forward[cur], peakfinder.weight()));
                                DEBUG("Pair made " << forward[cur] << " " << peakfinder.weight());
                            }   
                            cur++;
                        }
                    }
                }
			}
		}
		return result;
	}

    bool IsTipTip(EdgeId first, EdgeId second) const {
            return (this->graph().OutgoingEdgeCount(this->graph().EdgeEnd(first)) == 0 && this->graph().IncomingEdgeCount(this->graph().EdgeEnd(first)) == 1 && 
                this->graph().IncomingEdgeCount(this->graph().EdgeStart(second)) == 0 && this->graph().OutgoingEdgeCount(this->graph().EdgeStart(second)) == 1); 
    }

	virtual void ProcessEdgePair(EdgeId first, EdgeId second, const vector<PairInfo<EdgeId>>& raw_data, PairedInfoIndex<Graph> &result) const {
		if (make_pair(first, second) <= this->ConjugatePair(first, second)) {
			vector<size_t> forward = this->GetGraphDistancesLengths(first, second);
            vector<size_t> backward = this->GetGraphDistancesLengths(second, first);
            TRACE("Processing edge pair " << this->graph().int_id(first) << " " << this->graph().int_id(second));
            vector<PairInfo<EdgeId>> data = raw_data;
            //DEBUG("Extending paired information");
            //double weight_0 = this->WeightSum(data);
            //DEBUG("Extend left");
            //this->ExtendInfoLeft(first, second, data);
            //DEBUG("Extend right");
            //this->ExtendInfoRight(first, second, data);
            
            //DEBUG("Weight increased " << (WeightSum(data) - weight_0));
            

            //TODO: clean up
		    vector<pair<int, double> > estimated;
            if (forward.size() == 0 && backward.size() == 0 && IsTipTip(first, second)) 
            {
                estimated = FindEdgePairDistances(data);
                ++gap_distances;
            }
            else if (forward.size() + backward.size() > 0) { 
                DEBUG("Extending paired information");
                DEBUG("Extend left");
                this->ExtendInfoLeft(first, second, data);
                DEBUG("Extend right");
                this->ExtendInfoRight(first, second, data);
                
                estimated = this->base::EstimateEdgePairDistances(first, second, data, forward);
            }

            DEBUG(gap_distances << " distances between gap edge pairs have been found");

			vector<PairInfo<EdgeId>> res = this->ClusterResult(first, second, estimated);
			this->AddToResult(result, res);
			this->AddToResult(result, this->ConjugateInfos(res));

		}
	}

public:
	SmoothingDistanceEstimator(const Graph& graph, const PairedInfoIndex<Graph>& histogram, const GraphDistanceFinder<Graph>& dist_finder, boost::function<double(int)> weight_f, size_t linkage_distance, size_t max_distance, size_t threshold, double range_coeff, double delta_coeff, size_t cutoff, size_t min_peak_points, double inv_density, double percentage, double derivative_threshold) : 
        base(graph, histogram, dist_finder, weight_f, linkage_distance, max_distance),
        threshold_(threshold),
        range_coeff_(range_coeff), 
        delta_coeff_(delta_coeff), 
        cutoff_(cutoff), 
        min_peak_points_(min_peak_points), 
        inv_density_(inv_density), 
        percentage_(percentage), 
        derivative_threshold_(derivative_threshold),
        gap_distances(0)
    {
    }

    vector<pair<int, double> > FindEdgePairDistances(const vector<PairInfo<EdgeId> >& raw_data) const {
        size_t first_len = this->graph().length(raw_data[0].first);
        size_t second_len = this->graph().length(raw_data[0].second);
        TRACE("Lengths are " << first_len << " " << second_len);
        vector<PairInfo<EdgeId> > data;
        for (size_t i = 0 ; i < raw_data.size(); ++i) {
            if (math::ge(2. * rounded_d(raw_data[i]) + second_len, (double) first_len))
                if (rounded_d(raw_data[i]) >= (int) first_len)
                    data.push_back(raw_data[i]);
        }
        vector<pair<int, double> > result;
        double picture_weight = 0;
        for (size_t i = 0; i < data.size(); ++i) 
            picture_weight += data[i].weight;
        if (math::ls(picture_weight, 3.))
            return result;
		
        vector<PairInfo<EdgeId> > new_data;
        DataDivider data_divider(threshold_);
		const vector<Interval>& clusters = data_divider.DivideAndSmoothData<EdgeId>(data, new_data, this->weight_f_);
        DEBUG("Seeking for distances");
        TRACE("size " << new_data.size());
		
        for (size_t i = 0; i < clusters.size(); i++) {
            size_t begin = clusters[i].first;
            size_t end = clusters[i].second;
            TRACE("begin " << begin << " at " << rounded_d(new_data[begin]) <<  ", " << " end " << end << " at " << rounded_d(new_data[end - 1]));
            size_t data_length = rounded_d(new_data[end - 1]) - rounded_d(new_data[begin]) + 1;
            TRACE("data length " << data_length);
            if (end - begin > min_peak_points_) {
                PeakFinder<EdgeId> peakfinder(new_data, begin, end, round(data_length * range_coeff_), round(data_length * delta_coeff_), percentage_, derivative_threshold_);
                DEBUG("Processing window : " << rounded_d(new_data[begin]) << " " << rounded_d(new_data[end - 1]));

                TRACE("Smoothing via FFT");
                
                peakfinder.FFTSmoothing(cutoff_);

                TRACE("Listing peaks");

                const vector<pair<int, double> >& peaks = peakfinder.ListPeaks();

                TRACE("Listed");
                int index_of_max_weight = 0;
                for (auto iter = peaks.begin(); iter != peaks.end(); ++iter) {
                    //result.push_back(*iter);
                    TRACE("PEAKS " << iter->first << " " << iter->second);
                }
                if (peaks.size() == 0) 
                    continue;
                for (size_t i = 0; i < peaks.size(); ++i) 
                    if (math::ls(peaks[index_of_max_weight].second, peaks[i].second))
                        index_of_max_weight = i;
                if (index_of_max_weight >= 0) 
                    result.push_back(peaks[index_of_max_weight]);
            }
        }
        if (result.size() == 0)
            return result;
        int index_of_max_weight = 0;
        for (size_t i = 0; i < result.size(); ++i) 
            if (math::ls(result[index_of_max_weight].second, result[i].second))
                index_of_max_weight = i;
        vector<pair<int, double> > new_result;
        for (size_t i = 0; i < result.size(); ++i)
            if (result[i].second > .5 * result[index_of_max_weight].second)
                new_result.push_back(result[i]);
        return new_result;
    }

    private:
	    DECL_LOGGER("SmoothingDistanceEstimator")
};

}

#endif /* SMOOTHING_DISTANCE_ESTIMATION_HPP_ */

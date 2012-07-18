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
	int cutoff_;
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
    vector<pair<size_t, double> > EstimateEdgePairDistances(EdgeId first, EdgeId second, vector<PairInfo<EdgeId> > data, const vector<size_t> forward) const {
        vector<pair<size_t, double> > result;
        if (data.size() <= 1) 
            return result;

        DataDivider data_divider(threshold_);
		vector<Interval> clusters = data_divider.divideData<EdgeId>(data);
		size_t cur = 0;
        stringstream ss;
        for (size_t i = 0; i < forward.size(); i++){
            ss << forward[i] << " ";
        }
        DEBUG("Possible distances : " << ss.str());
        cout << " MINPEAKPOINTS " << min_peak_points_ << endl;
        cout << " cluster " << clusters.size() << endl;

		for (size_t i = 0; i < clusters.size(); i++) {
            int begin = clusters[i].first;
            int end = clusters[i].second;
            cout << end - begin << " " << min_peak_points_ << endl;
            if (end - begin >= (int) min_peak_points_) {
                cout << 1 << endl;
                size_t data_length = rounded_d(data[end - 1]) - rounded_d(data[begin]) + 1;
                while (cur < forward.size() && (int) forward[cur] < rounded_d(data[begin]))
                    cur++;
                if (cur == forward.size()) {
                    DEBUG("BREAKING " << rounded_d(data[begin]));
                    break;
                }
                cout << "Forward " << forward[cur] << " Dist " << rounded_d(data[end -1]) << endl;
                if ((int) forward[cur] > rounded_d(data[end - 1])) 
                    continue;
                else {
                    PeakFinder<EdgeId> peakfinder(data, begin, end, round(data_length * range_coeff_), round(data_length * delta_coeff_), percentage_, derivative_threshold_);
                    DEBUG("Processing window : " << rounded_d(data[begin]) << " " << rounded_d(data[end - 1]));
                    peakfinder.FFTSmoothing(cutoff_);
                    if ( (cur + 1) == forward.size() || (int) forward[cur + 1] > rounded_d(data[end - 1] )) {
                        if (round(inv_density_ * (end - begin)) > (int) data_length) {
                            result.push_back(make_pair(forward[cur], peakfinder.getNormalizedWeight()));       // default weight is one
                            DEBUG("Pair made " << forward[cur] << " " << peakfinder.getNormalizedWeight());
                        }
                        cur++;
                    } 
                    else {
                        while (cur < forward.size() && (int) forward[cur] <= rounded_d(data[end - 1] )) {
                            if (peakfinder.isPeak(forward[cur])) { 
                                result.push_back(make_pair(forward[cur], peakfinder.getNormalizedWeight()));
                                DEBUG("Pair made " << forward[cur] << " " << peakfinder.getNormalizedWeight());
                            }   
                            cur++;
                        }
                    }
                }
			}
		}
		return result;
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

	//vector<pair<int, double> > FindEdgePairDistances(vector<PairInfo<EdgeId> > data, vector<size_t> forward = NULL) {
        //vector<pair<size_t, double> > result;
        //if (data.size() <= 1)
            //return result;
		
        //vector<Interval> clusters = divideData(data);
		
        //size_t cur = 0;
        //std::stringstream ss;
        //for (size_t i = 0; i < forward.size(); i++){
            //ss << forward[i] << " ";
        //}
        //INFO("Possible distances : " << ss.str());

		//for (size_t i = 0; i < clusters.size(); i++) {
            //int begin = clusters[i].first;
            //int end = clusters[i].second;
            //size_t data_length = rounded_d(data[end - 1]) - rounded_d(data[begin]) + 1;
            //if (end - begin > min_peak_points_) {
                //while ((cur < forward.size()) && ( (int) forward[cur] < rounded_d(data[begin])))
					//cur++;
                //if (cur == forward.size()) 
                    //break;
                
                //if ((int) forward[cur] > rounded_d(data[end - 1])) 
                    //continue;
                //else {
                    //PeakFinder<EdgeId> peakfinder(data, begin, end);
                    //DEBUG("Processing window : " << rounded_d(data[begin]) << " " << rounded_d(data[end - 1]));
                    
                    //peakfinder.FFTSmoothing(cutoff_);

                    //vector<pair<int, double> > peaks = peakfinder.ListPeaks();

                    //for (auto iter = peaks.begin(); iter != peaks.end(); iter++) 
                        //result.push_back(*iter);
                //}
			//}
		//}
		//return result;
	//}

	//void Estimate(PairedInfoIndex<Graph> &result) {
		//for (auto iterator = this->histogram().begin(); iterator != this->histogram().end(); ++iterator) {
			//vector<PairInfo<EdgeId> > data = *iterator;
			//EdgeId first = data[0].first;
			//EdgeId second = data[0].second;
            //int firstNumber =  this->graph().int_ids().ReturnIntId(first);
            //int secondNumber =  this->graph().int_ids().ReturnIntId(second);

            //DEBUG("Estimating edges number : " << firstNumber << " " << secondNumber); 
            //vector<size_t> forward = this->GetGraphDistances(first, second);
			//vector<pair<size_t, double> > estimated = EstimateEdgePairDistances(data, forward);
			//vector<PairInfo<EdgeId> > clustered = this->ClusterResult(first, second, estimated);
			//this->AddToResult(result, clustered);
		//}
	//}
    private:
	    DECL_LOGGER("SmoothingDistanceEstimator")

};

}

#endif /* SMOOTHING_DISTANCE_ESTIMATION_HPP_ */

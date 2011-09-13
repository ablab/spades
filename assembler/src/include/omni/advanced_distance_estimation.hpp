#ifndef ADVANCED_DISTANCE_ESTIMATION_HPP_
#define ADVANCED_DISTANCE_ESTIMATION_HPP_

#include "paired_info.hpp"
#include "omni_utils.hpp"
#include "data_divider.hpp"
#include "peak_finder.hpp"
#include "distance_estimation.hpp"

namespace omnigraph {

template<class Graph>
class AdvancedDistanceEstimator: public DistanceEstimator<Graph> {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef std::pair<int, int> interval;
	IdTrackHandler<Graph> &int_ids_;


    double range_coeff_;
    double delta_coeff_;
	int cutoff_;
    size_t minpeakpoints_;    
    double inv_density_;
    double percentage_;
    double derivative_threshold_;


	vector<pair<size_t, double> > EstimateEdgePairDistances(vector<PairInfo<EdgeId> > data, vector<size_t> forward) {
        vector<pair<size_t, double> > result;
        if (data.size() <= 1) return result;
		std::vector<interval> clusters = divideData(data);
		size_t cur = 0;
        std::stringstream ss;
        for (size_t i = 0; i < forward.size(); i++){
            ss << forward[i] << " ";
        }
        DEBUG("Possible distances : " << ss.str());

		for (size_t i = 0; i < clusters.size(); i++) {
            size_t begin = clusters[i].first;
            size_t end = clusters[i].second;
            if (end - begin >= minpeakpoints_) {
                size_t data_length = rounded_d(data[end - 1]) - rounded_d(data[begin]) + 1;
                while ((cur<forward.size()) && (((int)forward[cur]) < rounded_d(data[begin])))
					cur++;
                if (cur == forward.size()) {
                    DEBUG("BREAKING " << rounded_d(data[begin]));
                    break;
                }
                if ((int) forward[cur] > rounded_d(data[end - 1])) continue;
                PeakFinder peakfinder(data, begin, end, data_length*range_coeff_, data_length*delta_coeff_, percentage_, derivative_threshold_);
				DEBUG("Processing window : " << rounded_d(data[begin]) << " " << rounded_d(data[end - 1]));
				peakfinder.FFTSmoothing(cutoff_);
                if ( ( (cur + 1) == forward.size()) || ( (int) forward[cur + 1] > rounded_d(data[end - 1]))){
                    if (inv_density_*(end - begin) > data_length){
                        result.push_back(make_pair(forward[cur], 1));
                        DEBUG("Pair made " << forward[cur]);
                    }
                    cur++;
                }else{
                
                    while (cur<forward.size() && ((int)forward[cur] <= rounded_d(data[end - 1]))) {
					    if (peakfinder.isPeak(forward[cur])){ 
                            result.push_back(make_pair(forward[cur], 1));
                            DEBUG("Pair made " << forward[cur]);
                        }   
					    cur++;
				    }
                }
			}
		}
		return result;
	}

public:
	AdvancedDistanceEstimator(Graph &graph, PairedInfoIndex<Graph> &histogram, IdTrackHandler<Graph> &int_ids, size_t insert_size, size_t read_length, size_t delta, size_t linkage_distance, size_t max_distance, size_t threshold, double range_coeff, double delta_coeff, 
    size_t cutoff, size_t minpeakpoints, double inv_density, double percentage, double derivative_threshold) : 
    DistanceEstimator<Graph>::DistanceEstimator(graph, histogram, insert_size, read_length, delta, linkage_distance, max_distance), 
    int_ids_(int_ids), range_coeff_(range_coeff), delta_coeff_(delta_coeff), cutoff_(cutoff), minpeakpoints_(minpeakpoints), inv_density_(inv_density), percentage_(percentage), derivative_threshold_(derivative_threshold){  
	        INFO("Advanced Estimator started");
            Threshold = threshold;
    }

	virtual ~AdvancedDistanceEstimator() {
	}

	vector<pair<int, double> > FindEdgePairDistances(vector<PairInfo<EdgeId> > data, vector<size_t> forward = NULL) {
        vector<pair<size_t, double> > result;
        if (data.size() <= 1) return result;
		std::vector<interval> clusters = divideData(data);
		size_t cur = 0;
        std::stringstream ss;
        for (size_t i = 0; i < forward.size(); i++){
            ss << forward[i] << " ";
        }
        INFO("Possible distances : " << ss.str());

		for (size_t i = 0; i < clusters.size(); i++) {
            size_t begin = clusters[i].first;
            size_t end = clusters[i].second;
            size_t data_length = rounded_d(data[end - 1]) - rounded_d(data[begin]) + 1;
            if (end - begin > minpeakpoints_) {
                while ((cur < forward.size()) && ( (int) forward[cur] < rounded_d(data[begin])))
					cur++;
                if (cur == forward.size()) {
                    break;
                }
                if ((int) forward[cur] > rounded_d(data[end - 1])) continue;
                PeakFinder peakfinder(data, begin, end);
				DEBUG("Processing window : " << rounded_d(data[begin]) << " " << rounded_d(data[end - 1]));
				
                peakfinder.FFTSmoothing(cutoff_);

                vector<pair<int, double> > peaks = peakfinder.ListPeaks();
                for (auto iter = peaks.begin(); iter != peaks.end(); iter++) result.push_back(*iter);
			}
		}
		return result;
	}

	virtual void Estimate(PairedInfoIndex<Graph> &result) {
		for (auto iterator = this->histogram_.begin(); iterator != this->histogram_.end(); ++iterator) {
			vector<PairInfo<EdgeId> > data = *iterator;
			EdgeId first = data[0].first;
			EdgeId second = data[0].second;
            int firstNumber =  int_ids_.ReturnIntId(first); 
            int secondNumber =  int_ids_.ReturnIntId(second); 

            DEBUG("Estimating edges number : " << firstNumber << " " << secondNumber); 
            vector<size_t> forward = this->GetGraphDistances(first, second);
			vector<pair<size_t, double> > estimated = EstimateEdgePairDistances(data, forward);
			vector<PairInfo<EdgeId> > clustered = this->ClusterResult(first, second, estimated);
			this->AddToResult(result, clustered);
		}
	}
};

}

#endif /* ADVANCED_DISTANCE_ESTIMATION_HPP_ */

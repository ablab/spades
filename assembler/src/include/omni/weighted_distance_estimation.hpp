//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef WEIGHTED_DISTANCE_ESTIMATION_HPP_
#define WEIGHTED_DISTANCE_ESTIMATION_HPP_

#include "xmath.h"
#include "paired_info.hpp"
#include "omni_utils.hpp"
#include "distance_estimation.hpp"

namespace omnigraph {

template<class Graph>
class WeightedDistanceEstimator: public DistanceEstimator<Graph> {
	typedef DistanceEstimator<Graph> base;
	typedef typename Graph::EdgeId EdgeId;

protected:
	const size_t max_distance_;

    boost::function<double(int)> weight_f_;

	virtual vector<pair<size_t, double>> EstimateEdgePairDistances(EdgeId first, EdgeId second,
			vector<PairInfo<EdgeId>> data,
			const vector<size_t> raw_forward) const {

        size_t first_len = this->graph().length(first);
        size_t second_len = this->graph().length(second);
        
		vector<pair<size_t, double>> result;
		int maxD = rounded_d(data.back());
		int minD = rounded_d(data.front());
		vector<size_t> forward;
		for (size_t i = 0; i < raw_forward.size(); ++i)
			if (minD - (int)max_distance_ <= (int)raw_forward[i] && (int)raw_forward[i] <= maxD + (int)max_distance_)
				forward.push_back(raw_forward[i]);
		if (forward.size() == 0)
			return result;
		
        size_t cur_dist = 0;
		vector<double> weights(forward.size());
		for (size_t i = 0; i < data.size(); i++) {
            if (math::ls(2 * data[i].d + second_len, (double) first_len))
                continue;
			while (cur_dist + 1 < forward.size()
					&& forward[cur_dist + 1] < data[i].d) {
				cur_dist++;
			}
			if (cur_dist + 1 < forward.size()
					&& math::ls(forward[cur_dist + 1] - data[i].d,
							data[i].d - (int) forward[cur_dist])) {
				cur_dist++;
				if (std::abs(forward[cur_dist] - data[i].d) < max_distance_)
					weights[cur_dist] += data[i].weight * weight_f_((int) forward[cur_dist] - data[i].d);
			} else if (cur_dist + 1 < forward.size()
					&& math::eq(forward[cur_dist + 1] - data[i].d,
							data[i].d - (int) forward[cur_dist])) {
				if (std::abs(forward[cur_dist] - data[i].d) < max_distance_)
					weights[cur_dist] += data[i].weight * 0.5 * weight_f_((int) forward[cur_dist] - data[i].d);
				cur_dist++;
				if (std::abs(forward[cur_dist] - data[i].d) < max_distance_)
					weights[cur_dist] += data[i].weight * 0.5 * weight_f_((int) forward[cur_dist] - data[i].d);
			} else {
				if (std::abs(forward[cur_dist] - data[i].d) < max_distance_)
					weights[cur_dist] += data[i].weight * weight_f_((int) forward[cur_dist] - data[i].d);
			}
		}
        
		for (size_t i = 0; i < forward.size(); i++) {
			if (math::gr(weights[i], 0.)) {
				result.push_back(make_pair(forward[i], weights[i]));
			}
		}
		return result;
	}


public:
	WeightedDistanceEstimator(const Graph &graph,
			const PairedInfoIndex<Graph>& histogram,
			const GraphDistanceFinder<Graph>& distance_finder, boost::function<double(int)> weight_f, 
			size_t linkage_distance, size_t max_distance) :
			base(graph, histogram, distance_finder, linkage_distance, max_distance), max_distance_(max_distance), weight_f_(weight_f) {
	}

};
    

}
#endif

//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef NAIVE_DISTANCE_ESTIMATION_HPP_
#define NAIVE_DISTANCE_ESTIMATION_HPP_

#include "xmath.h"
#include "paired_info.hpp"
#include "omni_utils.hpp"
#include "distance_estimation.hpp"

namespace omnigraph {

template<class Graph>
class NaiveDistanceEstimator: public DistanceEstimator<Graph> {
	typedef DistanceEstimator<Graph> base;
	typedef typename Graph::EdgeId EdgeId;

    boost::function<double(int)> weight_f_;

	virtual vector<pair<size_t, double>> EstimateEdgePairDistances(size_t first_len, size_t second_len,
			vector<PairInfo<EdgeId>> data,
			vector<size_t> raw_forward) const {
		vector<pair<size_t, double>> result;
		int maxD = rounded_d(data.back());
		int minD = rounded_d(data.front());
		vector<size_t> forward;
		for (size_t i = 0; i < raw_forward.size(); ++i)
			if (minD - (int)base::max_distance_ < (int)raw_forward[i] && (int)raw_forward[i] < maxD + (int)base::max_distance_)
				forward.push_back(raw_forward[i]);
		if (forward.size() == 0)
			return result;
		
		vector<double> weights(forward.size());
		for (size_t i = 0; i < data.size(); i++) {
			for (size_t dist = 0; dist < forward.size(); ++dist) {
                if (std::abs(forward[dist] - data[i].d) < base::max_distance_)
                    weights[dist] += data[i].weight * weight_f_(forward[dist] - data[i].d); 
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
	NaiveDistanceEstimator(const Graph &graph,
			const PairedInfoIndex<Graph>& histogram,
			const GraphDistanceFinder<Graph>& distance_finder, boost::function<double(int)> weight_f, 
			size_t linkage_distance, size_t max_distance) :
			base(graph, histogram, distance_finder, linkage_distance, max_distance), weight_f_(weight_f) {
	}

	virtual ~NaiveDistanceEstimator() {
	}

};
    

}
#endif

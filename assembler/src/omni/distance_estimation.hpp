#ifndef DISTANCE_ESTIMATION_HPP_
#define DISTANCE_ESTIMATION_HPP_

#include "paired_info.hpp"

namespace omnigraph {

template<class Graph>
class AbstractDistanceEstimator {
public:
	virtual ~AbstractDistanceEstimator() {
	}

	virtual void Estimate(PairedInfoIndex<Graph> &result) = 0;
};

template<class Graph>
class DistanceEstimator: AbstractDistanceEstimator<Graph> {
private:
	Graph &graph_;
	PairedInfoIndex<Graph> &histogram_;

	vector<size_t> GetGraphDistances(EdgeId first, EdgeId second) {

	}

	void EstimateEdgePairDistances(PairedInfoIndex<Graph> &result,
			vector<PairInfo<EdgeId> > data, vector<size_t> forward) {

	}

public:
	DistanceEstimator(Graph &graph, PairedInfoIndex<Graph> &histogram) :
		graph_(graph), histogram_(histogram) {
	}

	virtual ~DistanceEstimator() {
	}

	virtual void Estimate(PairedInfoIndex<Graph> &result) {
		for (auto iterator = histogram_.begin(); iterator != histogram_.end(); ++iterator) {
			vector < PairInfo<EdgeId> > data = *iterator;
			EdgeId first = data[0].first;
			EdgeId second = data[0].second;
			vector < size_t > forward = GetGraphDistances(first, second);
			EstimateEdgePairDistances(result, data, forward);
		}
	}
};

}

#endif /* DISTANCE_ESTIMATION_HPP_ */

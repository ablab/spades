#ifndef DISTANCE_ESTIMATION_HPP_
#define DISTANCE_ESTIMATION_HPP_

#include "paired_info.hpp"
#include "omni_utils.hpp"

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
	size_t insert_size_;
	size_t read_length_;
	size_t gap_;
	size_t delta_;

	const vector<size_t> GetGraphDistances(EdgeId first, EdgeId second) {
		DifferentDistancesCallback<Graph> callback(graph_);
		PathProcessor<Graph> path_processor(
				graph_,
				0. + gap_ - delta_ - graph_.length(first) - graph_.length(
						second), insert_size_ + delta_, graph_.EdgeEnd(first),
				graph_.EdgeStart(second), callback);
		path_processor.Process();
		return callback.distances();
	}

	void EstimateEdgePairDistances(PairedInfoIndex<Graph> &result,
			vector<PairInfo<EdgeId> > data, vector<size_t> forward) {
		sort(forward.begin(), forward.end());
		size_t cur = 0;
		for (size_t i = 0; i < forward.size(); i++) {
			double weight = 0;
			for (; cur < data.size(); cur++) {
				if (data[i].d < 0) {
					continue;
				}
				if (i + 1 < forward.size() && forward[i + 1]
						- data[cur].d < data[cur].d
						- forward[i]) {
					break;
				}
				weight += data[cur].weight;
			}
			if(weight > 0) {
				result.AddPairInfo(PairInfo(data[0].first, data[0].second, forward[i], weight));
			}
		}
	}

public:
	DistanceEstimator(Graph &graph, PairedInfoIndex<Graph> &histogram,
			size_t insert_size, size_t read_length, size_t delta) :
		graph_(graph), histogram_(histogram), insert_size_(insert_size),
				read_length_(read_length),
				gap_(insert_size - 2 * read_length_), delta_(delta) {
	}

	virtual ~DistanceEstimator() {
	}

	virtual void Estimate(PairedInfoIndex<Graph> &result) {
		for (auto iterator = histogram_.begin(); iterator != histogram_.end(); ++iterator) {
			vector<PairInfo<EdgeId> > data = *iterator;
			EdgeId first = data[0].first;
			EdgeId second = data[0].second;
			vector < size_t > forward = GetGraphDistances(first, second);
			EstimateEdgePairDistances(result, data, forward);
		}
	}
};

}

#endif /* DISTANCE_ESTIMATION_HPP_ */

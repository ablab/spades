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
protected:
	typedef typename Graph::EdgeId EdgeId;
	Graph &graph_;
	PairedInfoIndex<Graph> &histogram_;
	IdTrackHandler<Graph> &int_ids_;
	size_t insert_size_;
	size_t read_length_;
	size_t gap_;
	size_t delta_;
	size_t linkage_distance_;
	size_t max_distance_;

	const vector<size_t> GetGraphDistances(EdgeId first, EdgeId second) {
		DifferentDistancesCallback < Graph > callback(graph_);
		PathProcessor < Graph > path_processor(graph_,
				omnigraph::PairInfoPathLengthLowerBound(graph_.k(),
						graph_.length(first), graph_.length(second), gap_,
						delta_), omnigraph::PairInfoPathLengthUpperBound(
						graph_.k(), insert_size_, delta_)
				//				0. + gap_  + 2 * (graph_.k() + 1) - graph_.k() - delta_ - graph_.length(first) - graph_.length(
				//						second), insert_size_ - graph_.k() - 2 + delta_
				, graph_.EdgeEnd(first), graph_.EdgeStart(second), callback);
		path_processor.Process();
		auto result = callback.distances();
		for (size_t i = 0; i < result.size(); i++) {
			result[i] += graph_.length(first);
		}
		if (first == second) {
			result.push_back(0);
		}
		sort(result.begin(), result.end());
		return result;
	}

	vector<pair<size_t, double> > EstimateEdgePairDistances(vector<PairInfo<
			EdgeId> > data, vector<size_t> forward, bool debug = false) {
		vector < pair<size_t, double> > result;
		if (forward.size() == 0)
			return result;
		//        if (debug) for (size_t i = 0; i<forward.size(); i++) INFO("Distances " << forward[i]);
		size_t cur_dist = 0;
		vector<double> weights(forward.size());
		for (size_t i = 0; i < data.size(); i++) {
			if (data[i].d < 0)
				continue;
			while (cur_dist + 1 < forward.size() && forward[cur_dist + 1]
					< data[i].d) {
				cur_dist++;
			}
			if (cur_dist + 1 < forward.size() && math::ls(forward[cur_dist + 1]
					- data[i].d, data[i].d - forward[cur_dist])) {
				cur_dist++;
				if (std::abs(forward[cur_dist] - data[i].d) < max_distance_)
					weights[cur_dist] += data[i].weight;
				//                if (debug) INFO("Adding " << forward[cur_dist] << " " << data[i].d << " " << data[i].weight);
			} else if (cur_dist + 1 < forward.size() && math::eq(
					forward[cur_dist + 1] - data[i].d, data[i].d
							- forward[cur_dist])) {
				if (std::abs(forward[cur_dist] - data[i].d) < max_distance_)
					weights[cur_dist] += data[i].weight * 0.5;
				//                if (debug) INFO("Adding " << forward[cur_dist] << " " << data[i].d << " " << data[i].weight * 0.5);
				cur_dist++;
				if (std::abs(forward[cur_dist] - data[i].d) < max_distance_)
					weights[cur_dist] += data[i].weight * 0.5;
				//                if (debug) INFO("Adding " << forward[cur_dist] << " " << data[i].d << " " << data[i].weight * 0.5);
			} else {
				if (std::abs(forward[cur_dist] - data[i].d) < max_distance_)
					weights[cur_dist] += data[i].weight;
			}
		}
		for (size_t i = 0; i < forward.size(); i++) {
			if (weights[i] != 0) {
				result.push_back(make_pair(forward[i], weights[i]));
			}
		}
		return result;
	}

	vector<PairInfo<EdgeId> > ClusterResult(EdgeId edge1, EdgeId edge2, vector<
			pair<size_t, double> > estimated) {
		vector < PairInfo < EdgeId >> result;
		for (size_t i = 0; i < estimated.size(); i++) {
			size_t left = i;
			double weight = estimated[i].second;
			while (i + 1 < estimated.size() && estimated[i + 1].first
					- estimated[i].first <= linkage_distance_) {
				i++;
				weight += estimated[i].second;
			}
			double center = (estimated[left].first + estimated[i].first) * 0.5;
			double var = (estimated[i].first - estimated[left].first) * 0.5;
			PairInfo < EdgeId > new_info(edge1, edge2, center, weight, var);
			result.push_back(new_info);
		}
		return result;
	}

	void AddToResult(PairedInfoIndex<Graph> &result,
			vector<PairInfo<EdgeId> > clustered) {
		for (auto it = clustered.begin(); it != clustered.end(); ++it) {
			result.AddPairInfo(*it);
		}
	}

public:
	DistanceEstimator(Graph &graph, PairedInfoIndex<Graph> &histogram,
			IdTrackHandler<Graph> &int_ids, size_t insert_size,
			size_t read_length, size_t delta, size_t linkage_distance,
			size_t max_distance) :
		graph_(graph), histogram_(histogram), int_ids_(int_ids), insert_size_(
				insert_size), read_length_(read_length), gap_(insert_size - 2
				* read_length_), delta_(delta), linkage_distance_(
				linkage_distance), max_distance_(max_distance) {
	}

	virtual ~DistanceEstimator() {
	}

	virtual void Estimate(PairedInfoIndex<Graph> &result) {
		for (auto iterator = histogram_.begin(); iterator != histogram_.end(); ++iterator) {
			vector < PairInfo<EdgeId> > data = *iterator;
			EdgeId first = data[0].first;
			EdgeId second = data[0].second;
			vector < size_t > forward = GetGraphDistances(first, second);
			//			bool debug = (int_ids_.ReturnIntId(data[0].first) == 71456 && int_ids_.ReturnIntId(data[0].second) == 71195);

			vector < pair<size_t, double> > estimated
					= EstimateEdgePairDistances(data, forward, false);
			//            if (debug) for (size_t i = 0; i< estimated.size(); i++)
			//                INFO("Edges MY : " << estimated[i].first << " " << estimated[i].second);
			vector < PairInfo<EdgeId> > clustered = ClusterResult(first,
					second, estimated);
			//            if (debug) for (size_t i = 0; i<clustered.size(); i++)
			//                INFO("Edges MY clusterizing: " << clustered[i].d << " " << clustered[i].weight);
			AddToResult(result, clustered);
		}
	}
};

}

#endif /* DISTANCE_ESTIMATION_HPP_ */

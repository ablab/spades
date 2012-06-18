//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef TRICKY_DISTANCE_ESTIMATION_HPP_
#define TRICKY_DISTANCE_ESTIMATION_HPP_

#include "xmath.h"
#include "paired_info.hpp"
#include "omni_utils.hpp"
#include "distance_estimation.hpp"

namespace omnigraph {

template<class Graph>
class TrickyDistanceEstimator: AbstractDistanceEstimator<Graph> {
	typedef AbstractDistanceEstimator<Graph> base;
	typedef typename Graph::EdgeId EdgeId;

	const size_t max_distance_;

    boost::function<size_t(int)> weight_f_;

	vector<pair<size_t, double>> EstimateEdgePairDistances(size_t first_len, size_t second_len,
			vector<PairInfo<EdgeId>> data,
			vector<size_t> raw_forward) {
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
			if (2 * data[i].d + second_len < first_len)
				continue;
			while (cur_dist + 1 < forward.size()
					&& forward[cur_dist + 1] < data[i].d) {
				cur_dist++;
			}
			if (cur_dist + 1 < forward.size()
					&& math::ls(forward[cur_dist + 1] - data[i].d,
							data[i].d - forward[cur_dist])) {
				cur_dist++;
				if (std::abs(forward[cur_dist] - data[i].d) < max_distance_)
					weights[cur_dist] += data[i].weight * weight_f_(forward[cur_dist] - data[i].d);
			} else if (cur_dist + 1 < forward.size()
					&& math::eq(forward[cur_dist + 1] - data[i].d,
							data[i].d - forward[cur_dist])) {
				if (std::abs(forward[cur_dist] - data[i].d) < max_distance_)
					weights[cur_dist] += data[i].weight * 0.5 * weight_f_(forward[cur_dist] - data[i].d); 
				cur_dist++;
				if (std::abs(forward[cur_dist] - data[i].d) < max_distance_)
					weights[cur_dist] += data[i].weight * 0.5 * weight_f_(forward[cur_dist] - data[i].d);
			} else {
				if (std::abs(forward[cur_dist] - data[i].d) < max_distance_)
					weights[cur_dist] += data[i].weight * weight_f_(forward[cur_dist] - data[i].d);
			}
		}
		for (size_t i = 0; i < forward.size(); i++) {
			if (weights[i] != 0) {
				result.push_back(make_pair(forward[i], weights[i]));
			}
		}
		return result;
	}

	pair<EdgeId, EdgeId> ConjugatePair(EdgeId first, EdgeId second) {
		return make_pair(this->graph().conjugate(second), this->graph().conjugate(first));
	}

	PairInfo<EdgeId> ConjugateInfo(const PairInfo<EdgeId>& pair_info) {
		return PairInfo<EdgeId>(
				this->graph().conjugate(pair_info.second),
				this->graph().conjugate(pair_info.first),
				pair_info.d + this->graph().length(pair_info.second)
						- this->graph().length(pair_info.first),
				pair_info.weight, pair_info.variance);
	}

	vector<PairInfo<EdgeId>> ConjugateInfos(const vector<PairInfo<EdgeId>>& data) {
		vector<PairInfo<EdgeId>> answer;
		for (auto it = data.begin(); it != data.end(); ++it) {
			answer.push_back(ConjugateInfo(*it));
		}
		return answer;
	}

	void ProcessEdgePair(EdgeId first, EdgeId second, const vector<PairInfo<EdgeId>>& data, PairedInfoIndex<Graph> &result) {
		if (make_pair(first, second) <= ConjugatePair(first, second)) {
			vector<size_t> forward = this->GetGraphDistances(first, second);
			vector<pair<size_t, double> > estimated = EstimateEdgePairDistances(this->graph().length(first), this->graph().length(second),
				data, forward/*, false*/);
			vector<PairInfo<EdgeId>> res = ClusterResult(first, second, estimated);
			this->AddToResult(result, res);
			this->AddToResult(result, ConjugateInfos(res));
		}
	}

public:
	TrickyDistanceEstimator(const Graph &graph,
			const PairedInfoIndex<Graph>& histogram,
			const GraphDistanceFinder<Graph>& distance_finder, boost::function<size_t(int)> weight_f, 
			size_t linkage_distance, size_t max_distance) :
			base(graph, histogram, distance_finder, linkage_distance), max_distance_(
					max_distance), weight_f_(weight_f) {
	}

	virtual ~TrickyDistanceEstimator() {
	}

	virtual void Estimate(PairedInfoIndex<Graph> &result) {
		for (auto it = this->histogram().begin();
				it != this->histogram().end(); ++it) {
			ProcessEdgePair(it.first(), it.second(), *it, result);
		}
	}

    virtual void EstimateParallel(PairedInfoIndex<Graph> &result, size_t nthreads) {
        std::vector< std::pair<EdgeId, EdgeId> > edge_pairs;

        INFO("Collecting edge pairs");

        for (auto iterator = this->histogram().begin();
                iterator != this->histogram().end(); ++iterator) {

            edge_pairs.push_back(std::make_pair(iterator.first(), iterator.second()));
        }

        std::vector< PairedInfoIndex<Graph>* > buffer(nthreads);
        buffer[0] = &result;
        for (size_t i = 1; i < nthreads; ++i) {
            buffer[i] = new PairedInfoIndex<Graph>(this->graph(), result.GetMaxDifference());
        }

        INFO("Processing");
        #pragma omp parallel num_threads(nthreads)
        {
            #pragma omp for
            for (size_t i = 0; i < edge_pairs.size(); ++i)
            {
                EdgeId first = edge_pairs[i].first;
                EdgeId second = edge_pairs[i].second;
                ProcessEdgePair(first, second, this->histogram().GetEdgePairInfo(first, second), *buffer[omp_get_thread_num()]);
            }
        }

        INFO("Merging maps");
        for (size_t i = 1; i < nthreads; ++i) {
            buffer[0]->AddAll(*(buffer[i]));
            delete buffer[i];
        }
    }
};
    

}
#endif

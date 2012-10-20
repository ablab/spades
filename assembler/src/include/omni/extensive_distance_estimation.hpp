//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef EXTENSIVE_DISTANCE_ESTIMATION_HPP_
#define EXTENSIVE_DISTANCE_ESTIMATION_HPP_

#include "xmath.h"
#include "paired_info.hpp"
#include "omni_utils.hpp"
#include "distance_estimation.hpp"
#include <algorithm>

// No variation support in the original data

namespace omnigraph {

template<class Graph>
class ExtensiveDistanceEstimator: public WeightedDistanceEstimator<Graph> {
public:
	ExtensiveDistanceEstimator(const Graph &graph,
			const PairedInfoIndex<Graph>& histogram,
			const GraphDistanceFinder<Graph>& distance_finder, boost::function<double(int)> weight_f, 
            size_t linkage_distance, size_t max_distance) :
			base(graph, histogram, distance_finder, weight_f, linkage_distance, max_distance) {
	}

	virtual ~ExtensiveDistanceEstimator() {
	}

protected:
    typedef WeightedDistanceEstimator<Graph> base;
	typedef typename base::EdgeId EdgeId;

    void ExtendInfoLeft(EdgeId first, EdgeId second, vector<PairInfo<EdgeId>>& data, size_t max_shift) const {
        ExtendLeftDFS(first, second, data, 0, max_shift);
    }

    void ExtendInfoRight(EdgeId first, EdgeId second, vector<PairInfo<EdgeId>>& data, size_t max_shift) const {
        ExtendRightDFS(first, second, data, 0, max_shift);
    }

private:
	typedef typename Graph::VertexId VertexId;

	virtual void ProcessEdgePair(EdgeId first, EdgeId second, const vector<PairInfo<EdgeId>>& raw_data, PairedInfoIndex<Graph> &result) const {
		if (make_pair(first, second) <= this->ConjugatePair(first, second)) {
			const vector<size_t>& forward = this->GetGraphDistancesLengths(first, second);
            vector<PairInfo<EdgeId>> data(raw_data);
            DEBUG("Extending paired information");
            double weight_0 = WeightSum(data);
            DEBUG("Extend left");
            ExtendInfoLeft(first, second, data, 1000);
            DEBUG("Extend right");
            ExtendInfoRight(first, second, data, 1000);
            DEBUG("Weight increased " << (WeightSum(data) - weight_0));
            const vector<pair<int, double> >& estimated = this->EstimateEdgePairDistances(first, second, data, forward);
			const vector<PairInfo<EdgeId>>& res = this->ClusterResult(first, second, estimated);
			this->AddToResult(result, res);
			this->AddToResult(result, this->ConjugateInfos(res));
		}
	}

    double WeightSum(const vector<PairInfo<EdgeId>>& data) const {
        double answer = 0.;
        for (auto iter = data.begin(); iter != data.end(); ++iter) {
            answer += iter->weight;
        }
        return answer;
    }

    bool IsSorted(const vector<PairInfo<EdgeId>>& hist) const {
        for (size_t i = 0; i < hist.size() - 1; ++i) {
            if (hist[i].d > hist[i + 1].d) 
                return false;
        }
        return true;
    }

    void MergeInto(const vector<PairInfo<EdgeId>>& what, vector<PairInfo<EdgeId>>& where, int shift) const {
        // assuming they are sorted already

        if (what.size() == 0)
            return;
        if (where.size() == 0) {
            where = what;
            VERIFY(IsSorted(where));
            return;
        }
            
        // heuristics
        if (math::le(where.back().d,  what.front().d + shift)) {
            for (auto iter = what.begin(); iter != what.end(); ++iter) {
                PairInfo<EdgeId> to_be_added = *iter;
                to_be_added.first = where[0].first;
                to_be_added.second = where[0].second;
                to_be_added.d += shift;
                where.push_back(to_be_added);
            }
        }
        else {
            for (auto iter = what.begin(); iter != what.end(); ++iter) {
                PairInfo<EdgeId> to_be_added = *iter;
                to_be_added.first = where[0].first;
                to_be_added.second = where[0].second;
                to_be_added.d += shift;
                auto low_bound = lower_bound(where.begin(), where.end(), to_be_added);
                if (math::eq(to_be_added.d, low_bound->d)) {
                    low_bound->weight += to_be_added.weight;
                } 
                else {
                    where.insert(low_bound, to_be_added);
                }
            }
        }
        VERIFY(IsSorted(where));
    }

    vector<PairInfo<EdgeId>> FilterPositive(const vector<PairInfo<EdgeId>>& data, size_t first_len, size_t second_len) const {
        // assuming it is sorted
        if (data.size() == 0)
            return data;
        vector<PairInfo<EdgeId>> answer;
        for (auto iterator = data.begin(); iterator != data.end(); ++iterator) {
            if (math::ge(2. * iterator->d + second_len, (double) first_len))
                answer.push_back(*iterator);
        }
        return answer;
    }

    // left edge being extended to the left, shift is negative always
    void ExtendLeftDFS(EdgeId current, const EdgeId last, vector<PairInfo<EdgeId>>& data, int shift, size_t max_shift) const {
        VertexId start = this->graph().EdgeStart(current);
        if (current == last) 
            return;
        if (this->graph().OutgoingEdgeCount(start) > 1) 
            return; 
        const vector<EdgeId>& InEdges = this->graph().IncomingEdges(start);
        for (auto iterator = InEdges.begin(); iterator != InEdges.end(); ++iterator) {
            EdgeId next = *iterator;
            const vector<PairInfo<EdgeId>>& infos = this->histogram().GetEdgePairInfo(next, last);
            if (-shift < (int) max_shift)
                ExtendLeftDFS(next, last, data, shift - (int) this->graph().length(next), max_shift);
            const vector<PairInfo<EdgeId>>& filtered_infos = FilterPositive(infos, this->graph().length(next), this->graph().length(last));
            if (filtered_infos.size() > 0) 
                MergeInto(filtered_infos, data, shift - (int) this->graph().length(next));
        }
    }

    // right edge being extended to the right, shift is negative always
    void ExtendRightDFS(const EdgeId first, EdgeId current, vector<PairInfo<EdgeId>>& data, int shift, size_t max_shift) const {
        VertexId end = this->graph().EdgeEnd(current);
        if (current == first)
            return;
        if (this->graph().IncomingEdgeCount(end) > 1) 
            return;
        const vector<EdgeId>& OutEdges = this->graph().OutgoingEdges(end);
        for (auto iter = OutEdges.begin(); iter != OutEdges.end(); ++iter) {
            EdgeId next = *iter;
            const vector<PairInfo<EdgeId>>& infos = this->histogram().GetEdgePairInfo(first, next);
            if (-shift < (int) max_shift)
                ExtendRightDFS(first, next, data, shift - (int) this->graph().length(current), max_shift);
            const vector<PairInfo<EdgeId>>& filtered_infos = FilterPositive(infos, this->graph().length(first), this->graph().length(next));
            if (filtered_infos.size() > 0)
                MergeInto(filtered_infos, data, shift - (int) this->graph().length(current));
        }
    }

    DECL_LOGGER("ExtensiveDistanceEstimator")
};
    

}
#endif

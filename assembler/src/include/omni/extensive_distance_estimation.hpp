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
	typedef WeightedDistanceEstimator<Graph> base;
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

    //std::set<pair<EdgeId, EdgeId>> ExtendedLeft;

    //std::set<pair<EdgeId, EdgeId>> ExtendedRight;

    double WeightSum(const vector<PairInfo<EdgeId>>& data) const {
        double answer = 0.;
        for (auto iter = data.begin(); iter != data.end(); ++iter) {
            answer += iter->weight;
        }
        return answer;
    }
    bool isSorted(const vector<PairInfo<EdgeId>>& hist) const {
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
            VERIFY(isSorted(where));
            return;
        }
        //INFO("BEFORE SORTING");
        //for (auto iter = where.begin(); iter != where.end(); ++iter) {
            //cout << iter->d << " " << endl;   
        //}
            
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
                //INFO("Something is going to be added!");
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
        VERIFY(isSorted(where));

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
    void ExtendLeftDFS(EdgeId current, const EdgeId last, vector<PairInfo<EdgeId>>& data, int shift) const {
        VertexId start = this->graph().EdgeStart(current);
        
        if (current == last) 
            return;

        if (this->graph().OutgoingEdgeCount(start) > 1) 
            return; 

        const vector<EdgeId>& InEdges = this->graph().IncomingEdges(start);

        for (auto iterator = InEdges.begin(); iterator != InEdges.end(); ++iterator) {
            EdgeId next = *iterator;
            
            const vector<PairInfo<EdgeId>>& infos = this->histogram().GetEdgePairInfo(next, last);

            if (-shift < 10000)
                ExtendLeftDFS(next, last, data, shift - (int) this->graph().length(next));
 
            const vector<PairInfo<EdgeId>>& filtered_infos = FilterPositive(infos, this->graph().length(next), this->graph().length(last));

            if (filtered_infos.size() > 0) 
                MergeInto(filtered_infos, data, shift - (int) this->graph().length(next));
        }
    }

    // right edge being extended to the right, shift is negative always
    void ExtendRightDFS(const EdgeId first, EdgeId current, vector<PairInfo<EdgeId>>& data, int shift) const {
        VertexId end = this->graph().EdgeEnd(current);

        if (current == first)
            return;

        if (this->graph().IncomingEdgeCount(end) > 1) 
            return;
            
        const vector<EdgeId>& OutEdges = this->graph().OutgoingEdges(end);
        for (auto iter = OutEdges.begin(); iter != OutEdges.end(); ++iter) {
            EdgeId next = *iter;

            const vector<PairInfo<EdgeId>>& infos = this->histogram().GetEdgePairInfo(first, next);

            if (-shift < 10000)
                ExtendRightDFS(first, next, data, shift - (int) this->graph().length(current));

            const vector<PairInfo<EdgeId>>& filtered_infos = FilterPositive(infos, this->graph().length(first), this->graph().length(next));

            if (filtered_infos.size() > 0)
                MergeInto(filtered_infos, data, shift - (int) this->graph().length(current));
        }
    }


    void ExtendInfoLeft(const EdgeId first, const EdgeId second, vector<PairInfo<EdgeId>>& data) const {
        ExtendLeftDFS(first, second, data, 0);
    }

    void ExtendInfoRight(const EdgeId first, const EdgeId second, vector<PairInfo<EdgeId>>& data) const {
        ExtendRightDFS(first, second, data, 0);
    }

	void ProcessEdgePair(const EdgeId first, const EdgeId second, const vector<PairInfo<EdgeId>>& raw_data, PairedInfoIndex<Graph> &result) const {
		if (make_pair(first, second) <= this->ConjugatePair(first, second)) {
			vector<size_t> forward = this->GetGraphDistances(first, second);
            vector<PairInfo<EdgeId>> data = raw_data;
            DEBUG("Extending paired information");
            double weight_0 = WeightSum(data);
            DEBUG("Extend left");
            ExtendInfoLeft(first, second, data);
            DEBUG("Extend right");
            ExtendInfoRight(first, second, data);
            
            //INFO("Then negative");

            //vector<PairInfo<EdgeId>> tmp_data;
            //INFO("Extend left");
            //ExtendInfoLeft(second, first, tmp_data);
            //INFO("Extend right");
            //ExtendInfoRight(second, first, tmp_data);

            //INFO("EXTENDED");
            //vector<PairInfo<EdgeId>> reflected_data;
            //for (int i = tmp_data.size() - 1; i >= 0; --i) {
                ////int len_1 = this->graph().length(first);
                ////int len_2 = this->graph().length(second);
                //if (math::ls(tmp_data[i].d, 0.)) 
                    //break;
                //tmp_data[i].first = first;
                //tmp_data[i].second = second;
                //tmp_data[i].d = -(tmp_data[i].d); // - len_2 + len_1);
                //reflected_data.push_back(tmp_data[i]);
            //}

            //MergeInto(reflected_data, data, 0);

            DEBUG("Weight increased " << (WeightSum(data) - weight_0));

			vector<pair<size_t, double> > estimated = this->EstimateEdgePairDistances(first, second,
				data, forward);
			vector<PairInfo<EdgeId>> res = this->ClusterResult(first, second, estimated);
			this->AddToResult(result, res);
			this->AddToResult(result, this->ConjugateInfos(res));
		}
	}

public:
	ExtensiveDistanceEstimator(const Graph &graph,
			const PairedInfoIndex<Graph>& histogram,
			const GraphDistanceFinder<Graph>& distance_finder, boost::function<double(int)> weight_f, 
            size_t linkage_distance, size_t max_distance) :
			base(graph, histogram, distance_finder, weight_f, linkage_distance, max_distance) {
	}

	virtual ~ExtensiveDistanceEstimator() {
	}

private:
    DECL_LOGGER("ExtensiveDistanceEstimator")
};
    

}
#endif

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
                to_be_added.second = where[0].second;
                to_be_added.d += shift;
                where.push_back(to_be_added);
            }
        }
        else {
            for (auto iter = what.begin(); iter != what.end(); ++iter) {
                //INFO("Something is going to be added!");
                PairInfo<EdgeId> to_be_added = *iter;
                double dist_to_be_added = to_be_added.d + shift;

                to_be_added.second = where[0].second;
                to_be_added.d = dist_to_be_added;


                auto low_bound = lower_bound(where.begin(), where.end(), to_be_added);
                
                if (dist_to_be_added == low_bound->d) {
                    low_bound->weight += to_be_added.weight;
                } 
                else {
                    //INFO("INSERTING " << to_be_added);
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

    // right edge being extended to the left
    void ExtendLeftDFS(const EdgeId first, EdgeId current, vector<PairInfo<EdgeId>>& data, size_t length, const size_t second_len) const {
        VertexId start = this->graph().EdgeStart(current);
        
        int shift = ((int) length - (int) second_len);


        if (current == first) 
            return;

        if (this->graph().OutgoingEdgeCount(start) > 1) 
            return; 

        const vector<EdgeId>& InEdges = this->graph().IncomingEdges(start);

        for (auto iterator = InEdges.begin(); iterator != InEdges.end(); ++iterator) {
            EdgeId next = *iterator;
            
            const vector<PairInfo<EdgeId>>& infos = this->histogram().GetEdgePairInfo(first, next);

            if (shift < 1000)
                ExtendLeftDFS(first, next, data, length + this->graph().length(next), second_len);
 
            const vector<PairInfo<EdgeId>>& filtered_infos = FilterPositive(infos, this->graph().length(first), this->graph().length(next));

            if (filtered_infos.size() > 0) 
                MergeInto(filtered_infos, data, shift + this->graph().length(next));
        }
    }

    // right edge being extended to the right
    void ExtendRightDFS(const EdgeId first, EdgeId current, vector<PairInfo<EdgeId>>& data, size_t length, const size_t second_len) const {
        VertexId end = this->graph().EdgeEnd(current);
        
        int shift = -((int) length - (int) second_len);

        if (current == first)
            return;

        if (this->graph().IncomingEdgeCount(end) > 1) 
            return;
            
        const vector<EdgeId>& OutEdges = this->graph().OutgoingEdges(end);
        for (auto iter = OutEdges.begin(); iter != OutEdges.end(); ++iter) {
            EdgeId next = *iter;

            const vector<PairInfo<EdgeId>>& infos = this->histogram().GetEdgePairInfo(first, next);

            if ( (-shift) < 1000)
                ExtendRightDFS(first, next, data, length + this->graph().length(next), second_len);

            const vector<PairInfo<EdgeId>>& filtered_infos = FilterPositive(infos, this->graph().length(first), this->graph().length(next));

            if (filtered_infos.size() > 0)
                MergeInto(filtered_infos, data, shift - this->graph().length(next));
        }
    }


    void ExtendInfoLeft(const EdgeId first, const EdgeId second, vector<PairInfo<EdgeId>>& data) const {
        size_t second_len = this->graph().length(second);
        ExtendLeftDFS(first, second, data, second_len, second_len);
    }

    void ExtendInfoRight(const EdgeId first, const EdgeId second, vector<PairInfo<EdgeId>>& data) const {
        size_t second_len = this->graph().length(second);
        ExtendRightDFS(first, second, data, second_len, second_len);
    }

	void ProcessEdgePair(const EdgeId first, const EdgeId second, const vector<PairInfo<EdgeId>>& raw_data, PairedInfoIndex<Graph> &result) const {
		//if (make_pair(first, second) <= ConjugatePair(first, second)) {
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

			vector<pair<size_t, double> > estimated = EstimateEdgePairDistances(first, second,
				data, forward);
			vector<PairInfo<EdgeId>> res = ClusterResult(first, second, estimated);
			this->AddToResult(result, res);
			this->AddToResult(result, ConjugateInfos(res));
		//}
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

	virtual void Estimate(PairedInfoIndex<Graph> &result) const {
		for (auto it = this->histogram().begin();
				it != this->histogram().end(); ++it) {
			ProcessEdgePair(it.first(), it.second(), *it, result);
		}
	}

    virtual void EstimateParallel(PairedInfoIndex<Graph> &result, size_t nthreads) const {
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

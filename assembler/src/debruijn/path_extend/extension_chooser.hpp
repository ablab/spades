//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * extension.hpp
 *
 *  Created on: Mar 5, 2012
 *      Author: andrey
 */

#ifndef EXTENSION_HPP_
#define EXTENSION_HPP_

#include "weight_counter.hpp"
#include "pe_utils.hpp"
#include <iostream>
#include <fstream>

namespace path_extend {

typedef std::multimap<double, EdgeWithDistance> AlternativeConteiner;


class PathAnalyzer {

protected:
    Graph& g_;

public:
    PathAnalyzer(Graph& g): g_(g) {
    }

    int ExcludeTrivial(const BidirectionalPath& path, std::set<int>& edges, int from = -1) {
        int edgeIndex = (from == -1) ? path.Size() - 1 : from;
        if ((int) path.Size() <= from) {
            return edgeIndex;
        }
        VertexId currentVertex = g_.EdgeEnd(path[edgeIndex]);
        while (edgeIndex >= 0 && g_.CheckUniqueIncomingEdge(currentVertex)) {
            EdgeId e = g_.GetUniqueIncomingEdge(currentVertex);
            currentVertex = g_.EdgeStart(e);

            edges.insert(edgeIndex);
            --edgeIndex;
        }
        return edgeIndex;
    }

    int ExcludeTrivialWithBulges(const BidirectionalPath& path, std::set<int>& edges) {
        edges.clear();

        if (path.Empty()) {
            return 0;
        }

        int lastEdge = path.Size() - 1;
        do {
            lastEdge = ExcludeTrivial(path, edges, lastEdge);

            if (lastEdge >= 0) {
                VertexId v = g_.EdgeEnd(path[lastEdge]);
                VertexId u = g_.EdgeStart(path[lastEdge]);
                auto bulgeCandidates = g_.IncomingEdges(v);
                bool bulge = true;

                for (auto iter = bulgeCandidates.begin(); iter != bulgeCandidates.end(); ++iter) {
                    if (g_.EdgeStart(*iter) != u) {
                        bulge = false;
                        break;
                    }
                }

                if (!bulge) {
                    break;
                }

                --lastEdge;
            }
        } while (lastEdge >= 0);

        return lastEdge;
    }
};


class ExtensionChooserListener {

public:

    virtual void ExtensionChosen(double weight) = 0;

    virtual void ExtensionChosen(AlternativeConteiner& alts) = 0;

    virtual ~ExtensionChooserListener() {

    }
};


class ExtensionChooser {

public:
    typedef std::vector<EdgeWithDistance> EdgeContainer;

protected:
    Graph& g_;

    WeightCounter * wc_;

    PathAnalyzer analyzer_;

    double priorityCoefficient_;

    bool excludeTrivial_;
    bool excludeTrivialWithBulges_;

    std::vector<ExtensionChooserListener *> listeners_;

public:
    ExtensionChooser(Graph& g, WeightCounter * wc = 0, double priority = 0.0): g_(g), wc_(wc), analyzer_(g), priorityCoefficient_(priority),
        excludeTrivial_(true), excludeTrivialWithBulges_(true), listeners_() {
    }

    virtual ~ExtensionChooser() {

    }

    virtual EdgeContainer Filter(BidirectionalPath& path, EdgeContainer& edges) = 0;

    bool isExcludeTrivial() const
    {
        return excludeTrivial_;
    }

    double CountWeight(BidirectionalPath& path, EdgeId e) {
        return wc_->CountWeight(path, e);
    }

    bool isExcludeTrivialWithBulges() const
    {
        return excludeTrivialWithBulges_;
    }

    void setExcludeTrivial(bool excludeTrivial)
    {
        this->excludeTrivial_ = excludeTrivial;
    }

    void setExcludeTrivialWithBulges(bool excludeTrivialWithBulges)
    {
        this->excludeTrivialWithBulges_ = excludeTrivialWithBulges;
    }

    PairedInfoLibraries& getLibs() {
        return wc_->getLibs();
    }

    void Subscribe(ExtensionChooserListener * listener) {
        listeners_.push_back(listener);
    }

    void NotifyAll(double weight) {
        for (auto iter = listeners_.begin(); iter != listeners_.end(); ++iter) {
            (*iter)->ExtensionChosen(weight);
        }
    }

    void NotifyAll(AlternativeConteiner& alts) {
        for (auto iter = listeners_.begin(); iter != listeners_.end(); ++iter) {
            (*iter)->ExtensionChosen(alts);
        }
    }

    bool WeighConterBased() const {
        return wc_ != 0;
    }
	protected:
	DECL_LOGGER("ExtensionChooser");
};


class JointExtensionChooser: public ExtensionChooser {

protected:
    ExtensionChooser * first_;

    ExtensionChooser * second_;

public:
    JointExtensionChooser(Graph& g, ExtensionChooser * first, ExtensionChooser * second): ExtensionChooser(g),
        first_(first), second_(second)
    {
    }

    virtual EdgeContainer Filter(BidirectionalPath& path, EdgeContainer& edges) {
        EdgeContainer e1 = first_->Filter(path, edges);
        return second_->Filter(path, e1);
    }
};


class TrivialExtensionChooser: public ExtensionChooser {

public:
    TrivialExtensionChooser(Graph& g): ExtensionChooser(g)  {
    }

    virtual EdgeContainer Filter(BidirectionalPath& path, EdgeContainer& edges) {
        if (edges.size() == 1) {
             return edges;
        }
        return EdgeContainer();
    }
};


class TrivialExtensionChooserWithPI: public ExtensionChooser {

public:
    TrivialExtensionChooserWithPI(Graph& g, WeightCounter * wc): ExtensionChooser(g, wc) {
    }

    virtual EdgeContainer Filter(BidirectionalPath& path, EdgeContainer& edges) {
        wc_->GetExcludedEdges().clear();
        if (edges.size() == 1) {
                double weight = wc_->CountWeight(path, edges.back().e_);
                NotifyAll(weight);

                if (wc_->IsExtensionPossible(weight)) {
                    return edges;
                }
        }
        return EdgeContainer();
    }
};


class SimpleExtensionChooser: public ExtensionChooser {

protected:

    void RemoveTrivial(BidirectionalPath& path) {
        wc_->GetExcludedEdges().clear();

        if (excludeTrivialWithBulges_) {
            analyzer_.ExcludeTrivialWithBulges(path, wc_->GetExcludedEdges());
        }
        else if (excludeTrivial_) {
            analyzer_.ExcludeTrivial(path, wc_->GetExcludedEdges());
        }
    }

	void RemoveTrivialAndCommon(BidirectionalPath& path, EdgeId first,
			EdgeId second) {
		RemoveTrivial(path);
		if (path.Size() == 0) {
			return;
		}
		int index = path.Size() - 1;
		while (index >= 0) {
			bool common_edge = wc_->PairInfoExist(path[index], first,
					path.LengthAt(index))
					and wc_->PairInfoExist(path[index], second,
							path.LengthAt(index));
			bool ideal1 = wc_->CountIdealInfo(path[index], first,
					path.LengthAt(index)) > 0.0;
			bool ideal2 = wc_->CountIdealInfo(path[index], second,
					path.LengthAt(index)) > 0.0;
			if (common_edge or ideal1 != ideal2) {
				wc_->GetExcludedEdges().insert(index);
				DEBUG("excluded trivial and common " << index);
			}
			index--;

		}
	}

	void find_weights(BidirectionalPath& path, EdgeContainer& edges, AlternativeConteiner& weights) {
		for (auto iter = edges.begin(); iter != edges.end(); ++iter) {
			double weight = wc_->CountWeight(path, iter->e_);
			weights.insert(std::make_pair(weight, *iter));
			DEBUG("Candidate " << g_.int_id(iter->e_) << " weight " << weight);
			path.getLoopDetector().AddAlternative(iter->e_, weight);

		}
		NotifyAll(weights);
	}

	void find_possible_edges(AlternativeConteiner& weights, EdgeContainer& top, double maxWeight) {
		auto possibleEdge = weights.lower_bound(maxWeight / priorityCoefficient_);
		for (auto iter = possibleEdge; iter != weights.end(); ++iter) {
			top.push_back(iter->second);
		}
	}

	EdgeContainer find_result(BidirectionalPath& path, EdgeContainer& edges) {
		AlternativeConteiner weights;
		find_weights(path, edges, weights);
		EdgeContainer top;
		auto maxWeight = (--weights.end())->first;
		find_possible_edges(weights, top, maxWeight);
		EdgeContainer result;
		if (top.size() >= 1 && wc_->IsExtensionPossible(maxWeight)) {
			result = top;
		}
		return result;
	}
public:
    SimpleExtensionChooser(Graph& g, WeightCounter * wc, double priority): ExtensionChooser(g, wc, priority){

    }

    virtual EdgeContainer Filter(BidirectionalPath& path, EdgeContainer& edges) {
        if (edges.empty()) {
            return edges;
        }
        RemoveTrivial(path);
        path.Print();
        EdgeContainer result = find_result(path, edges);
        if (result.size() > 1){
        	DEBUG("result size MORE 1");
        	EdgeId first = result.at(0).e_;
        	EdgeId second = result.at(1).e_;
        	RemoveTrivialAndCommon(path, first, second);
        	result = find_result(path, edges);
        	if (result.size() == 1){
        		DEBUG("CHANGE RESULT");
        	}
        }
        return result;
    }

};



class ScaffoldingExtensionChooser: public ExtensionChooser {

    static bool compare(pair<int,double> a, pair<int,double> b)
    {
        if (a.first < b.first) return true;
        else return false;
    }

public:
	ScaffoldingExtensionChooser(Graph& g, WeightCounter * wc, double priority): ExtensionChooser(g, wc, priority) {

    }

	double AddInfoFromEdge(const std::vector<int>& distances, const std::vector<double>& weights, std::vector<pair<int,double> >& histogram, const BidirectionalPath& path, size_t j, double threshold)
	{
		double mean = 0.0;
		double sum  = 0.0;

		for (size_t l = 0; l < distances.size(); ++ l){
			if (distances[l] > max(0, (int) path.LengthAt(j) - (int) g_.k()) && weights[l] >= threshold) {
                mean += ((distances[l] - (int) path.LengthAt(j)) * weights[l]);
                sum += weights[l];
                histogram.push_back(make_pair(distances[l] - path.LengthAt(j), weights[l]));
			}
		}
		return mean / sum;
	}


    int CountMean(vector< pair<int,double> >& histogram)
    {
        double dist = 0.0;
        double sum = 0;
        for (size_t i = 0; i < histogram.size(); ++ i) {
			 sum += histogram[i].second;
		}
        for (size_t i = 0; i < histogram.size(); ++ i) {
			 dist += (histogram[i].first * histogram[i].second / sum);
		}
        return round(dist);
    }

    int CountDev(vector< pair<int,double> >& histogram, int avg)
    {
        double dev = 0.0;
        double sum = 0;
        for (size_t i = 0; i < histogram.size(); ++ i) {
             sum += histogram[i].second;
        }
        for (size_t i = 0; i < histogram.size(); ++ i) {
             dev += (((double) (histogram[i].first - avg)) * ((double) (histogram[i].first - avg)) * ((double) histogram[i].second));
        }
        return round(sqrt(dev / sum));
    }

    vector< pair<int,double> > FilterHistogram(vector< pair<int,double> >& histogram, int start, int end, int threshold)
    {
        vector< pair<int,double> > res;

        for (size_t i = 0; i < histogram.size(); ++i){
            if (histogram[i].first >= start && histogram[i].first <= end && histogram[i].second >= threshold) {
                res.push_back(histogram[i]);
            }
        }

        return res;
    }

    double CountAvrgDists(BidirectionalPath& path, EdgeId e, std::vector<pair<int,double> > & histogram)
    {
		std::vector<int> distances;
		std::vector<double> weights;

		double max_weight = 0.0;
		//bool print = true;
		for (size_t j = 0; j < path.Size(); ++ j) {
			wc_->GetDistances(path.At(j), e, distances, weights);

			for (size_t l = 0; l < weights.size(); ++ l){
				if (weights[l] > max_weight) {
				    max_weight = weights[l];
				}
			}

			if (distances.size() > 0) {
				AddInfoFromEdge(distances, weights, histogram, path, j, 0);
			}
			distances.clear();
		}
		return max_weight;
    }

    void FindBestFittedEdges(BidirectionalPath& path, EdgeContainer& edges, EdgeContainer& result)
    {
		std::vector<pair<int,double> > histogram;
		for (size_t i = 0; i < edges.size(); ++i){
			histogram.clear();
			double max_w = CountAvrgDists(path, edges[i].e_, histogram);

			for (int j = 0; j < 2; ++j) {
                int mean = CountMean(histogram);
                int dev = CountDev(histogram, mean);
                double cutoff = min(max_w * cfg::get().pe_params.param_set.scaffolder_options.rel_cutoff, (double) cfg::get().pe_params.param_set.scaffolder_options.cutoff);
                histogram = FilterHistogram(histogram, mean - (5 - j)  * dev, mean + (5 - j) * dev, cutoff);
			}

			double sum = 0.0;
			for (size_t j = 0; j < histogram.size(); ++j) {
			    sum += histogram[j].second;
			}

			if (sum > cfg::get().pe_params.param_set.scaffolder_options.sum_threshold) {
				sort(histogram.begin(), histogram.end(), compare);
				int gap = CountMean(histogram);

				if (wc_->CountIdealInfo(path, edges[i].e_, gap) > 0.0) {
					result.push_back(EdgeWithDistance(edges[i].e_, gap));
                }
			}
		}
    }


    void FindBestFittedEdgesForClustered(BidirectionalPath& path, EdgeContainer& edges, EdgeContainer& result)
    {
        std::vector<pair<int,double> > histogram;
        for (size_t i = 0; i < edges.size(); ++i){
            histogram.clear();
            CountAvrgDists(path, edges[i].e_, histogram);

            double sum = 0.0;
            for (size_t j = 0; j < histogram.size(); ++j) {
                sum += histogram[j].second;
            }

            if (sum > cfg::get().pe_params.param_set.scaffolder_options.cl_threshold) {
                sort(histogram.begin(), histogram.end(), compare);
                int gap = CountMean(histogram);

                if (wc_->CountIdealInfo(path, edges[i].e_, gap) > 0.0) {
                    result.push_back(EdgeWithDistance(edges[i].e_, gap));
                }
            }
        }
    }

    virtual EdgeContainer Filter(BidirectionalPath& path, EdgeContainer& edges) {
        if (edges.empty()) {
            return edges;
        }
        EdgeContainer result;

        if (cfg::get().pe_params.param_set.scaffolder_options.cluster_info) {
            FindBestFittedEdgesForClustered(path, edges, result);
        } else {
            FindBestFittedEdges(path, edges, result);
        }

        return result;
    }
};


bool ReverseComparePairBySecond(const boost::tuple<EdgeId, int,int> & a, const boost::tuple<EdgeId, int, int>& b) {
    return get<1>(a) > get<1>(b);
}

//
//class PathsDrivenExtensionChooser: public ExtensionChooser {
//
//protected:
//
//    GraphCoverageMap coverageMap_;
//
//public:
//    PathsDrivenExtensionChooser(Graph& g, PathContainer& pc): ExtensionChooser(g, 0, .0), coverageMap_(g, pc) {
//
//    }
//
//    virtual EdgeContainer Filter(BidirectionalPath& path, EdgeContainer& edges) {
//        if (edges.empty()) {
//            return edges;
//        }
//        DEBUG("We in Filter of PathsDrivenExtension");
//        set<EdgeId> candidatesSet;
//        for (auto it = edges.begin(); it != edges.end(); ++it) {
//            candidatesSet.insert(it->e_);
//        }
//
//        auto supportingPaths = coverageMap_.GetCoveringPaths(path.Back());
//        vector<pair<BidirectionalPath*, size_t> > supportPathCandidates;
//        set<EdgeId> filteredCandidatesSet;
//        EdgeId nullEdge;
//        for (auto it = supportingPaths.begin(); it != supportingPaths.end(); ++it) {
//            auto positions = (*it)->FindAll(path.Back());
//            for (size_t i = 0; i < positions.size(); ++i) {
//                if (positions[i] < (*it)->Size() - 1 && candidatesSet.count((*it)->At(positions[i] + 1)) > 0) {
//                    supportPathCandidates.push_back(make_pair(*it, positions[i]));
//                    filteredCandidatesSet.insert((*it)->At(positions[i] + 1));
//                } else if (positions[i] == (*it)->Size() - 1){
//                	supportPathCandidates.push_back(make_pair(*it, positions[i]));
//                	filteredCandidatesSet.insert(nullEdge);
//
//                }
//            }
//        }
//
//        if (filteredCandidatesSet.size() > 1) {
//            DEBUG("Several extensions are supported, calculating scores now");
//
//            vector< boost::tuple<EdgeId, int, int > > trustedCandidates;
//            for (size_t i = 0; i < supportPathCandidates.size(); ++i) {
//                int coveredEdges = 0;
//                BidirectionalPath * supportingPath = supportPathCandidates[i].first;
//                int backPos =  supportPathCandidates[i].second;
//                DEBUG("Supporting path #" << i);
//                supportingPath->Print();
//
//                while ((int) path.Size() - 1 - coveredEdges >= 0 && backPos - coveredEdges  >= 0) {
//                    if (path[(int) path.Size() - 1 - coveredEdges] != supportingPath->At(backPos - coveredEdges)) {
//                        break;
//                    }
//                    ++coveredEdges;
//                }
//                if (backPos < (int)supportingPath->Size() - 1){
//                	int unCoveredEdges = std::min((int)path.Size() - 1 - coveredEdges, backPos - coveredEdges);
//                	trustedCandidates.push_back(boost::make_tuple(supportingPath->At(backPos + 1), coveredEdges, unCoveredEdges));
//                	DEBUG("This path supports " << g_.int_id(supportingPath->At(backPos + 1)) <<
//                        " by " << coveredEdges << " edges");
//                } else {
//                	trustedCandidates.push_back(boost::make_tuple(nullEdge, coveredEdges, 0));
//                	DEBUG("This path supports end "  <<
//                	                        " by " << coveredEdges << " edges");
//                }
//            }
//
//            sort(trustedCandidates.begin(), trustedCandidates.end(), ReverseComparePairBySecond);
//            filteredCandidatesSet.clear();
//            size_t i = 0;
//            while (i < trustedCandidates.size() &&
//                    get<1>(trustedCandidates[0]) == get<1>(trustedCandidates[i])) {
//                filteredCandidatesSet.insert(get<0>(trustedCandidates[i]));
//                ++i;
//            }
//            if (get<0>(trustedCandidates[0]).get() == 0){
//            	filteredCandidatesSet.clear();
//            	size_t index = 0;
//            	while (get<0>(trustedCandidates[index]).get() == 0 or get<2>(trustedCandidates[index]) >= 0){
//            		index++;
//            	}//TODO: if not one variant:unCovered =0 and not null should return two ore more version
//            	if (index < trustedCandidates.size()){
//            		int first_good_index = index;
//            		while (index < trustedCandidates.size() and (get<1>(trustedCandidates[index]) == get<1>(trustedCandidates[first_good_index]))){
//            			if (get<0>(trustedCandidates[index]).get() != 0 or get<2>(trustedCandidates[index]) < 0){
//            				filteredCandidatesSet.insert(get<0>(trustedCandidates[index]));
//            			}
//            			index++;
//            		}
//            		DEBUG("NEW filterefCandidate " << get<1>(trustedCandidates[index]) << " "<<get<2>(trustedCandidates[index]));
//            	}
//            	DEBUG("ONLY with END PATH, trustedCandidates "<< trustedCandidates.size());
//            }
//            DEBUG("Found " << filteredCandidatesSet.size() << " trusted extension(s), supported paths " << i << " best score " << get<1>(trustedCandidates[0]));
//        } else if (filteredCandidatesSet.size() == 1){
//            DEBUG("Only one extension is supported: " << g_.int_id(*(filteredCandidatesSet.begin())));
//        } else {
//            DEBUG("NO extensions is supported" );
//        }
//
//        EdgeContainer result;
//        for (auto it = edges.begin(); it != edges.end(); ++it) {
//            if (filteredCandidatesSet.count(it->e_) > 0) {
//                result.push_back(*it);
//            }
//        }
//        DEBUG("result size " << result.size());
//        return result;
//    }
//};

class LongReadsExtensionChooser: public ExtensionChooser {

protected:

    GraphCoverageMap coverageMap_;


public:
    LongReadsExtensionChooser(Graph& g, PathContainer& pc): ExtensionChooser(g, 0, .0), coverageMap_(g, pc) {

    }

    virtual EdgeContainer Filter(BidirectionalPath& path, EdgeContainer& edges) {
        if (edges.empty()) {
            return edges;
        }
        DEBUG("We in Filter of PathsDrivenExtension");
        map<EdgeId, double> weights_cands;
        for (auto it = edges.begin(); it != edges.end(); ++it) {
            weights_cands.insert(make_pair(it->e_, 0));
        }
        set<EdgeId> filtered_cands;
        auto supporting_paths = coverageMap_.GetCoveringPaths(path.Back());
        for (auto it = supporting_paths.begin(); it != supporting_paths.end(); ++it) {
        	auto positions = (*it)->FindAll(path.Back());
            for (size_t i = 0; i < positions.size(); ++i) {
            	if (positions[i] > 0 && positions[i] < (*it)->Size() - 1 && covered_path(path, **it, positions[i] ) && unique_back_path(**it, positions[i])){
            		EdgeId next = (*it)->At(positions[i] + 1);
            		weights_cands[next] = weights_cands[next] + (*it)->getWeight();
            		filtered_cands.insert(next);
            	}
            }
        }

        if (filtered_cands.size() > 1) {
        	vector<pair<EdgeId, double> > sorted_candidates = to_vector(weights_cands);
        	DEBUG("First extension is supported" <<g_.int_id(sorted_candidates[0].first) << " weight " << sorted_candidates[0].second);
        	DEBUG("First extension is supported" <<g_.int_id(sorted_candidates[1].first) << " weight " << sorted_candidates[1].second);
        	if (sorted_candidates[0].second > 1.5 * sorted_candidates[1].second && sorted_candidates[0].second > 1.5){
        		filtered_cands.clear();
        		filtered_cands.insert(sorted_candidates[0].first);
        	}
        } else if (filtered_cands.size() == 1 && weights_cands[ *(filtered_cands.begin())] >1){
        	EdgeId candidate = *(filtered_cands.begin());
            DEBUG("Only one extension is supported: " << g_.int_id(candidate) << " with weight " << weights_cands[candidate]);
        } else {
            DEBUG("NO extensions is supported" );
        }
        EdgeContainer result;
        for (auto it = edges.begin(); it != edges.end(); ++it) {
            if (filtered_cands.find(it->e_) != filtered_cands.end()) {
                result.push_back(*it);
            }
        }
        DEBUG("result size " << result.size());
        return result;
    }

private:
    bool unique_back_path( BidirectionalPath& cov_path, size_t pos) {
    	DEBUG("running unique_path");
    	int int_pos = pos;
    	while (int_pos >=0) {
    		DEBUG(int_pos <<" "<< g_.int_id(cov_path.At(int_pos) ) );
    		if (unique_edge(cov_path.At(int_pos)) && g_.length(cov_path.At(int_pos)) > g_.k())
    			return true;
    		int_pos--;
    	}
    	DEBUG("path has no unique prefix");
    	return false;
    }

    bool unique_edge(EdgeId e){
    	//TODO: hack
    	DEBUG("running unique_edge :"<< g_.int_id(e));
    	if (g_.length(e) > cfg::get().rr.max_repeat_length) return true;
        auto supporting_paths = coverageMap_.GetCoveringPaths(e);
        for (auto it = supporting_paths.begin(); it != supporting_paths.end(); ++it) {
        	auto positions = (*it)->FindAll(e);
        	if (positions.size() > 1) return false;
        	for (auto jit = it; jit!= supporting_paths.end(); jit ++) {
            	auto j_positions = (*jit)->FindAll(e);
            	if (j_positions.size() > 1) return false;
            	double w1 = (*it)->getWeight();
            	double w2 = (*it)->getWeight();
            	if (w1 < 1.5 || w2 < 1.5 || w1/w2 > 4 || w2/w1 > 4)
            		 continue;
				if (!consistent_path(**it, **jit, positions[0], j_positions[0]))
					return false;
        	}
        }
        DEBUG("edge was unique!");
        return true;
    }

    bool covered_path(BidirectionalPath& path, BidirectionalPath& cov_path, size_t pos){
    	int cur_pos1 = path.Size() - 1;
    	int cur_pos2 = pos;
    	while (cur_pos1 >= 0 && cur_pos2 >=0){
    		if (path.At(cur_pos1) == cov_path.At(cur_pos2)){
    			cur_pos1--;
    			cur_pos2--;
    		} else {
    			return false;
    		}
    	}
    	return true;
    }

    bool consistent_path(BidirectionalPath& path, BidirectionalPath& cov_path, size_t pos, size_t cov_pos){
    	DEBUG("running consistent path");
    	path.Print();
    	cov_path.Print();
    	DEBUG(pos <<" " << path.Size());
    	DEBUG(cov_pos <<" " << cov_path.Size());
    	int cur_pos1 = pos;
    	int cur_pos2 = cov_pos;
    	while (cur_pos1 >= 0 && cur_pos2 >=0){
    		if (path.At(cur_pos1) == cov_path.At(cur_pos2)){
    			cur_pos1--;
    			cur_pos2--;
    		} else {
    			DEBUG("path are not consitent at "<< cur_pos1 <<" "<<cur_pos2);
    			return false;
    		}
    	}
    	cur_pos1 = pos;
		cur_pos2 = cov_pos;
		while (cur_pos1 < path.Size() && cur_pos2 < cov_path.Size()){
			if (path.At(cur_pos1) == cov_path.At(cur_pos2)){
				cur_pos1++;
				cur_pos2++;
			} else {
				DEBUG("path are not consitent at "<< cur_pos1 <<" "<<cur_pos2);
				return false;
			}
		}
		DEBUG("paths are consistent");
    	return true;
    }

    vector<pair<EdgeId, double> > to_vector(map<EdgeId, double>& candidates){
    	vector<pair<EdgeId, double> > result;
    	while (candidates.size() > 0){
    		double max = 0;
    		EdgeId max_edge = candidates.begin()->first;
    		for (auto iter = candidates.begin(); iter != candidates.end(); ++iter){
    			if (iter->second > max){
    				max = iter->second;
    				max_edge = iter->first;
    			}
    		}
    		result.push_back(make_pair(max_edge, max));
    		candidates.erase(max_edge);
    	}
    	return result;
    }
};


}


#endif /* EXTENSION_HPP_ */

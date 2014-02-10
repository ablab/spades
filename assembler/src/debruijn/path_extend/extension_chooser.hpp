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

#include <cfloat>
#include <iostream>
#include <fstream>
#include "weight_counter.hpp"
#include "pe_utils.hpp"
#include "next_path_searcher.hpp"

namespace path_extend {

typedef std::multimap<double, EdgeWithDistance> AlternativeConteiner;


class PathAnalyzer {

protected:
    const Graph& g_;

public:
    PathAnalyzer(const Graph& g): g_(g) {
    }

    int ExcludeTrivial(const BidirectionalPath& path, std::map<size_t, double>& edges, int from = -1) {
        int edgeIndex = (from == -1) ? (int) path.Size() - 1 : from;
        if ((int) path.Size() <= from) {
            return edgeIndex;
        }
        VertexId currentVertex = g_.EdgeEnd(path[edgeIndex]);
        while (edgeIndex >= 0 && g_.CheckUniqueIncomingEdge(currentVertex)) {
            EdgeId e = g_.GetUniqueIncomingEdge(currentVertex);
            currentVertex = g_.EdgeStart(e);

            edges.insert(make_pair((size_t)edgeIndex, 0.0));
            --edgeIndex;
        }
        return edgeIndex;
    }

    int ExcludeTrivialWithBulges(const BidirectionalPath& path, std::map<size_t, double>& edges) {
        edges.clear();

        if (path.Empty()) {
            return 0;
        }

        int lastEdge = (int) path.Size() - 1;
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

protected:
    DECL_LOGGER("ExtensionChooser")
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
    const Graph& g_;

    WeightCounter * wc_;

    PathAnalyzer analyzer_;

    double prior_coeff_;

    bool excludeTrivial_;
    bool excludeTrivialWithBulges_;
    DECL_LOGGER("ExtensionChooser");

    std::vector<ExtensionChooserListener *> listeners_;

public:
    ExtensionChooser(const Graph& g, WeightCounter * wc = 0, double priority = 0.0): g_(g), wc_(wc), analyzer_(g), prior_coeff_(priority),
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

    void ClearExcludedEdges() {
        wc_->GetExcludedEdges().clear();
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
    void RemoveTrivial(BidirectionalPath& path){
    	wc_->GetExcludedEdges().clear();
        if (excludeTrivialWithBulges_)
        {
            analyzer_.ExcludeTrivialWithBulges(path, wc_->GetExcludedEdges());
        }
        else if (excludeTrivial_)
        {
            analyzer_.ExcludeTrivial(path, wc_->GetExcludedEdges());
        }
    }


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

    virtual EdgeContainer Filter(BidirectionalPath& /*path*/, EdgeContainer& edges) {
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

	void RemoveTrivialAndCommon(BidirectionalPath& path, EdgeContainer& edges) {
        ClearExcludedEdges();
        if (edges.size() < 2) {
            return;
        }
        RemoveTrivial(path);
        int index = (int) path.Size() - 1;
        std::map<size_t, double>& excluded_edges = wc_->GetExcludedEdges();
        while (index >= 0) {
            if (excluded_edges.find(index) != excluded_edges.end()) {
                index--;
                continue;
            }
            EdgeId path_edge = path[index];
            double min_ideal_w = wc_->CountIdealInfo(path_edge, edges.at(0).e_,
                                                     path.LengthAt(index));
            bool common = true;
            for (size_t i = 0; i < edges.size(); ++i) {
                double ideal_weight = wc_->CountIdealInfo(path_edge,
                                                          edges.at(i).e_,
                                                          path.LengthAt(index));
                min_ideal_w = std::min(min_ideal_w, ideal_weight);
                if (!wc_->PairInfoExist(path_edge, edges.at(i).e_,
                                        (int) path.LengthAt(index))) {
                    common = false;
                }
            }
            if (common) {
                DEBUG("common info from " << index);
                excluded_edges.insert(make_pair((size_t) index, 0.0));
            } else {
                excluded_edges.insert(make_pair((size_t) index, min_ideal_w));
            }
            index--;
        }
        stringstream not_excl;
        not_excl << "not excluded edges ";
        for (size_t i = 0; i < path.Size(); ++i) {
            if (excluded_edges.find(i) != excluded_edges.end() && excluded_edges[i] > 0.0) {
                not_excl << i << " " << excluded_edges[i] << " , ";
            }
        }
        DEBUG(not_excl.str());
    }

	void FindWeights(BidirectionalPath& path, EdgeContainer& edges,
			AlternativeConteiner& weights) {
		for (auto iter = edges.begin(); iter != edges.end(); ++iter) {
			double weight = wc_->CountWeight(path, iter->e_);
			weights.insert(std::make_pair(weight, *iter));
			DEBUG("Candidate " << g_.int_id(iter->e_) << " weight " << weight << " length " << g_.length(iter->e_));
		}
		NotifyAll(weights);
	}

	void FindPossibleEdges(AlternativeConteiner& weights, EdgeContainer& top,
			double max_weight) {
		auto possibleEdge = weights.lower_bound(max_weight / prior_coeff_);
		for (auto iter = possibleEdge; iter != weights.end(); ++iter) {
			top.push_back(iter->second);
		}
	}

	EdgeContainer FindFilteredEdges(BidirectionalPath& path,
			EdgeContainer& edges) {
		AlternativeConteiner weights;
		FindWeights(path, edges, weights);
		EdgeContainer top;
		auto maxWeight = (--weights.end())->first;
		FindPossibleEdges(weights, top, maxWeight);
		EdgeContainer result;
		if (top.size() >= 1 && wc_->IsExtensionPossible(maxWeight)) {
			result = top;
		}
		return result;
	}
public:
	SimpleExtensionChooser(const Graph& g, WeightCounter * wc, double priority) :
			ExtensionChooser(g, wc, priority) {

	}

	virtual EdgeContainer Filter(BidirectionalPath& path,
			EdgeContainer& edges) {
	    DEBUG("Paired-end extension chooser");
		if (edges.empty()) {
			return edges;
		}
		RemoveTrivial(path);
		path.Print();
		EdgeContainer result = edges;
		bool first_time = true;
		bool changed = true;
		if (first_time || (result.size() > 1 && changed)) {
		    first_time = false;
			RemoveTrivialAndCommon(path, result);
			EdgeContainer new_result = FindFilteredEdges(path, result);
			if (new_result.size() == result.size()) {
				changed = false;
			}
			result = new_result;
		}
		if (result.size() == 1) {
            DEBUG("Paired-end extension chooser helped");
        }
		return result;
	}

};



class ScaffoldingExtensionChooser: public ExtensionChooser {

    bool cluster_info_;

    static bool compare(pair<int,double> a, pair<int,double> b)
    {
        if (a.first < b.first) return true;
        else return false;
    }

public:
	ScaffoldingExtensionChooser(const Graph& g, WeightCounter * wc, double priority, bool cluster_info = true): ExtensionChooser(g, wc, priority), cluster_info_(cluster_info) {
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
        return (int) round(dist);
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
        return (int) round(sqrt(dev / sum));
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

			for (size_t l = 0; l < weights.size(); ++l){
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
                histogram = FilterHistogram(histogram, mean - (5 - j)  * dev, mean + (5 - j) * dev, (int) round(cutoff));
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
                    DEBUG("scaffolding " << g_.int_id(edges[i].e_) << " gap " << gap);
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

        if (cluster_info_) {
            FindBestFittedEdgesForClustered(path, edges, result);
        } else {
            FindBestFittedEdges(path, edges, result);
        }

        return result;
    }
};

inline bool EdgeWithWeightCompareReverse(const pair<EdgeId, double>& p1,
                                      const pair<EdgeId, double>& p2) {
    return p1.second > p2.second;
}

class UniqueEdgeAnalyzer {
protected:
    DECL_LOGGER("ExtensionChooser")

public:
    UniqueEdgeAnalyzer(const Graph& g, const GraphCoverageMap& cov_map,
                       double filter_threshold, double prior_threshold)
            : g_(g),
              cov_map_(cov_map),
              filter_threshold_(filter_threshold),
              prior_threshold_(prior_threshold),
              unique_edges_founded_(false) { }

    bool IsUnique(EdgeId e) {
        if (!unique_edges_founded_) {
            FindAllUniqueEdges();
        }
        unique_edges_founded_ = true;
        return unique_edges_.count(e) > 0;
    }

private:
    bool UniqueEdge(EdgeId e) {
        if (g_.length(e) > cfg::get().max_repeat_length)
            return true;
        DEBUG("Analyze unique edge " << g_.int_id(e));
        if (cov_map_.size() == 0) {
            return false;
        }
        auto cov_paths = cov_map_.GetCoveringPaths(e);
        for (auto it1 = cov_paths.begin(); it1 != cov_paths.end(); ++it1) {
            auto pos1 = (*it1)->FindAll(e);
            if (pos1.size() > 1) {
                DEBUG("***not unique " << g_.int_id(e) << " len " << g_.length(e) << "***");
                return false;
            }
            for (auto it2 = it1; it2 != cov_paths.end(); it2++) {
                auto pos2 = (*it2)->FindAll(e);
                if (pos2.size() > 1) {
                    DEBUG("***not unique " << g_.int_id(e) << " len " << g_.length(e) << "***");
                    return false;
                }
                if (!ConsistentPath(**it1, pos1[0], **it2, pos2[0])) {
                    if (CheckInconsistence(**it1, pos1[0], **it2, pos2[0],
                                           cov_paths)) {
                        DEBUG("***not unique " << g_.int_id(e) << " len " << g_.length(e) << "***");
                        return false;
                    }
                }
            }
        }
        DEBUG("***edge " << g_.int_id(e) << " is unique.***");
        return true;
    }

    bool ConsistentPath(const BidirectionalPath& path1, size_t pos1,
                        const BidirectionalPath& path2, size_t pos2) const {
        return EqualBegins(path1, pos1, path2, pos2, false)
                && EqualEnds(path1, pos1, path2, pos2, false);
    }
    bool SignificantlyDiffWeights(double w1, double w2) const {
        if (w1 > filter_threshold_ and w2 > filter_threshold_) {
            if (w1 > w2 * prior_threshold_ or w2 > w1 * prior_threshold_) {
                return true;
            }
            return false;
        }
        return true;
    }

    bool CheckInconsistence(
            const BidirectionalPath& path1, size_t pos1,
            const BidirectionalPath& path2, size_t pos2,
            const std::set<BidirectionalPath*>& cov_paths) const {
        size_t first_diff_pos1 = FirstNotEqualPosition(path1, pos1, path2, pos2, false);
        size_t first_diff_pos2 = FirstNotEqualPosition(path2, pos2, path1, pos1, false);
        if (first_diff_pos1 != -1UL && first_diff_pos2 != -1UL) {
            const BidirectionalPath cand1 = path1.SubPath(first_diff_pos1,
                                                          pos1 + 1);
            const BidirectionalPath cand2 = path2.SubPath(first_diff_pos2,
                                                          pos2 + 1);
            std::pair<double, double> weights = GetSubPathsWeights(cand1, cand2,
                                                                   cov_paths);
            DEBUG("Not equal begin " << g_.int_id(path1.At(first_diff_pos1)) << " weight " << weights.first << "; " << g_.int_id(path2.At(first_diff_pos2)) << " weight " << weights.second);
            if (!SignificantlyDiffWeights(weights.first, weights.second)) {
                DEBUG("not significantly different");
                return true;
            }
        }
        size_t last_diff_pos1 = LastNotEqualPosition(path1, pos1, path2, pos2, false);
        size_t last_diff_pos2 = LastNotEqualPosition(path2, pos2, path1, pos1, false);
        if (last_diff_pos1 != -1UL) {
            const BidirectionalPath cand1 = path1.SubPath(pos1,
                                                          last_diff_pos1 + 1);
            const BidirectionalPath cand2 = path2.SubPath(pos2,
                                                          last_diff_pos2 + 1);
            std::pair<double, double> weights = GetSubPathsWeights(cand1, cand2,
                                                                   cov_paths);
            DEBUG("Not equal end " << g_.int_id(path1.At(last_diff_pos1)) << " weight " << weights.first << "; " << g_.int_id(path2.At(last_diff_pos2)) << " weight " << weights.second);
            if (!SignificantlyDiffWeights(weights.first, weights.second)) {
                DEBUG("not significantly different");
                return true;
            }
        }
        return false;
    }
    std::pair<double, double> GetSubPathsWeights(
            const BidirectionalPath& cand1, const BidirectionalPath& cand2,
            const std::set<BidirectionalPath*>& cov_paths) const {
        double weight1 = 0.0;
        double weight2 = 0.0;
        for (auto iter = cov_paths.begin(); iter != cov_paths.end(); ++iter) {
            BidirectionalPath* path = *iter;
            if (ContainSubPath(*path, cand1)) {
                weight1 += path->GetWeight();
            } else if (ContainSubPath(*path, cand2)) {
                weight2 += path->GetWeight();
            }
        }
        return std::make_pair(weight1, weight2);
    }

    bool ContainSubPath(const BidirectionalPath& path,
                        const BidirectionalPath& subpath) const {
        for (size_t i = 0; i < path.Size(); ++i) {
            if (path.CompareFrom(i, subpath))
                return true;
        }
        return false;
    }

    void FindAllUniqueCoverageEdges(){
       if (cfg::get().ds.single_cell) {
           return;
       }
       double sum_cov = 0;
       size_t sum_len = 0;
       size_t total_len = 0;
       for (auto iter = g_.ConstEdgeBegin(); !iter.IsEnd(); ++iter) {
           total_len += g_.length(*iter);
           if (g_.length(*iter) >= cfg::get().max_repeat_length) {
               sum_cov += g_.coverage(*iter) * (double)g_.length(*iter);
               sum_len += g_.length(*iter);
           }
       }
       if (sum_len * 4 < total_len) return;
       sum_cov /= (double)sum_len;
       DEBUG("average coverage of long edges: " << sum_cov) ;
       for (auto iter = g_.ConstEdgeBegin(); !iter.IsEnd(); ++iter) {
           if (g_.length(*iter) > 500 && (double)g_.coverage(*iter) < 1.2 * sum_cov) {
               if (unique_edges_.find(*iter) == unique_edges_.end()) {
                   unique_edges_.insert(*iter);
                   unique_edges_.insert(g_.conjugate(*iter));
                   DEBUG("Added coverage based unique edge " << g_.int_id(*iter) << " len "<< g_.length(*iter) << " " << g_.coverage(*iter));
               }
           }
       }
   }


    void FindAllUniqueEdges() {
       DEBUG("Looking for unique edges");
       for (auto iter = g_.ConstEdgeBegin(); !iter.IsEnd(); ++iter) {
           if (UniqueEdge(*iter)) {
               unique_edges_.insert(*iter);
               unique_edges_.insert(g_.conjugate(*iter));
           }
       }
       unique_edges_founded_ = true;
       DEBUG("coverage based uniqueness started");
       FindAllUniqueCoverageEdges();
       DEBUG("Unique edges are found");
    }

    const Graph& g_;
    const GraphCoverageMap& cov_map_;
    double filter_threshold_;
    double prior_threshold_;bool unique_edges_founded_;
    std::set<EdgeId> unique_edges_;
};

class SimpleScaffolding {
public:
    SimpleScaffolding(const Graph& g) : g_(g) {}

    BidirectionalPath FindMaxCommonPath(const vector<BidirectionalPath*>& paths,
                                        size_t max_diff_len) {
        BidirectionalPath max_end(g_);
        for (auto it1 = paths.begin(); it1 != paths.end(); ++it1) {
            BidirectionalPath* p1 = *it1;
            //DEBUG("begin find subpath for");
            //p1->PrintInString();
            for (size_t i = 0; i < p1->Size(); ++i) {
                if (p1->Length() - p1->LengthAt(i) > max_diff_len) {
                    //DEBUG("more max diff " << i << " len " << p1->Length() - p1->LengthAt(i))
                    break;
                }
                bool contain_all = true;
                for (size_t i1 = i + 1; i1 <= p1->Size() && contain_all; ++i1) {
                    BidirectionalPath subpath = p1->SubPath(i, i1);
                    //DEBUG("analyze subpath ");
                    //subpath.PrintInString();
                    for (auto it2 = paths.begin();  it2 != paths.end() && contain_all; ++it2) {
                        BidirectionalPath* p2 = *it2;
                        vector<size_t> positions2 = p2->FindAll(subpath.At(0));
                        bool contain = false;
                        for (size_t ipos2 = 0; ipos2 < positions2.size(); ++ipos2) {
                            size_t pos2 = positions2[ipos2];
                            if (p2->Length() - p2->LengthAt(pos2) <= max_diff_len
                                    && EqualEnds(subpath, 0, *p2, pos2, false)) {
                                contain = true;
                                break;
                            } else {
                                //DEBUG("not consist with path " << pos2);
                                //p2->PrintInString();
                            }
                        }
                        if (!contain) {
                            //DEBUG("not contains all " << i);
                            contain_all = false;
                        }
                    }
                    if (contain_all && (i1 - i) >= max_end.Size()) {
                        //DEBUG("common end " << i << " " << i1)
                        max_end.Clear();
                        max_end.PushBack(subpath);
                    }
                }
            }
        }
        return max_end;
    }

private:
    const Graph& g_;

};
class LongReadsExtensionChooser : public ExtensionChooser {
public:
    LongReadsExtensionChooser(const Graph& g, PathContainer& pc,
                              double filtering_threshold,
                              double weight_priority_threshold,
                              double unique_edge_priority_threshold)
            : ExtensionChooser(g, 0, .0),
              filtering_threshold_(filtering_threshold),
              weight_priority_threshold_(weight_priority_threshold),
              cov_map_(g, pc),
              unique_edge_analyzer_(g, cov_map_, filtering_threshold, unique_edge_priority_threshold),
              simple_scaffolding_(g) {

    }

    /* Choose extension as correct only if we have reads that traverse a unique edge from the path and this extension.
     * Edge is unique if all reads mapped to this edge are consistent.
     * Two reads are consistent if they can form one path in the graph.
     */
    virtual EdgeContainer Filter(BidirectionalPath& path,
                                 EdgeContainer& edges) {
        if (edges.empty()) {
            return edges;
        }DEBUG("We in Filter of LongReadsExtensionChooser");
        path.Print();
        map<EdgeId, double> weights_cands;
        for (auto it = edges.begin(); it != edges.end(); ++it) {
            weights_cands.insert(make_pair(it->e_, 0.0));
        }
        set<EdgeId> filtered_cands;
        map<EdgeId, set<BidirectionalPath*> > support_paths_ends;
        auto support_paths = cov_map_.GetCoveringPaths(path.Back());
        for (auto it = support_paths.begin(); it != support_paths.end(); ++it) {
            auto positions = (*it)->FindAll(path.Back());
            for (size_t i = 0; i < positions.size(); ++i) {
                if ((int) positions[i] < (int) (*it)->Size() - 1
                        && EqualBegins(path, (int) path.Size() - 1, **it,
                                       positions[i], true)) {

                    if (UniqueBackPath(**it, positions[i])) {
                        EdgeId next = (*it)->At(positions[i] + 1);
                        weights_cands[next] += (*it)->GetWeight();
                        filtered_cands.insert(next);
                        if (support_paths_ends.count(next) == 0){
                            support_paths_ends[next] = set<BidirectionalPath*>();
                        }
                        support_paths_ends[next].insert(new BidirectionalPath((*it)->SubPath(positions[i] + 1)));
                    }
                }
            }
        }
        DEBUG("Candidates");
        for (auto iter = weights_cands.begin(); iter != weights_cands.end();
                ++iter) {
            DEBUG("Candidate " << g_.int_id(iter->first) << " weight " << iter->second);
        }
        vector<pair<EdgeId, double> > sort_res = MapToSortVector(weights_cands);
        DEBUG("sort res " << sort_res.size() << " tr " << weight_priority_threshold_);
        if (sort_res.size() < 1 || sort_res[0].second < filtering_threshold_) {
            filtered_cands.clear();
        } else if (sort_res.size() > 1
                && sort_res[0].second
                        > weight_priority_threshold_ * sort_res[1].second) {
            filtered_cands.clear();
            filtered_cands.insert(sort_res[0].first);
        } else if (sort_res.size() > 1) {
            for (size_t i = 0; i < sort_res.size(); ++i) {
                if (sort_res[i].second * weight_priority_threshold_ < sort_res[0].second){
                    filtered_cands.erase(sort_res[i].first);
                }
            }
        }
        DEBUG("fil size " << filtered_cands.size());
        if (filtered_cands.size() > 1) {
            filtered_cands = SimpleScaffold(filtered_cands, support_paths_ends, sort_res);
        }
        EdgeContainer result;
        for (auto it = edges.begin(); it != edges.end(); ++it) {
            if (filtered_cands.find(it->e_) != filtered_cands.end()) {
                result.push_back(*it);
            }
        }
        if (result.size() != 1) {
            DEBUG("Long reads doesn't help =(");
        }
        return result;
    }

private:
    set<EdgeId> SimpleScaffold(set<EdgeId>& candidates, map<EdgeId, set<BidirectionalPath*> >& ends, vector<pair<EdgeId, double> > sort_weights) {
        vector<BidirectionalPath*> all_ends;
        for (auto it = ends.begin(); it != ends.end(); ++it) {
            if (candidates.find(it->first) == candidates.end()) {
                continue;
            }
            all_ends.insert(all_ends.end(), it->second.begin(), it->second.end());
        }
        DEBUG("Scaffolding");
        BidirectionalPath end = simple_scaffolding_.FindMaxCommonPath(all_ends, g_.k() + 40); //TODO: move to config
        end.Print();
        if (end.Size() > 0){
            set<EdgeId> res;
            res.insert(sort_weights[0].first);
            DEBUG("scaffolding helps");
            return res;
        }
        return candidates;
    }
    bool UniqueBackPath(const BidirectionalPath& path, size_t pos) {
        int int_pos = (int) pos;
        while (int_pos >= 0) {
            if (unique_edge_analyzer_.IsUnique(path.At(int_pos)) > 0)
                return true;
            int_pos--;
        }
        return false;
    }

    vector<pair<EdgeId, double> > MapToSortVector(
            map<EdgeId, double>& map) const {
        vector<pair<EdgeId, double> > result1(map.begin(), map.end());
        std::sort(result1.begin(), result1.end(), EdgeWithWeightCompareReverse);
        return result1;
    }

    double filtering_threshold_;
    double weight_priority_threshold_;
    const GraphCoverageMap cov_map_;
    UniqueEdgeAnalyzer unique_edge_analyzer_;
    SimpleScaffolding simple_scaffolding_;
};

class MatePairExtensionChooser : public ExtensionChooser {
public:
    MatePairExtensionChooser(const Graph& g, PairedInfoLibrary& lib,
                              const PathContainer& paths)
            : ExtensionChooser(g, 0, .0),
              g_(g),
              lib_(lib),
              search_dist_(lib.GetISMax()),
              weight_counter_(g, lib, 10),
              cov_map_(g_, paths),
              path_searcher_(g_, cov_map_, lib_.GetISMax(), PathsWeightCounter(g, lib, 30)),
              unique_edge_analyzer_(g, cov_map_, 0., 1000.),
              simple_scaffolder_(g) {
    }
    virtual EdgeContainer Filter(BidirectionalPath& path,
                                 EdgeContainer& init_edges) {
        DEBUG("mp chooser");
        path.Print();
        if (path.Length() < lib_.GetISMin()) {
            return EdgeContainer();
        }
        EdgeContainer edges = init_edges;
        if (IsBulge(init_edges)){
        	 edges = TryResolveBulge(path, init_edges);
        	 DEBUG("bulge " << edges.size())
        }
        map<EdgeId, BidirectionalPath*> best_paths;
        for (size_t iedge = 0; iedge < edges.size(); ++iedge) {
            set<BidirectionalPath*> following_paths = path_searcher_.FindNextPaths(path, edges[iedge].e_);
            vector<BidirectionalPath*> max_weighted = MaxWeightedPath(path, following_paths);
            if (max_weighted.size() == 0) {
                DEBUG("too much paths or tip");
                DeleteMapWithPaths(best_paths);
                DeleteNextPaths(following_paths);
                best_paths.clear();
                break;
            } else {
                best_paths[edges[iedge].e_] = new BidirectionalPath(*max_weighted[0]);
            }
            DeleteNextPaths(following_paths);
        }

        set<BidirectionalPath*> next_paths;
        if (edges.size() == 0) {
        	DEBUG("scaffolding edges size " << edges.size())
            next_paths = path_searcher_.FindNextPaths(path, path.Back());
        } else if (best_paths.size() == edges.size()) {
            for (size_t iedge = 0; iedge < edges.size(); ++iedge) {
                if (best_paths.count(edges[iedge].e_) > 0){
                    next_paths.insert(best_paths[edges[iedge].e_]);
                }
            }
        }
        EdgeContainer result = ChooseBest(path, next_paths);
        if (result.size() != 1) {
            DEBUG("scaffold tree");
            result = ScaffoldTree(path);
        }
        DeleteNextPaths(next_paths);
        if (result.size() != 1) {
            DEBUG("nobody can extend " << g_.int_id(path.Back()));
        }
        return result;
    }
private:
    EdgeContainer ScaffoldTree(BidirectionalPath& path) {
        DEBUG("try scaffold tree");
        set<BidirectionalPath*> next_paths = path_searcher_.ScaffoldTree(path);
        VERIFY(next_paths.size() <= 1);
        if (next_paths.size() != 1) {
            DeleteNextPaths(next_paths);
            return EdgeContainer();
        }
        DEBUG("Path to add dirty ");
        BidirectionalPath* res = NULL;
        for (BidirectionalPath* p : next_paths) {
            res = p;
            p->Print();
            for (size_t i = 0; i < p->Size() - 1; ++i) {
                path.PushBack(p->At(i), p->GapAt(i));
            }
        }
        EdgeContainer result = EdgeContainer();
        if (res != NULL) {
            DEBUG("Finally the path, and an edge to add " << g_.int_id(res->Back()) << " : " << g_.length(res->Back()));
            path.Print();
            result = EdgeContainer(1, EdgeWithDistance(res->Back(), res->GapAt(res->Size() - 1)));
        }
        DeleteNextPaths(next_paths);
        return result;
    }

	bool IsBulge(EdgeContainer& edges) {
		if (edges.size() == 0)
			return false;
		for (EdgeWithDistance e : edges) {
			if (!InBuble(e.e_, g_))
				return false;
		}
		return true;
	}

    map<EdgeId, double> FindBulgeWeights(BidirectionalPath& p, EdgeContainer& edges) const {
        map<EdgeId, double> result;
        for (size_t i = 0; i < edges.size(); ++i) {
            result[edges[i].e_] = 0.0;
        }
        for (size_t i = 0; i < p.Size(); ++i) {
            bool common = true;
            bool common_ideal = true;
            for (EdgeWithDistance e : edges) {
                common_ideal = common_ideal && weight_counter_.HasIdealPI(p.At(i), e.e_, (int) p.LengthAt(i));
                common = common && weight_counter_.HasPI(p.At(i), e.e_, (int) p.LengthAt(i));
            }
            if (!common_ideal || common) {
                continue;
            }
            for (size_t j = 0; j < edges.size(); ++j) {
                result[edges[j].e_] += weight_counter_.PI(p.At(i), edges[j].e_, (int) p.LengthAt(i));
            }
        }
        return result;
    }

    EdgeContainer TryResolveBulge(BidirectionalPath& p, EdgeContainer& edges) const {
        map<EdgeId, double> weights = FindBulgeWeights(p, edges);
        double max_w = 0.0;
        EdgeWithDistance max = *edges.begin();
        for (EdgeWithDistance e : edges) {
            double w = weights[e.e_];
            DEBUG("bulge " << g_.int_id(e.e_) << " w = " << w);
            if (w > max_w) {
                max_w = w;
                max = e;
            }
        }
        for (EdgeWithDistance e : edges) {
            if (e.e_ == max.e_) {
                continue;
            }
            if (math::eq(weights[e.e_], max_w)) {
                return edges;
            }
        }
        EdgeContainer result;
        result.push_back(max);
        return result;
    }

    void DeleteNextPaths(set<BidirectionalPath*>& paths) {
        for (auto i = paths.begin(); i != paths.end(); ++i) {
            delete (*i);
        }
    }
    void DeleteMapWithPaths(map<EdgeId, BidirectionalPath*> m) {
        for (auto i = m.begin(); i != m.end(); ++i){
            delete i->second;
        }
    }
    EdgeContainer ChooseBest(const BidirectionalPath& path,
                             set<BidirectionalPath*>& next_paths) {
        DEBUG("Try to choose from best paths...");
        vector<BidirectionalPath*> best_path = MaxWeightedPath(path,
                                                               next_paths);
        EdgeContainer result;
        if (best_path.size() == 1) {
            DEBUG("result is good");
            result.push_back(EdgeWithDistance((*best_path.begin())->At(0),
                                     (*best_path.begin())->GapAt(0)));
        } else if (best_path.size() > 1) {
            result = TryToScaffold(path, best_path);
        }
        return result;
    }
    bool HasUniqueEdges(const BidirectionalPath& init_path,
                        const BidirectionalPath& path1,
                        const set<size_t>& unique_init_edges) {
        for (size_t i1 = 0; i1 < init_path.Size(); ++i1) {
            if (unique_init_edges.find(i1) == unique_init_edges.end()){
                continue;
            }
            for (size_t i2 = 0; i2 < path1.Size(); ++i2) {
                int gap = (int)init_path.LengthAt(i1) + (int)path1.Length() - (int)path1.LengthAt(i2);
                if (unique_edge_analyzer_.IsUnique(init_path.At(i1))
                        && unique_edge_analyzer_.IsUnique(path1.At(i2))
                        && weight_counter_.HasPI(init_path.At(i1), path1.At(i2), gap)) {
                    DEBUG("has unique edge " << g_.int_id(init_path.At(i1)) << " " <<g_.int_id(path1.At(i2)));
                    return true;
                }
            }
        }
        return false;
    }

    bool SignificallyDifferentEdges(const BidirectionalPath& init_path,
                                    const BidirectionalPath& path1,
                                    const map<size_t, double>& pi1,
                                    double /*w1*/,
                                    const BidirectionalPath& path2,
                                    const map<size_t, double>& pi2,
                                    double /*w2*/,
                                    const set<size_t>& unique_init_edges) {
        double not_common_w1 = 0.0;
        double common_w = 0.0;
        for (auto iter = pi1.begin(); iter != pi1.end(); ++iter) {
            auto iter2 = pi2.find(iter->first);
            double w = 0.0;
            if (iter2 != pi2.end() && !math::eq(iter2->second, 0.0)) {
                w = min(iter2->second, iter->second);
            }
            not_common_w1 += iter->second - w;
            common_w += w;
        }

        if (common_w < 0.8 * (not_common_w1 + common_w)||
                (HasUniqueEdges(init_path, path1, unique_init_edges) && !HasUniqueEdges(init_path, path2, unique_init_edges))) {
            DEBUG("common_w " << common_w  << " sum*0.8  = " << 0.8 * (not_common_w1 + common_w))
            return true;
        }
        return false;
    }
    set<size_t> FindNotCommonEdges(const BidirectionalPath& path, const std::map<BidirectionalPath*, map<size_t, double> >& all_pi) {
        set<size_t> res;
        for (size_t i = 0; i < path.Size(); ++i) {
            bool not_common = true;
            bool exist = false;
            for (auto iter = all_pi.begin(); iter != all_pi.end(); ++iter){
                const map<size_t, double>& info = iter->second;
                if (info.count(i) > 0 && math::gr(info.at(i), 0.0)) {
                    if (exist) {
                        not_common = false;
                        break;
                    } else {
                        exist = true;
                    }
                }
            }
            if (not_common && exist) {
                res.insert(i);
            }
        }
        return res;
    }
    void DeleteSmallWeights(const BidirectionalPath& path,
                            set<BidirectionalPath*>& paths,
                            std::map<BidirectionalPath*, map<size_t, double> >& all_pi) {
        double max_weight = 0.0;
        BidirectionalPath* max_path = NULL;
        for (auto iter = paths.begin(); iter != paths.end(); ++iter) {
            if ((*iter)->GetWeight() > max_weight) {
                max_weight = max(max_weight, (*iter)->GetWeight());
                max_path = *iter;
            }
        }
        set<BidirectionalPath*> to_del;
        for (auto iter = paths.begin(); iter != paths.end(); ++iter) {
            if (math::gr(max_weight, (*iter)->GetWeight() * 1.5) &&
                    SignificallyDifferentEdges(path, *max_path, all_pi.find(max_path)->second, max_weight,
                                               *    *iter, all_pi.find(*iter)->second, (*iter)->GetWeight(),
                                               FindNotCommonEdges(path, all_pi)))
                to_del.insert(*iter);
        }
        for (BidirectionalPath* p : to_del){
            paths.erase(p);
            all_pi.erase(p);
        }
    }

    void DeleteCommonPi(const BidirectionalPath& path,
            const std::map<BidirectionalPath*, map<size_t, double> >& all_pi) {
        weight_counter_.ClearCommonWeight();
        for (size_t i = 0; i < path.Size(); ++i) {
            double common = DBL_MAX;
            for (auto iter = all_pi.begin(); iter != all_pi.end(); ++iter) {
                if (iter->second.count(i) == 0) {
                    common = 0.0;
                    break;
                } else {
                    common = min(common, iter->second.at(i));
                }
            }
            weight_counter_.SetCommonWeightFrom(i, common);
        }
    }

    void CountAllPairInfo(
            const BidirectionalPath& path,
            const set<BidirectionalPath*>& next_paths,
            std::map<BidirectionalPath*, map<size_t, double> >& result) {
        result.clear();
        if (next_paths.size() == 0){
            return;
        }
        size_t common_begin = 0;
        BidirectionalPath* p = *next_paths.begin();
        bool contain_all = true;
        while (contain_all && common_begin < p->Size()) {
            EdgeId e = p->At(common_begin);
            for (BidirectionalPath* next : next_paths) {
                if (common_begin >= next->Size() || next->At(common_begin) != e) {
                    contain_all = false;
                    break;
                }
            }
            if (!contain_all) {
                break;
            }
            common_begin++;
        }
        DEBUG("common begin " << common_begin);
        for (BidirectionalPath* next : next_paths) {
            result[next] = weight_counter_.FindPairInfoFromPath(path, 0, path.Size(), *next, common_begin, next->Size());
        }
    }

    void CountWeightsAndFilter(
            const BidirectionalPath& path,
            set<BidirectionalPath*>& next_paths, bool delete_small_w) {
        std::map<BidirectionalPath*, map<size_t, double> > all_pi;
        CountAllPairInfo(path, next_paths, all_pi);
        DeleteCommonPi(path, all_pi);
        for (BidirectionalPath* next : next_paths) {
            next->SetWeight((float)weight_counter_.CountPairInfo(path, 0, path.Size(),
                                                         *next, 0,
                                                         next->Size()));
        }
        if (delete_small_w) {
            DeleteSmallWeights(path, next_paths, all_pi);
        }
    }

    struct PathWithWeightSort {
        PathWithWeightSort(MatePairExtensionChooser& mp_chooser, const BidirectionalPath& path, std::map<BidirectionalPath*, map<size_t, double> >& all_pi)
                : mp_chooser_(mp_chooser),
                  path_(path),
                  all_pi_(all_pi) {
            not_common_ = mp_chooser_.FindNotCommonEdges(path_, all_pi_);
        }

        bool operator()(const BidirectionalPath* p1, const BidirectionalPath* p2) {
            if (mp_chooser_.HasUniqueEdges(path_, *p1, not_common_) && !mp_chooser_.HasUniqueEdges(path_, *p2, not_common_)) {
                return true;
            }
            if (mp_chooser_.HasUniqueEdges(path_, *p2, not_common_) && !mp_chooser_.HasUniqueEdges(path_, *p1, not_common_)) {
                return false;
            }
            if (!math::eq(p1->GetWeight(), p2->GetWeight())) {
                return math::gr(p1->GetWeight(), p2->GetWeight());
            }
            if (!math::eq(p1->GetWeight(), p2->GetWeight())) {
                return math::gr(p1->GetWeight(), p2->GetWeight());
            }
            if (p1->Length() != p2->Length()) {
                return p1->Length() > p2->Length();
            }
            return p1->Size() > p2->Size();
        }
        MatePairExtensionChooser& mp_chooser_;
        const BidirectionalPath& path_;
        std::map<BidirectionalPath*, map<size_t, double> >& all_pi_;
        set<size_t> not_common_;
    };

    vector<BidirectionalPath*> SortResult(const BidirectionalPath& path,
            set<BidirectionalPath*>& next_paths) {
        DEBUG("sort result");
        std::map<BidirectionalPath*, map<size_t, double> > all_pi;
        CountAllPairInfo(path, next_paths, all_pi);
        CountWeightsAndFilter(path, next_paths, false);
        vector<BidirectionalPath*> to_sort(next_paths.begin(), next_paths.end());
        DEBUG("sort weight and filter");
        for (size_t i = 0; i < to_sort.size(); ++i){
               INFO("i " << i << " from " << to_sort.size());
               to_sort[i]->PrintInfo();
        }
        PathWithWeightSort comparator(*this, path, all_pi);
        std::sort(to_sort.begin(), to_sort.end(), comparator);
        DEBUG("end sort result");
        return to_sort;
    }

    vector<BidirectionalPath*> MaxWeightedPath(
            const BidirectionalPath& path,
            const set<BidirectionalPath*>& following_paths) {
        set<BidirectionalPath*> result(following_paths);
        set<BidirectionalPath*> prev_result;
        while (prev_result.size() != result.size()) {
            prev_result = result;
            DEBUG("iteration with paths " << result.size());
            CountWeightsAndFilter(path, result, true);
            if (result.size() == 0)
                result = prev_result;
            if (result.size() == 1)
                break;
        }
        if (result.size() == 0) {
            DEBUG("bad case");
            return vector<BidirectionalPath*>();
        }
        return SortResult(path, result);
    }

    BidirectionalPath ChooseFromEnds(const BidirectionalPath& path,
                                     const vector<BidirectionalPath*>& paths,
                                     const BidirectionalPath& end) {
        DEBUG("choose from ends " << paths.size());
        end.Print();
        vector<BidirectionalPath*> new_paths;
        vector<BidirectionalPath*> paths_to_cover;
        for (BidirectionalPath* p : paths) {
            int from = 0;
            int pos = p->FindFirst(end, from);
            while (pos > -1) {
                BidirectionalPath* new_p = new BidirectionalPath(path);
                BidirectionalPath* new_end = new BidirectionalPath(p->SubPath(0, pos + end.Size()));
                new_p->PushBack(*new_end);
                new_paths.push_back(new_p);
                paths_to_cover.push_back(new_end);
                from = pos + 1;
                pos = p->FindFirst(end, from);
            }
        }
        BidirectionalPath max = **new_paths.begin();
        size_t covered_edges_max = 0;
        size_t min_size = max.Size();
        for (BidirectionalPath* p: new_paths) {
            size_t cov_edges = 0;
            for (BidirectionalPath* e : paths_to_cover) {
                vector<size_t> poses = p->FindAll(e->Back());
                for (size_t pos : poses) {
                    if (EqualBegins(*p, pos, *e, e->Size() - 1, true)) {
                        cov_edges++;
                        break;
                    }
                }
            }
            if (cov_edges > covered_edges_max || (cov_edges == covered_edges_max && min_size > p->Size())) {
                DEBUG("cov_e " << cov_edges << " s " << p->Size());
                max.Clear();
                max.PushBack(*p);
                covered_edges_max = cov_edges;
                min_size = max.Size();
            }
        }
        for (BidirectionalPath* p: new_paths) {
            delete p;
        }
        for (BidirectionalPath* p : paths_to_cover) {
            delete p;
        }

        DEBUG("res");
        max.SubPath(path.Size()).Print();
        return max.SubPath(path.Size());
    }

    EdgeContainer TryToScaffold(const BidirectionalPath& path,
                                const vector<BidirectionalPath*>& paths) {
        DEBUG("Simple Scaffolding")
        if (paths.size() == 0) {
            return EdgeContainer();
        }
        for (BidirectionalPath* p : paths) {
            p->Print();
        }

        BidirectionalPath max_end = simple_scaffolder_.FindMaxCommonPath(paths, search_dist_);
        if (max_end.Size() == 0) {
            return EdgeContainer();
        }
        BidirectionalPath result_end = ChooseFromEnds(path, paths, max_end);
        int begin = result_end.FindFirst(max_end);
        weight_counter_.ClearCommonWeight();
        double common = weight_counter_.CountPairInfo(
                path, 0, path.Size(), result_end,
                begin, begin + max_end.Size(), false);
        double not_common = weight_counter_.CountPairInfo(
                path, 0, path.Size(), result_end, 0,
                begin, false);
        DEBUG("common " << common << " not common " << not_common << " max common " << begin << " " << begin + max_end.Size());
        result_end.Print();
        EdgeContainer result;
        size_t to_add = begin;
        size_t gap_length = result_end.Length() - result_end.LengthAt(to_add);
        DEBUG(" edge to add " << g_.int_id(result_end.At(to_add)) << " with length " << gap_length);
        result.push_back(EdgeWithDistance(result_end.At(to_add), gap_length));
        return result;
    }

    const Graph& g_;
    PairedInfoLibrary& lib_;
    size_t search_dist_;
    PathsWeightCounter weight_counter_;
    const GraphCoverageMap cov_map_;
    NextPathSearcher path_searcher_;
    UniqueEdgeAnalyzer unique_edge_analyzer_;
    SimpleScaffolding simple_scaffolder_;
};
}
#endif /* EXTENSION_HPP_ */

//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

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

typedef std::multimap<double, EdgeWithDistance> AlternativeContainer;

class PathAnalyzer {
    const Graph& g_;

public:
    PathAnalyzer(const Graph& g): g_(g) {
    }

    int ExcludeTrivial(const BidirectionalPath& path, std::set<size_t>& edges, int from = -1) const {
        int edgeIndex = (from == -1) ? (int) path.Size() - 1 : from;
        if ((int) path.Size() <= from) {
            return edgeIndex;
        }
        VertexId currentVertex = g_.EdgeEnd(path[edgeIndex]);
        while (edgeIndex >= 0 && g_.CheckUniqueIncomingEdge(currentVertex)) {
            EdgeId e = g_.GetUniqueIncomingEdge(currentVertex);
            currentVertex = g_.EdgeStart(e);

            edges.insert((size_t) edgeIndex);
            --edgeIndex;
        }
        return edgeIndex;
    }

    int ExcludeTrivialWithBulges(const BidirectionalPath& path, std::set<size_t>& edges) const {

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
    DECL_LOGGER("PathAnalyzer")
};


class ExtensionChooserListener {

public:

    virtual void ExtensionChosen(double weight) = 0;

    virtual void ExtensionChosen(const AlternativeContainer& alts) = 0;

    virtual ~ExtensionChooserListener() {

    }
};


class ExtensionChooser {

public:
    typedef std::vector<EdgeWithDistance> EdgeContainer;

protected:
    const Graph& g_;
    shared_ptr<WeightCounter> wc_;
    std::vector<ExtensionChooserListener *> listeners_;

private:
    double weight_threshold_;
    PathAnalyzer analyzer_;

    bool excludeTrivial_;
    bool excludeTrivialWithBulges_;


public:
    ExtensionChooser(const Graph& g, shared_ptr<WeightCounter> wc = nullptr, double weight_threshold = -1.): 
        g_(g), wc_(wc), 
        weight_threshold_(weight_threshold), analyzer_(g), 
        excludeTrivial_(true), excludeTrivialWithBulges_(true) {
    }

    virtual ~ExtensionChooser() {

    }

    virtual EdgeContainer Filter(const BidirectionalPath& path, const EdgeContainer& edges) const = 0;

    bool isExcludeTrivial() const
    {
        return excludeTrivial_;
    }

    bool isExcludeTrivialWithBulges() const
    {
        return excludeTrivialWithBulges_;
    }

    void setExcludeTrivial(bool excludeTrivial) {
        this->excludeTrivial_ = excludeTrivial;
    }

    void setExcludeTrivialWithBulges(bool excludeTrivialWithBulges) {
        this->excludeTrivialWithBulges_ = excludeTrivialWithBulges;
    }

    bool CheckThreshold(double weight) const {
        return math::ge(weight, weight_threshold_);
    }

    const PairedInfoLibraries& getLibs() {
        VERIFY(wc_);
        return wc_->getLibs();
    }

    void Subscribe(ExtensionChooserListener * listener) {
        listeners_.push_back(listener);
    }

    void NotifyAll(double weight) const {
        for (auto iter = listeners_.begin(); iter != listeners_.end(); ++iter) {
            (*iter)->ExtensionChosen(weight);
        }
    }

    void NotifyAll(const AlternativeContainer& alts) const {
        for (auto iter = listeners_.begin(); iter != listeners_.end(); ++iter) {
            (*iter)->ExtensionChosen(alts);
        }
    }

    bool WeighConterBased() const {
        return wc_ != nullptr;
    }

    const WeightCounter& wc() const {
        VERIFY(wc_);
        return *wc_;
    }

protected:
    void RemoveTrivial(const BidirectionalPath& path, std::set<size_t>& to_exclude) const {
        if (excludeTrivialWithBulges_) {
            analyzer_.ExcludeTrivialWithBulges(path, to_exclude);
        } else if (excludeTrivial_) {
            analyzer_.ExcludeTrivial(path, to_exclude);
        }
    }

private:
    DECL_LOGGER("ExtensionChooser");
};


class JointExtensionChooser: public ExtensionChooser {

protected:
    shared_ptr<ExtensionChooser> first_;

    shared_ptr<ExtensionChooser> second_;

public:
    JointExtensionChooser(Graph& g, shared_ptr<ExtensionChooser> first, shared_ptr<ExtensionChooser> second): ExtensionChooser(g),
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
    TrivialExtensionChooserWithPI(Graph& g, shared_ptr<WeightCounter> wc, double weight_threshold): 
            ExtensionChooser(g, wc, weight_threshold) {
    }

    virtual EdgeContainer Filter(BidirectionalPath& path, EdgeContainer& edges) {
        if (edges.size() == 1) {
            double weight = wc_->CountWeight(path, edges.back().e_, std::set<size_t>());
            NotifyAll(weight);

            if (CheckThreshold(weight)) {
                return edges;
            }
        }
        return EdgeContainer();
    }
};

class ExcludingExtensionChooser: public ExtensionChooser {

    double prior_coeff_;

    AlternativeContainer FindWeights(const BidirectionalPath& path, const EdgeContainer& edges, const std::set<size_t>& to_exclude) const {
        AlternativeContainer weights;
        for (auto iter = edges.begin(); iter != edges.end(); ++iter) {
            double weight = wc_->CountWeight(path, iter->e_, to_exclude);
            weights.insert(std::make_pair(weight, *iter));
            DEBUG("Candidate " << g_.int_id(iter->e_) << " weight " << weight << " length " << g_.length(iter->e_));
        }
        NotifyAll(weights);
        return weights;
    }

    EdgeContainer FindPossibleEdges(const AlternativeContainer& weights, 
            double max_weight) const {
        EdgeContainer top;
        auto possible_edge = weights.lower_bound(max_weight / prior_coeff_);
        for (auto iter = possible_edge; iter != weights.end(); ++iter) {
            top.push_back(iter->second);
        }
        return top;
    }

    EdgeContainer FindFilteredEdges(const BidirectionalPath& path,
            const EdgeContainer& edges, const std::set<size_t>& to_exclude) const {
        AlternativeContainer weights = FindWeights(path, edges, to_exclude);
        auto max_weight = (--weights.end())->first;
        EdgeContainer top = FindPossibleEdges(weights, max_weight);
        EdgeContainer result;
        if (top.size() >= 1 && CheckThreshold(max_weight)) {
            result = top;
        }
        return result;
    }

protected:

    virtual void ExcludeEdges(const BidirectionalPath& path, const EdgeContainer& edges, std::set<size_t>& to_exclude) const = 0;

public:
    ExcludingExtensionChooser(const Graph& g, shared_ptr<WeightCounter> wc, double weight_threshold, double priority) :
            ExtensionChooser(g, wc, weight_threshold), prior_coeff_(priority) {

    }

    virtual EdgeContainer Filter(const BidirectionalPath& path,
            const EdgeContainer& edges) const {
        DEBUG("Paired-end extension chooser");
        if (edges.empty()) {
            return edges;
        }
        std::set<size_t> to_exclude;
        RemoveTrivial(path, to_exclude);
        path.Print();
        EdgeContainer result = edges;
        ExcludeEdges(path, result, to_exclude);
        result = FindFilteredEdges(path, result, to_exclude);
        if (result.size() == 1) {
            DEBUG("Paired-end extension chooser helped");
        }
        return result;
    }

private:
    DECL_LOGGER("ExcludingExtensionChooser");

};

class SimpleExtensionChooser: public ExcludingExtensionChooser {
protected:
	virtual void ExcludeEdges(const BidirectionalPath& path, const EdgeContainer& edges, std::set<size_t>& to_exclude) const {
        if (edges.size() < 2) {
            return;
        }
        int index = (int) path.Size() - 1;
        while (index >= 0) {
            if (to_exclude.count(index)) {
                index--;
                continue;
            }
            EdgeId path_edge = path[index];

            bool common = true;
            for (size_t i = 0; i < edges.size(); ++i) {
                if (!wc_->PairInfoExist(path_edge, edges.at(i).e_,
                                        (int) path.LengthAt(index))) {
                    common = false;
                }
            }
            if (common) {
                DEBUG("common info from " << index);
                to_exclude.insert((size_t) index);
            } 

            index--;
        }
    }

public:

	SimpleExtensionChooser(const Graph& g, shared_ptr<WeightCounter> wc, double weight_threshold, double priority) :
	    ExcludingExtensionChooser(g, wc, weight_threshold, priority) {
	}

private:
    DECL_LOGGER("SimpleExtensionChooser");
};

class LongEdgeExtensionChooser: public ExcludingExtensionChooser {
protected:
    virtual void ExcludeEdges(const BidirectionalPath& path, const EdgeContainer& edges, std::set<size_t>& to_exclude) const {
        if (edges.size() < 2) {
            return;
        }
        int index = (int) path.Size() - 1;
        while (index >= 0) {
            if (to_exclude.count(index)) {
                index--;
                continue;
            }
            EdgeId path_edge = path[index];
            //FIXME configure!
            if (path.graph().length(path_edge) < 200)
                to_exclude.insert((size_t) index);
            index--;
        }
    }
public:
    LongEdgeExtensionChooser(const Graph& g, shared_ptr<WeightCounter> wc, double weight_threshold, double priority) :
        ExcludingExtensionChooser(g, wc, weight_threshold, priority) {
    }
};

class ScaffoldingExtensionChooser : public ExtensionChooser {
    double raw_weight_threshold_;
    double cl_weight_threshold_;
    const double is_scatter_coeff_ = 3.0;

    void AddInfoFromEdge(const std::vector<int>& distances, const std::vector<double>& weights, 
                         std::vector<pair<int, double>>& histogram, size_t len_to_path_end) const {
        for (size_t l = 0; l < distances.size(); ++l) {
            //todo commented out condition seems unnecessary and should be library dependent! do we need "max(0" there?
            if (/*distances[l] > max(0, (int) len_to_path_end - int(1000)) && */math::ge(weights[l], raw_weight_threshold_)) {
                histogram.push_back(make_pair(distances[l] - (int) len_to_path_end, weights[l]));
            }
        }
    }

    int CountMean(const vector<pair<int, double> >& histogram) const {
        double dist = 0.0;
        double sum = 0.0;
        for (size_t i = 0; i < histogram.size(); ++i) {
            dist += histogram[i].first * histogram[i].second;
            sum += histogram[i].second;
        }
        dist /= sum;
        return (int) round(dist);
    }

    void CountAvrgDists(const BidirectionalPath& path, EdgeId e, std::vector<pair<int, double>> & histogram) const {
        for (size_t j = 0; j < path.Size(); ++j) {
            std::vector<int> distances;
            std::vector<double> weights;
            wc_->GetDistances(path.At(j), e, distances, weights);
            if (distances.size() > 0) {
                AddInfoFromEdge(distances, weights, histogram, path.LengthAt(j));
            }
        }
    }

	double CountIdealInfo(const BidirectionalPath& p, EdgeId e, size_t gap) const {
		double w = 0.0;
		for (int i = (int) p.Size() - 1; i >= 0; --i) {
			w += wc_->CountIdealInfo(p[i], e, gap + p.LengthAt(i));
		}
		return w;
	}

    void FindBestFittedEdgesForClustered(const BidirectionalPath& path, const set<EdgeId>& edges, EdgeContainer& result) const {
        for (EdgeId e : edges) {
            std::vector<pair<int, double>> histogram;
            CountAvrgDists(path, e, histogram);
            double sum = 0.0;
            for (size_t j = 0; j < histogram.size(); ++j) {
                sum += histogram[j].second;
            }
            if (sum <= cl_weight_threshold_) {
                continue;
            }
            int gap = CountMean(histogram);
            if (CountIdealInfo(path, e, gap) > 0.0) {
                DEBUG("scaffolding " << g_.int_id(e) << " gap " << gap);
                result.push_back(EdgeWithDistance(e, gap));
            }
        }
    }

    bool IsTip(EdgeId e) const {
        return g_.IncomingEdgeCount(g_.EdgeStart(e)) == 0;
    }

    set<EdgeId> FindCandidates(const BidirectionalPath& path) const {
        set<EdgeId> jumping_edges;
        PairedInfoLibraries libs = wc_->getLibs();
        for (auto lib : libs) {
            //todo lib (and FindJumpEdges) knows its var so it can be counted there
            int is_scatter = int(math::round(double(lib->GetIsVar()) * is_scatter_coeff_));
            for (int i = (int) path.Size() - 1; i >= 0 && path.LengthAt(i) - g_.length(path.At(i)) <= lib->GetISMax(); --i) {
                set<EdgeId> jump_edges_i;
                lib->FindJumpEdges(path.At(i), jump_edges_i,
                                   std::max(0, (int)path.LengthAt(i) - is_scatter),
                                   //FIXME do we need is_scatter here?
                                   int((path.LengthAt(i) + lib->GetISMax() + is_scatter)),
                                   0);
                for (EdgeId e : jump_edges_i) {
                    if (IsTip(e)) {
                        jumping_edges.insert(e);
                    }
                }
            }
        }
        return jumping_edges;
    }

public:

    ScaffoldingExtensionChooser(const Graph& g, shared_ptr<WeightCounter> wc, double is_scatter_coeff) :
        ExtensionChooser(g, wc), raw_weight_threshold_(0.0),
        cl_weight_threshold_(cfg::get().pe_params.param_set.scaffolder_options.cl_threshold),
        is_scatter_coeff_(is_scatter_coeff) {
    }

    EdgeContainer Filter(const BidirectionalPath& path, const EdgeContainer& edges) const override {
        if (edges.empty()) {
            return edges;
        }
        set<EdgeId> candidates = FindCandidates(path);
        EdgeContainer result;
        FindBestFittedEdgesForClustered(path, candidates, result);
        return result;
    }
private:
    DECL_LOGGER("ScaffoldingExtensionChooser");
};

inline bool EdgeWithWeightCompareReverse(const pair<EdgeId, double>& p1,
                                      const pair<EdgeId, double>& p2) {
    return p1.second > p2.second;
}

class UniqueEdgeAnalyzer {
public:
    UniqueEdgeAnalyzer(const Graph& g, const GraphCoverageMap& cov_map,
                       double filter_threshold, double prior_threshold)
            : g_(g),
              cov_map_(cov_map),
              filter_threshold_(filter_threshold),
              prior_threshold_(prior_threshold) { 
        FindAllUniqueEdges();
    }

    bool IsUnique(EdgeId e) const {
        return unique_edges_.count(e) > 0;
    }

private:
    bool UniqueEdge(EdgeId e) const {
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
            const BidirectionalPathSet& cov_paths) const {
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
            const BidirectionalPathSet& cov_paths) const {
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

    void FindAllUniqueCoverageEdges() {
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
       DEBUG("coverage based uniqueness started");
       FindAllUniqueCoverageEdges();
       DEBUG("Unique edges are found");
    }

    const Graph& g_;
    const GraphCoverageMap& cov_map_;
    double filter_threshold_;
    double prior_threshold_;
    std::set<EdgeId> unique_edges_;

    DECL_LOGGER("UniqueEdgeAnalyzer")
};

class SimpleScaffolding {
public:
    SimpleScaffolding(const Graph& g) : g_(g) {}

    BidirectionalPath FindMaxCommonPath(const vector<BidirectionalPath*>& paths,
                                        size_t max_diff_len) const {
        BidirectionalPath max_end(g_);
        for (auto it1 = paths.begin(); it1 != paths.end(); ++it1) {
            BidirectionalPath* p1 = *it1;
            for (size_t i = 0; i < p1->Size(); ++i) {
                if (p1->Length() - p1->LengthAt(i) > max_diff_len) {
                    break;
                }
                bool contain_all = true;
                for (size_t i1 = i + 1; i1 <= p1->Size() && contain_all; ++i1) {
                    BidirectionalPath subpath = p1->SubPath(i, i1);
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
                            }
                        }
                        if (!contain) {
                            contain_all = false;
                        }
                    }
                    if (contain_all && (i1 - i) >= max_end.Size()) {
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
            : ExtensionChooser(g),
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
    EdgeContainer Filter(const BidirectionalPath& path,
                                 const EdgeContainer& edges) const override {
        if (edges.empty()) {
            return edges;
        }DEBUG("We in Filter of LongReadsExtensionChooser");
        path.Print();
        map<EdgeId, double> weights_cands;
        for (auto it = edges.begin(); it != edges.end(); ++it) {
            weights_cands.insert(make_pair(it->e_, 0.0));
        }
        set<EdgeId> filtered_cands;
        map<EdgeId, BidirectionalPathSet > support_paths_ends;
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
                            support_paths_ends[next] = BidirectionalPathSet();
                        }
                        support_paths_ends[next].insert(new BidirectionalPath((*it)->SubPath(positions[i] + 1)));
                    }
                }
            }
        }
        DEBUG("Candidates");
        for (auto iter = weights_cands.begin(); iter != weights_cands.end(); ++iter) {
            DEBUG("Candidate " << g_.int_id(iter->first) << " weight " << iter->second);
        }
        vector<pair<EdgeId, double> > sort_res = MapToSortVector(weights_cands);
        DEBUG("sort res " << sort_res.size() << " tr " << weight_priority_threshold_);
        if (sort_res.size() < 1 || sort_res[0].second < filtering_threshold_) {
            filtered_cands.clear();
        } else if (sort_res.size() > 1
                && sort_res[0].second > weight_priority_threshold_ * sort_res[1].second) {
            filtered_cands.clear();
            filtered_cands.insert(sort_res[0].first);
        } else if (sort_res.size() > 1) {
            for (size_t i = 0; i < sort_res.size(); ++i) {
                if (sort_res[i].second * weight_priority_threshold_ < sort_res[0].second) {
                    filtered_cands.erase(sort_res[i].first);
                }
            }
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
    bool UniqueBackPath(const BidirectionalPath& path, size_t pos) const {
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

    DECL_LOGGER("LongReadsExtensionChooser");
};

class MatePairExtensionChooser : public ExtensionChooser {
public:
    MatePairExtensionChooser(const Graph& g, shared_ptr<PairedInfoLibrary> lib,
                              const PathContainer& paths, size_t max_number_of_paths_to_search)
            : ExtensionChooser(g),
              g_(g),
              lib_(lib),
              search_dist_(lib->GetISMax()),
              weight_counter_(g, lib, 10),
              cov_map_(g_, paths),
              path_searcher_(g_, cov_map_, lib_->GetISMax(), PathsWeightCounter(g, lib, 30), max_number_of_paths_to_search),
              unique_edge_analyzer_(g, cov_map_, 0., 1000.),
              simple_scaffolder_(g) {
    }

    //Attention! Uses const_cast to modify path!!!
    EdgeContainer Filter(const BidirectionalPath& path,
                         const EdgeContainer& init_edges) const override {
        DEBUG("mp chooser");
        path.Print();
        if (path.Length() < lib_->GetISMin()) {
            return EdgeContainer();
        }
        EdgeContainer edges = TryResolveBulge(path, init_edges);
        map<EdgeId, BidirectionalPath*> best_paths;
        for (size_t iedge = 0; iedge < edges.size(); ++iedge) {
            BidirectionalPathSet following_paths = path_searcher_.FindNextPaths(path, edges[iedge].e_);
            vector<BidirectionalPath*> max_weighted = MaxWeightedPath(path, following_paths);
            if (max_weighted.size() == 0) {
                DEBUG("too much paths or tip");
                DeleteMapWithPaths(best_paths);
                DeletePaths(following_paths);
                best_paths.clear();
                break;
            } else {
                best_paths[edges[iedge].e_] = new BidirectionalPath(*max_weighted[0]);
            }
            DeletePaths(following_paths);
        }

        BidirectionalPathSet next_paths;
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
            result = ScaffoldTree(const_cast<BidirectionalPath&>(path));
        }
        DeletePaths(next_paths);
        if (result.size() != 1) {
            DEBUG("nobody can extend " << g_.int_id(path.Back()));
        }
        return result;
    }

private:
    EdgeContainer ScaffoldTree(BidirectionalPath& path) const {
        DEBUG("try scaffold tree");
        vector<BidirectionalPath*> next_paths = path_searcher_.ScaffoldTree(path);
        VERIFY(next_paths.size() <= 1);
        EdgeContainer result;
        if (!next_paths.empty() && next_paths.back()->Size() > 0) {
            BidirectionalPath* res = next_paths.back();
            for (size_t i = 0; i < res->Size() - 1; ++i) {
                path.PushBack(res->At(i), res->GapAt(i));
            }
            result = EdgeContainer(1, EdgeWithDistance(res->Back(), res->GapAt(res->Size() - 1)));
        }
        DeletePaths(next_paths);
        return result;
    }

	bool IsBulge(const EdgeContainer& edges) const {
		if (edges.size() == 0)
			return false;
		for (EdgeWithDistance e : edges) {
			if (!InBuble(e.e_, g_))
				return false;
		}
		return true;
	}

    map<EdgeId, double> FindBulgeWeights(const BidirectionalPath& p, const EdgeContainer& edges) const {
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

    EdgeContainer TryResolveBulge(const BidirectionalPath& p, const EdgeContainer& edges) const {
        if (!IsBulge(edges))
            return edges;
        map<EdgeId, double> weights = FindBulgeWeights(p, edges);
        double max_w = 0.0;
        EdgeContainer result;
        for (EdgeWithDistance e : edges) {
            double w = weights[e.e_];
            DEBUG("bulge " << g_.int_id(e.e_) << " w = " << w);
            if (math::gr(w, max_w)) {
                max_w = w;
                result.clear();
                result.push_back(e);
            } else if (math::eq(w, max_w)) {
                result.push_back(e);
            }
        }
        if (result.size() != 1) {
            result = edges;
        }
        return result;
    }

    EdgeContainer ChooseBest(const BidirectionalPath& path, const BidirectionalPathSet& next_paths) const {
        DEBUG("Try to choose from best paths...");
        vector<BidirectionalPath*> best_path = MaxWeightedPath(path, next_paths);
        EdgeContainer result;
        if (best_path.size() == 1) {
            result.push_back(EdgeWithDistance((*best_path.begin())->At(0), (*best_path.begin())->GapAt(0)));
        } else if (best_path.size() > 1) {
            result = TryToScaffold(path, best_path);
        }
        return result;
    }

    bool HasPIFromUniqueEdges(const BidirectionalPath& p1, const BidirectionalPath& p2, const set<size_t>& p1_unique_edges) const {
        for (size_t i1 = 0; i1 < p1.Size(); ++i1) {
            if (p1_unique_edges.find(i1) == p1_unique_edges.end()) {
                continue;
            }
            for (size_t i2 = 0; i2 < p2.Size(); ++i2) {
                int gap = (int) p1.LengthAt(i1) + (int) p2.Length() - (int) p2.LengthAt(i2);
                if (unique_edge_analyzer_.IsUnique(p2.At(i2)) && weight_counter_.HasPI(p1.At(i1), p2.At(i2), gap)) {
                    DEBUG("has unique edge " << g_.int_id(p1.At(i1)) << " " << g_.int_id(p2.At(i2)));
                    return true;
                }
            }
        }
        return false;
    }

    bool SignificallyDifferentEdges(const BidirectionalPath& init_path, const BidirectionalPath& path1, const map<size_t, double>& pi1,
                                    const BidirectionalPath& path2, const map<size_t, double>& pi2, const set<size_t>& unique_init_edges) const {
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
        if (common_w < 0.8 * (not_common_w1 + common_w)
                || (HasPIFromUniqueEdges(init_path, path1, unique_init_edges) && !HasPIFromUniqueEdges(init_path, path2, unique_init_edges))) {
            DEBUG("common_w " << common_w << " sum * 0.8  = " << 0.8 * (not_common_w1 + common_w))
            return true;
        }
        return false;
    }

    set<size_t> FindNotCommonEdges(const BidirectionalPath& path, const BidirectionalPathMap< map<size_t, double> >& all_pi) const {
        set<size_t> res;
        for (size_t i = 0; i < path.Size(); ++i) {
            if (!unique_edge_analyzer_.IsUnique(path.At(i))) {
                continue;
            }
            size_t pi_count = 0;
            for (auto iter = all_pi.begin(); iter != all_pi.end(); ++iter) {
                const map<size_t, double>& info = iter->second;
                if (info.count(i) > 0 && math::gr(info.at(i), 0.0)) {
                    pi_count++;
                }
            }
            if (pi_count == 1)
                res.insert(i);
        }
        return res;
    }

    void DeleteSmallWeights(const BidirectionalPath& path, BidirectionalPathSet& paths, BidirectionalPathMap< map<size_t, double> >& all_pi) const {
        double max_weight = 0.0;
        BidirectionalPath* max_path = NULL;
        for (auto iter = paths.begin(); iter != paths.end(); ++iter) {
            if ((*iter)->GetWeight() >= max_weight) {
                max_weight = max(max_weight, (*iter)->GetWeight());
                max_path = *iter;
            }
        }
        BidirectionalPathSet to_del;
        for (auto iter = paths.begin(); iter != paths.end(); ++iter) {
            if (math::gr(max_weight, (*iter)->GetWeight() * 1.5) //TODO: move 1.5 to config
                    && SignificallyDifferentEdges(path, *max_path, all_pi.find(max_path)->second, **iter, all_pi.find(*iter)->second,
                                                  FindNotCommonEdges(path, all_pi)))
                to_del.insert(*iter);
        }
        for (BidirectionalPath* p : to_del) {
            paths.erase(p);
            all_pi.erase(p);
        }
    }

    void DeleteCommonPi(const BidirectionalPath& p, BidirectionalPathMap< map<size_t, double> >& all_pi) const {
        weight_counter_.ClearCommonWeight();
        for (size_t i = 0; i < p.Size(); ++i) {
            double common = DBL_MAX;
            for (auto iter = all_pi.begin(); iter != all_pi.end(); ++iter) {
                common = iter->second.count(i) == 0 ? 0.0 : min(common, iter->second.at(i));
            }
            weight_counter_.SetCommonWeightFrom(i, common);
        }
    }

    size_t FindCommonBegin(const BidirectionalPathSet& paths) const {
        if (paths.size() == 0) {
            return 0;
        }
        size_t common_begin = 0;
        BidirectionalPath* p = *paths.begin();
        while (common_begin < p->Size()) {
            EdgeId e = p->At(common_begin);
            for (BidirectionalPath* next : paths) {
                if (common_begin >= next->Size() || next->At(common_begin) != e) {
                    return common_begin;
                }
            }
            common_begin++;
        }
        return common_begin;
    }

    void CountAllPairInfo(const BidirectionalPath& path, const BidirectionalPathSet& next_paths,
                BidirectionalPathMap<map<size_t, double>>& result) const {
        result.clear();
        size_t common_begin = FindCommonBegin(next_paths);
        DEBUG("common begin " << common_begin);
        for (BidirectionalPath* next : next_paths) {
            result[next] = weight_counter_.FindPairInfoFromPath(path, 0, path.Size(), *next, common_begin, next->Size());
        }
    }

    void CountWeightsAndFilter(const BidirectionalPath& path, BidirectionalPathSet& next_paths, bool delete_small_w) const {
        BidirectionalPathMap<map<size_t, double> > all_pi;
        CountAllPairInfo(path, next_paths, all_pi);
        DeleteCommonPi(path, all_pi);
        for (BidirectionalPath* next : next_paths) {
            next->SetWeight((float) weight_counter_.CountPairInfo(path, 0, path.Size(), *next, 0, next->Size()));
        }
        if (delete_small_w) {
            DeleteSmallWeights(path, next_paths, all_pi);
        }
    }

    struct PathWithWeightSort {
        PathWithWeightSort(const MatePairExtensionChooser& mp_chooser, const BidirectionalPath& path, BidirectionalPathMap< map<size_t, double> >& all_pi)
                : mp_chooser_(mp_chooser),
                  path_(path),
                  not_common_(mp_chooser_.FindNotCommonEdges(path_, all_pi)) {
        }

        bool operator()(const BidirectionalPath* p1, const BidirectionalPath* p2) {
            if (mp_chooser_.HasPIFromUniqueEdges(path_, *p1, not_common_) && !mp_chooser_.HasPIFromUniqueEdges(path_, *p2, not_common_)) {
                return true;
            }
            if (mp_chooser_.HasPIFromUniqueEdges(path_, *p2, not_common_) && !mp_chooser_.HasPIFromUniqueEdges(path_, *p1, not_common_)) {
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
        const MatePairExtensionChooser& mp_chooser_;
        const BidirectionalPath& path_;
        const set<size_t> not_common_;
    };

    vector<BidirectionalPath*> SortResult(const BidirectionalPath& path, BidirectionalPathSet& next_paths) const {
        BidirectionalPathMap< map<size_t, double> > all_pi;
        CountAllPairInfo(path, next_paths, all_pi);
        CountWeightsAndFilter(path, next_paths, false);
        vector<BidirectionalPath*> to_sort(next_paths.begin(), next_paths.end());
        PathWithWeightSort comparator(*this, path, all_pi);
        std::sort(to_sort.begin(), to_sort.end(), comparator);
        return to_sort;
    }

    vector<BidirectionalPath*> MaxWeightedPath(const BidirectionalPath& path, const BidirectionalPathSet& following_paths) const {
        BidirectionalPathSet result(following_paths);
        BidirectionalPathSet prev_result;
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

    BidirectionalPath ChooseFromEnds(const BidirectionalPath& path, const vector<BidirectionalPath*>& paths, const BidirectionalPath& end) const { //TODO" rewrite
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
        for (BidirectionalPath* p : new_paths) {
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
        for (BidirectionalPath* p : new_paths) {
            delete p;
        }
        for (BidirectionalPath* p : paths_to_cover) {
            delete p;
        }
        BidirectionalPath result = max.SubPath(path.Size());
        DEBUG("res");
        result.Print();
        return result;
    }

    int CheckPairInfo(const BidirectionalPath& path, const BidirectionalPath& result_end, int to_add) const {
        while (to_add < (int)result_end.Size()) {
            map<size_t, double> weights = weight_counter_.FindPairInfoFromPath(path, 0, path.Size(), result_end, to_add, to_add + 1);
            double weight_to_edge = 0.0;
            for (auto iter = weights.begin(); iter != weights.end(); ++iter) {
                weight_to_edge += iter->second;
            }
            if (math::gr(weight_to_edge, 0.0)) {
                break;
            }
            to_add++;
        }
        return to_add;
    }

    EdgeContainer TryToScaffold(const BidirectionalPath& path, const vector<BidirectionalPath*>& paths) const {
        if (paths.size() == 0) {
            return EdgeContainer();
        }
        DEBUG("Simple Scaffolding")
        for (BidirectionalPath* p : paths) {
            p->Print();
        }
        BidirectionalPath max_end = simple_scaffolder_.FindMaxCommonPath(paths, search_dist_);
        if (max_end.Size() == 0) {
            return EdgeContainer();
        }
        BidirectionalPath result_end = ChooseFromEnds(path, paths, max_end);
        int to_add = result_end.FindFirst(max_end);
        result_end.Print();
        EdgeContainer result;
        to_add = CheckPairInfo(path, result_end, to_add);
        if (to_add < 0 || to_add >= (int) result_end.Size()) {
            return EdgeContainer();
        }
        size_t gap_length = result_end.Length() - result_end.LengthAt(to_add);
        DEBUG(" edge to add " << g_.int_id(result_end.At(to_add)) << " with length " << gap_length);
        result.push_back(EdgeWithDistance(result_end.At(to_add), gap_length));
        return result;
    }

    const Graph& g_;
    shared_ptr<PairedInfoLibrary> lib_;
    size_t search_dist_;
    mutable PathsWeightCounter weight_counter_;
    const GraphCoverageMap cov_map_;
    NextPathSearcher path_searcher_;
    UniqueEdgeAnalyzer unique_edge_analyzer_;
    SimpleScaffolding simple_scaffolder_;

    DECL_LOGGER("MatePairExtensionChooser");
};
}
#endif /* EXTENSION_HPP_ */

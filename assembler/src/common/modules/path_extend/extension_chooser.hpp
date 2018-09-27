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
#include <map>
#include "weight_counter.hpp"
#include "pe_utils.hpp"

//#include "scaff_supplementary.hpp"

namespace path_extend {

typedef std::multimap<double, EdgeWithDistance> AlternativeContainer;

class PathAnalyzer {
protected:
    const Graph& g_;

public:
    PathAnalyzer(const Graph& g): g_(g) {
    }

    void RemoveTrivial(const BidirectionalPath& path, std::set<size_t>& to_exclude, bool exclude_bulges = true) const {
        if (exclude_bulges) {
            ExcludeTrivialWithBulges(path, to_exclude);
        } else {
            ExcludeTrivial(path, to_exclude);
        }
    }

protected:
    virtual int ExcludeTrivial(const BidirectionalPath& path, std::set<size_t>& edges, int from = -1) const {
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

    virtual int ExcludeTrivialWithBulges(const BidirectionalPath& path, std::set<size_t>& edges) const {

        if (path.Empty()) {
            return 0;
        }

        int lastEdge = (int) path.Size() - 1;
        do {
            lastEdge = ExcludeTrivial(path, edges, lastEdge);
            bool bulge = true;

            if (lastEdge >= 0) {
                VertexId v = g_.EdgeEnd(path[lastEdge]);
                VertexId u = g_.EdgeStart(path[lastEdge]);
                auto bulgeCandidates = g_.IncomingEdges(v);

                for (const auto& candidate: bulgeCandidates) {
                    if (g_.EdgeStart(candidate) != u) {
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


class PreserveSimplePathsAnalyzer: public PathAnalyzer {

public:
    PreserveSimplePathsAnalyzer(const Graph &g) : PathAnalyzer(g) {
    }

    int ExcludeTrivial(const BidirectionalPath& path, std::set<size_t>& edges, int from = -1) const override {
        int edgeIndex = PathAnalyzer::ExcludeTrivial(path, edges, from);

        //Preserving simple path
        if (edgeIndex == -1) {
            edges.clear();
            return (from == -1) ? (int) path.Size() - 1 : from;;
        }
        return edgeIndex;
    }

    int ExcludeTrivialWithBulges(const BidirectionalPath& path, std::set<size_t>& edges) const override {

        if (path.Empty()) {
            return 0;
        }

        int lastEdge = (int) path.Size() - 1;
        bool has_bulge = false;
        do {
            lastEdge = PathAnalyzer::ExcludeTrivial(path, edges, lastEdge);

            if (lastEdge >= 0) {
                VertexId v = g_.EdgeEnd(path[lastEdge]);
                VertexId u = g_.EdgeStart(path[lastEdge]);
                auto bulgeCandidates = g_.IncomingEdges(v);
                has_bulge = true;

                for (auto iter = bulgeCandidates.begin(); iter != bulgeCandidates.end(); ++iter) {
                    if (g_.EdgeStart(*iter) != u) {
                        has_bulge = false;
                        break;
                    }
                }

                --lastEdge;
            }
        } while (lastEdge >= 0);

        //Preserving simple path
        if (!has_bulge && lastEdge == -1) {
            edges.clear();
            lastEdge = (int) path.Size() - 1;
        }

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
    //FIXME memory leak?!
    std::vector<ExtensionChooserListener *> listeners_;

    double weight_threshold_;

public:
    ExtensionChooser(const Graph& g, shared_ptr<WeightCounter> wc = nullptr, double weight_threshold = -1.):
        g_(g), wc_(wc), 
        weight_threshold_(weight_threshold) {
    }

    virtual ~ExtensionChooser() {

    }

    virtual EdgeContainer Filter(const BidirectionalPath& path, const EdgeContainer& edges) const = 0;

    bool CheckThreshold(double weight) const {
        return math::ge(weight, weight_threshold_);
    }

    void Subscribe(ExtensionChooserListener * listener) {
        listeners_.push_back(listener);
    }

    void NotifyAll(double weight) const {
        for (auto listener_ptr : listeners_) {
            listener_ptr->ExtensionChosen(weight);
        }
    }

    void NotifyAll(const AlternativeContainer& alts) const {
        for (auto listener_ptr : listeners_) {
            listener_ptr->ExtensionChosen(alts);
        }
    }

    bool WeightCounterBased() const {
        return wc_ != nullptr;
    }

    shared_ptr<WeightCounter> wc() const {
        return wc_;
    }

protected:
    bool HasIdealInfo(EdgeId e1, EdgeId e2, size_t dist) const {
        return math::gr(wc_->PairedLibrary().IdealPairedInfo(e1, e2, (int) dist), 0.);
    }

    bool HasIdealInfo(const BidirectionalPath& p, EdgeId e, size_t gap) const {
        for (int i = (int) p.Size() - 1; i >= 0; --i)
            if (HasIdealInfo(p[i], e, gap + p.LengthAt(i)))
                return true;
        return false;
    }

private:
    DECL_LOGGER("ExtensionChooser");
};


class JointExtensionChooser: public ExtensionChooser {
    shared_ptr<ExtensionChooser> first_;
    shared_ptr<ExtensionChooser> second_;

public:
    JointExtensionChooser(const Graph& g,
                          shared_ptr<ExtensionChooser> first,
                          shared_ptr<ExtensionChooser> second): ExtensionChooser(g),
        first_(first), second_(second) {
    }

    EdgeContainer Filter(const BidirectionalPath& path, const EdgeContainer& edges) const override {
        EdgeContainer answer;
        auto r1 = first_->Filter(path, edges);
        auto r2 = second_->Filter(path, edges);
        for (auto ewd1 : r1) {
            for (auto ewd2 : r2) {
                if (ewd1.e_ == ewd2.e_) {
                    VERIFY(ewd1.d_ == ewd2.d_);
                    answer.push_back(ewd1);
                }
            }
        }
        return answer;
    }
};


class TrivialExtensionChooser: public ExtensionChooser {

public:
    TrivialExtensionChooser(Graph& g): ExtensionChooser(g)  {
    }

    EdgeContainer Filter(const BidirectionalPath& /*path*/, const EdgeContainer& edges) const override {
        if (edges.size() == 1) {
             return edges;
        }
        return EdgeContainer();
    }
};


class SimpleCoverageExtensionChooser: public ExtensionChooser {
    const SSCoverageStorage& coverage_storage_;
    //less than 1
    double coverage_delta_;
    //larger than 1
    double inverted_coverage_delta_;

    double min_upper_coverage_;

public:
    SimpleCoverageExtensionChooser(const SSCoverageStorage& coverage_storage, const Graph& g,
                                   double coverage_delta, double min_upper_coverage = 0) :
        ExtensionChooser(g), coverage_storage_(coverage_storage),
        coverage_delta_(coverage_delta),
        inverted_coverage_delta_(0),
        min_upper_coverage_(min_upper_coverage) {
        VERIFY(math::le(coverage_delta_, 1.0));
        VERIFY(!math::eq(coverage_delta_, 0.0));
        inverted_coverage_delta_ = 1.0 / coverage_delta_;
    }

    EdgeContainer Filter(const BidirectionalPath& path, const EdgeContainer& edges) const override {
        if (edges.size() != 2)
            return EdgeContainer();

        size_t index = path.Size() - 1;
        while (index > 0) {
            if (g_.IncomingEdgeCount(g_.EdgeStart(path[index])) == 2)
                break;
            index--;
        }

        if (index == 0) {
            return EdgeContainer();
        }
        DEBUG("Split found at " << index);
        EdgeId path_edge_at_split = path[index - 1];

        return Filter(path, edges, math::ls(coverage_storage_.GetCoverage(path_edge_at_split), coverage_storage_.GetCoverage(path_edge_at_split, true)));
    }

private:
    EdgeContainer Filter(const BidirectionalPath& path, const EdgeContainer& edges, bool reverse) const {
        DEBUG("COVERAGE extension chooser");
        VERIFY(edges.size() == 2);
        if (!IsEnoughCoverage(edges.front().e_, edges.back().e_, reverse)) {
            DEBUG("Candidates are not covered enough: e1 = " << coverage_storage_.GetCoverage(edges.front().e_, reverse) <<
                ", e2 = " << coverage_storage_.GetCoverage(edges.back().e_, reverse));
            return EdgeContainer();
        }

        if (IsCoverageSimilar(edges.front().e_, edges.back().e_, reverse)) {
            DEBUG("Candidates coverage is too similar: e1 = " << coverage_storage_.GetCoverage(edges.front().e_, reverse) <<
                ", e2 = " << coverage_storage_.GetCoverage(edges.back().e_, reverse));
            return EdgeContainer();
        }

        size_t index = path.Size() - 1;
        while (index > 0) {
            if (g_.IncomingEdgeCount(g_.EdgeStart(path[index])) == 2)
                break;
            index--;
        }

        EdgeContainer result;
        if (index > 0) {
            DEBUG("Split found at " << index);
            EdgeId path_edge_at_split = path[index - 1];
            EdgeId other_edge_at_split = GetOtherEdgeAtSplit(g_.EdgeEnd(path_edge_at_split), path_edge_at_split);
            VERIFY(other_edge_at_split != EdgeId());

            if (IsCoverageSimilar(path_edge_at_split, other_edge_at_split, reverse)) {
                DEBUG("Path edge and alternative is too similar: path = " << coverage_storage_.GetCoverage(path_edge_at_split, reverse) <<
                    ", other = " << coverage_storage_.GetCoverage(other_edge_at_split, reverse));

                return EdgeContainer();
            }
            if (!IsEnoughCoverage(path_edge_at_split, other_edge_at_split, reverse)) {
                DEBUG("Path edge and alternative  coverage is too low: path = " << coverage_storage_.GetCoverage(path_edge_at_split, reverse) <<
                    ", other = " << coverage_storage_.GetCoverage(other_edge_at_split, reverse));

                return EdgeContainer();
            }

            EdgeId candidate1 = edges.front().e_;
            EdgeId candidate2 = edges.back().e_;

            if (math::gr(coverage_storage_.GetCoverage(path_edge_at_split, reverse), coverage_storage_.GetCoverage(other_edge_at_split, reverse))) {
                DEBUG("path coverage is high, edge " << g_.int_id(path_edge_at_split) << ", path cov = "
                          << coverage_storage_.GetCoverage(path_edge_at_split, reverse) << ", other " << coverage_storage_.GetCoverage(other_edge_at_split, reverse));

                result.emplace_back(math::gr(coverage_storage_.GetCoverage(candidate1, reverse), coverage_storage_.GetCoverage(candidate2, reverse)) ? candidate1 : candidate2, 0);
            } else {
                DEBUG("path coverage is low, edge " << g_.int_id(path_edge_at_split) << ", path cov = "
                          << coverage_storage_.GetCoverage(path_edge_at_split, reverse) << ", other " << coverage_storage_.GetCoverage(other_edge_at_split, reverse));

                result.emplace_back(math::ls(coverage_storage_.GetCoverage(candidate1, reverse), coverage_storage_.GetCoverage(candidate2, reverse)) ? candidate1 : candidate2, 0);
            }

            if (!IsCoverageSimilar(path_edge_at_split, result.front().e_, reverse)) {
                DEBUG("Coverage is NOT similar: path = " << coverage_storage_.GetCoverage(path_edge_at_split, reverse) <<
                    ", candidate = " << coverage_storage_.GetCoverage(result.front().e_, reverse))
                result.clear();
            }
            else {
                DEBUG("Coverage is similar: path = " << coverage_storage_.GetCoverage(path_edge_at_split, reverse) <<
                    ", candidate = " << coverage_storage_.GetCoverage(result.front().e_, reverse))
                DEBUG("Coverage extension chooser helped, adding " << g_.int_id(result.front().e_));
            }
        }

        VERIFY(result.size() <= 1);
        return result;
    }

    bool IsEnoughCoverage(EdgeId e1, EdgeId e2, bool reverse) const {
        double cov1 = coverage_storage_.GetCoverage(e1, reverse);
        double cov2 = coverage_storage_.GetCoverage(e2, reverse);
        return math::ge(max(cov1, cov2), min_upper_coverage_) || math::eq(min(cov1, cov2), 0.0);
    }

    bool IsCoverageSimilar(EdgeId e1, EdgeId e2, bool reverse) const {
        double cov1 = coverage_storage_.GetCoverage(e1, reverse);
        double cov2 = coverage_storage_.GetCoverage(e2, reverse);

        if (math::eq(cov2, 0.0) || math::eq(cov1, 0.0)) {
            return false;
        }

        double diff = cov1 / cov2;
        if (math::ls(diff, 1.0))
            return math::gr(diff, coverage_delta_);
        else
            return math::ls(diff, inverted_coverage_delta_);
    }

    EdgeId GetOtherEdgeAtSplit(VertexId split, EdgeId e) const {
        VERIFY(g_.IncomingEdgeCount(split) == 2);
        for (auto other : g_.IncomingEdges(split)) {
            if (e != other)
                return other;
        }
        return EdgeId();
    }

    DECL_LOGGER("SimpleCoverageExtensionChooser");

};



class ExcludingExtensionChooser: public ExtensionChooser {
    PathAnalyzer analyzer_;
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
        VERIFY(!weights.empty());
        auto max_weight = (--weights.end())->first;
        EdgeContainer top = FindPossibleEdges(weights, max_weight);
        DEBUG("Top-scored edges " << top.size());
        EdgeContainer result;
        if (CheckThreshold(max_weight)) {
            result = top;
        }
        return result;
    }

protected:

    virtual void ExcludeEdges(const BidirectionalPath& path,
                              const EdgeContainer& /*edges*/,
                              std::set<size_t>& to_exclude) const {
        analyzer_.RemoveTrivial(path, to_exclude);
    }


public:
    ExcludingExtensionChooser(const Graph& g, shared_ptr<WeightCounter> wc, PathAnalyzer analyzer, double weight_threshold, double priority) :
            ExtensionChooser(g, wc, weight_threshold), analyzer_(analyzer), prior_coeff_(priority) {

    }

    virtual EdgeContainer Filter(const BidirectionalPath& path,
            const EdgeContainer& edges) const {
        DEBUG("Paired-end extension chooser");
        if (edges.empty()) {
            return edges;
        }
        std::set<size_t> to_exclude;
        path.PrintDEBUG();
        EdgeContainer result = edges;
        ExcludeEdges(path, result, to_exclude);
        DEBUG("Excluded " << to_exclude.size() << " edges")
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
    void ExcludeEdges(const BidirectionalPath& path, const EdgeContainer& edges, std::set<size_t>& to_exclude) const override {
        ExcludingExtensionChooser::ExcludeEdges(path, edges, to_exclude);

        if (edges.size() < 2) {
            return;
        }
        //excluding based on absence of ideal info
        int index = (int) path.Size() - 1;
        while (index >= 0) {
            if (to_exclude.count(index)) {
                index--;
                continue;
            }
            EdgeId path_edge = path[index];

            for (size_t i = 0; i < edges.size(); ++i) {
                if (!HasIdealInfo(path_edge,
                           edges.at(i).e_,
                           path.LengthAt(index))) {
                    DEBUG("Excluding edge because of no ideal info #" << index)
                    to_exclude.insert((size_t) index);
                }
            }

            index--;
        }
        
        //excluding based on presense of ambiguous paired info
        map<size_t, unsigned> edge_2_extension_cnt;
        for (size_t i = 0; i < edges.size(); ++i) {
            for (size_t e : wc_->PairInfoExist(path, edges.at(i).e_)) {
                edge_2_extension_cnt[e] += 1;
            }
        }

        for (auto e_w_ec : edge_2_extension_cnt) {
            if (e_w_ec.second == edges.size()) {
                DEBUG("Excluding edge because of ambiguous paired info #" << e_w_ec.first)
                to_exclude.insert(e_w_ec.first);
            }
        }
    }

public:

    SimpleExtensionChooser(const Graph& g, shared_ptr<WeightCounter> wc, double weight_threshold, double priority) :
        ExcludingExtensionChooser(g, wc, PathAnalyzer(g), weight_threshold, priority) {
    }

private:
    DECL_LOGGER("SimpleExtensionChooser");
};

//TODO this class should not exist with better configuration of excluding conditions
class IdealBasedExtensionChooser : public ExcludingExtensionChooser {
protected:
    void ExcludeEdges(const BidirectionalPath &path, const EdgeContainer &edges,
                      std::set<size_t> &to_exclude) const override {
        //commented for a reason
        //ExcludingExtensionChooser::ExcludeEdges(path, edges, to_exclude);
        //if (edges.size() < 2) {
        //    return;
        //}
        VERIFY(to_exclude.empty());
        //excluding based on absence of ideal info
        for (int index = (int) path.Size() - 1; index >= 0; index--) {
            EdgeId path_edge = path[index];

            for (size_t i = 0; i < edges.size(); ++i) {
                if (!HasIdealInfo(path_edge,
                                  edges.at(i).e_,
                                  path.LengthAt(index))) {
                    to_exclude.insert(size_t(index));
                }
            }
        }
    }

public:

    IdealBasedExtensionChooser(const Graph &g,
                               shared_ptr<WeightCounter> wc,
                               double weight_threshold,
                               double priority) :
        ExcludingExtensionChooser(g, wc, PathAnalyzer(g), weight_threshold, priority) {
    }

private:
    DECL_LOGGER("IdealBasedExtensionChooser");
};

class RNAExtensionChooser: public ExcludingExtensionChooser {
protected:
    void ExcludeEdges(const BidirectionalPath& path, const EdgeContainer& edges, std::set<size_t>& to_exclude) const override {
        ExcludingExtensionChooser::ExcludeEdges(path, edges, to_exclude);
        if (edges.size() < 2) {
            return;
        }
        size_t i = path.Size() - 1;
        PathAnalyzer analyzer(g_);
        while (i > 0) {
            if (g_.IncomingEdgeCount(g_.EdgeStart(path[i])) > 1)
                break;
            to_exclude.insert(i);
            --i;
            }

        if (i == 0)
            to_exclude.clear();
    }

public:

    RNAExtensionChooser(const Graph& g, shared_ptr<WeightCounter> wc, double weight_threshold, double priority) :
        ExcludingExtensionChooser(g, wc, PreserveSimplePathsAnalyzer(g), weight_threshold, priority) {
    }

private:
    DECL_LOGGER("SimpleExtensionChooser");
};

class LongEdgeExtensionChooser: public ExcludingExtensionChooser {
protected:
    virtual void ExcludeEdges(const BidirectionalPath& path, const EdgeContainer& edges, std::set<size_t>& to_exclude) const {
        ExcludingExtensionChooser::ExcludeEdges(path, edges, to_exclude);
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
        ExcludingExtensionChooser(g, wc, PathAnalyzer(g), weight_threshold, priority) {
    }
};

class ScaffoldingExtensionChooser : public ExtensionChooser {
    typedef ExtensionChooser base;
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

    void GetDistances(EdgeId e1, EdgeId e2, std::vector<int>& dist,
            std::vector<double>& w) const {
        wc_->PairedLibrary().CountDistances(e1, e2, dist, w);
    }

    void CountAvrgDists(const BidirectionalPath& path, EdgeId e, std::vector<pair<int, double>> & histogram) const {
        for (size_t j = 0; j < path.Size(); ++j) {
            std::vector<int> distances;
            std::vector<double> weights;
            GetDistances(path.At(j), e, distances, weights);
            if (distances.size() > 0) {
                AddInfoFromEdge(distances, weights, histogram, path.LengthAt(j));
            }
        }
    }

    void FindBestFittedEdgesForClustered(const BidirectionalPath& path, const set<EdgeId>& edges, EdgeContainer& result) const {
        for (EdgeId e : edges) {
            DEBUG("Analyzing edge " << g_.int_id(e))
            std::vector<pair<int, double>> histogram;
            CountAvrgDists(path, e, histogram);
            double sum = 0.0;
            for (size_t j = 0; j < histogram.size(); ++j) {
                sum += histogram[j].second;
            }
            DEBUG("Weight for scaffolding = " << sum << ", threshold = " << cl_weight_threshold_)
            if (math::ls(sum, cl_weight_threshold_)) {
                continue;
            }

            int gap = CountMean(histogram);
            DEBUG("Gap = " << gap)
            if (HasIdealInfo(path, e, gap)) {
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
        const auto& lib = wc_->PairedLibrary();
        //todo lib (and FindJumpEdges) knows its var so it can be counted there
        int is_scatter = int(math::round(lib.GetIsVar() * is_scatter_coeff_));
        for (int i = (int) path.Size() - 1; i >= 0 && path.LengthAt(i) - g_.length(path.At(i)) <= lib.GetISMax(); --i) {
            set<EdgeId> jump_edges_i;
            lib.FindJumpEdges(path.At(i), jump_edges_i,
                               std::max(0, (int)path.LengthAt(i) - is_scatter),
                               //FIXME do we need is_scatter here?
                               int((path.LengthAt(i) + lib.GetISMax() + is_scatter)),
                               0);
            for (EdgeId e : jump_edges_i) {
                if (IsTip(e)) {
                    jumping_edges.insert(e);
                }
            }
        }
        return jumping_edges;
    }

public:


    ScaffoldingExtensionChooser(const Graph& g, shared_ptr<WeightCounter> wc,
                                double cl_weight_threshold,
                                double is_scatter_coeff) :
        ExtensionChooser(g, wc), raw_weight_threshold_(0.0),
        cl_weight_threshold_(cl_weight_threshold),
        is_scatter_coeff_(is_scatter_coeff) {
    }

    EdgeContainer Filter(const BidirectionalPath& path, const EdgeContainer& edges) const override {
        DEBUG("Extension chooser filter, threshold = " << cl_weight_threshold_)
        if (edges.empty()) {
            return edges;
        }
        set<EdgeId> candidates = FindCandidates(path);
        for (auto e : path) {
            if (candidates.find(e) != candidates.end()) {
                DEBUG(g_.int_id(e) << " is removed due to presence in the path")
                candidates.erase(e);
            }
        }
        DEBUG("Found candidates:" << candidates.size())
        EdgeContainer result;
        FindBestFittedEdgesForClustered(path, candidates, result);
        DEBUG("Detected possible edges to scaffold: " << result.size())
        return result;
    }

private:
    DECL_LOGGER("ScaffoldingExtensionChooser");
};

inline bool EdgeWithWeightCompareReverse(const pair<EdgeId, double>& p1,
                                      const pair<EdgeId, double>& p2) {
    return p1.second > p2.second;
}

class LongReadsUniqueEdgeAnalyzer {
    DECL_LOGGER("LongReadsUniqueEdgeAnalyzer")
public:
    LongReadsUniqueEdgeAnalyzer(const Graph& g, const GraphCoverageMap& cov_map,
                                double filter_threshold, double prior_threshold,
                                size_t max_repeat_length, bool uneven_depth)
            : g_(g),
              cov_map_(cov_map),
              filter_threshold_(filter_threshold),
              prior_threshold_(prior_threshold),
              max_repeat_length_(max_repeat_length),
              uneven_depth_(uneven_depth) {

        FindAllUniqueEdges();
    }

    bool IsUnique(EdgeId e) const {
        return unique_edges_.count(e) > 0;
    }

private:
    bool UniqueEdge(EdgeId e) const {
        if (g_.length(e) > max_repeat_length_)
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
                    DEBUG("Checking inconsistency");
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
       VERIFY(!uneven_depth_);
       double sum_cov = 0;
       size_t sum_len = 0;
       size_t total_len = 0;
       for (auto iter = g_.ConstEdgeBegin(); !iter.IsEnd(); ++iter) {
           total_len += g_.length(*iter);
           if (g_.length(*iter) >= max_repeat_length_) {
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
        if (!uneven_depth_)
            FindAllUniqueCoverageEdges();
        DEBUG("Unique edges are found");
    }

    const Graph& g_;
    const GraphCoverageMap& cov_map_;
    double filter_threshold_;
    double prior_threshold_;
    std::set<EdgeId> unique_edges_;
    size_t max_repeat_length_;
    bool uneven_depth_;
};

//class SimpleScaffolding {
//public:
//    SimpleScaffolding(const Graph& g) : g_(g) {}
//
//    BidirectionalPath FindMaxCommonPath(const vector<BidirectionalPath*>& paths,
//                                        size_t max_diff_len) const {
//        BidirectionalPath max_end(g_);
//        for (auto it1 = paths.begin(); it1 != paths.end(); ++it1) {
//            BidirectionalPath* p1 = *it1;
//            for (size_t i = 0; i < p1->Size(); ++i) {
//                if (p1->Length() - p1->LengthAt(i) > max_diff_len) {
//                    break;
//                }
//                bool contain_all = true;
//                for (size_t i1 = i + 1; i1 <= p1->Size() && contain_all; ++i1) {
//                    BidirectionalPath subpath = p1->SubPath(i, i1);
//                    for (auto it2 = paths.begin();  it2 != paths.end() && contain_all; ++it2) {
//                        BidirectionalPath* p2 = *it2;
//                        vector<size_t> positions2 = p2->FindAll(subpath.At(0));
//                        bool contain = false;
//                        for (size_t ipos2 = 0; ipos2 < positions2.size(); ++ipos2) {
//                            size_t pos2 = positions2[ipos2];
//                            if (p2->Length() - p2->LengthAt(pos2) <= max_diff_len
//                                    && EqualEnds(subpath, 0, *p2, pos2, false)) {
//                                contain = true;
//                                break;
//                            }
//                        }
//                        if (!contain) {
//                            contain_all = false;
//                        }
//                    }
//                    if (contain_all && (i1 - i) >= max_end.Size()) {
//                        max_end.Clear();
//                        max_end.PushBack(subpath);
//                    }
//                }
//            }
//        }
//        return max_end;
//    }
//
//private:
//    const Graph& g_;
//};

class LongReadsExtensionChooser : public ExtensionChooser {
public:
    LongReadsExtensionChooser(const Graph& g,
                              const GraphCoverageMap& read_paths_cov_map,
                              double filtering_threshold,
                              double weight_priority_threshold,
                              double unique_edge_priority_threshold,
                              size_t min_significant_overlap,
                              size_t max_repeat_length,
                              bool uneven_depth)
            : ExtensionChooser(g),
              filtering_threshold_(filtering_threshold),
              weight_priority_threshold_(weight_priority_threshold),
              min_significant_overlap_(min_significant_overlap),
              cov_map_(read_paths_cov_map),
              unique_edge_analyzer_(g, cov_map_, filtering_threshold,
                                    unique_edge_priority_threshold,
                                    max_repeat_length, uneven_depth)
    {
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
        path.PrintDEBUG();
        map<EdgeId, double> weights_cands;
        for (auto it = edges.begin(); it != edges.end(); ++it) {
            weights_cands.insert(make_pair(it->e_, 0.0));
        }
        set<EdgeId> filtered_cands;
        auto support_paths = cov_map_.GetCoveringPaths(path.Back());
        DEBUG("Found " << support_paths.size() << " covering paths!!!");
        for (auto it = support_paths.begin(); it != support_paths.end(); ++it) {
            auto positions = (*it)->FindAll(path.Back());
            for (size_t i = 0; i < positions.size(); ++i) {
                if ((int) positions[i] < (int) (*it)->Size() - 1
                        && EqualBegins(path, (int) path.Size() - 1, **it,
                                       positions[i], false)) {
                    DEBUG("Checking unique path_back for " << (*it)->GetId());

                    if (UniqueBackPath(**it, positions[i])) {
                        DEBUG("Success");

                        EdgeId next = (*it)->At(positions[i] + 1);
                        weights_cands[next] += (*it)->GetWeight();
                        filtered_cands.insert(next);
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
            if (unique_edge_analyzer_.IsUnique(path.At(int_pos)) > 0 && g_.length(path.At(int_pos)) >= min_significant_overlap_)
                return true;
            int_pos--;
        }
        return false;
    }

    vector<pair<EdgeId, double> > MapToSortVector(const map<EdgeId, double>& map) const {
        vector<pair<EdgeId, double> > result(map.begin(), map.end());
        std::sort(result.begin(), result.end(), EdgeWithWeightCompareReverse);
        return result;
    }

    double filtering_threshold_;
    double weight_priority_threshold_;
    size_t min_significant_overlap_;
    const GraphCoverageMap& cov_map_;
    LongReadsUniqueEdgeAnalyzer unique_edge_analyzer_;

    DECL_LOGGER("LongReadsExtensionChooser");
};


class CoordinatedCoverageExtensionChooser: public ExtensionChooser {
public:
    CoordinatedCoverageExtensionChooser(const Graph& g, 
            CoverageAwareIdealInfoProvider& coverage_provider, 
            size_t max_edge_length_in_repeat, double delta, size_t min_path_len) :
            ExtensionChooser(g), provider_(coverage_provider), 
            max_edge_length_in_repeat_(max_edge_length_in_repeat), delta_(delta), min_path_len_(min_path_len) {
    }

    EdgeContainer Filter(const BidirectionalPath& path,
            const EdgeContainer& edges) const override {

        if (edges.size() < 2) {
            DEBUG("If unique candidate has not been accepted by previous choosers better not to touch it");
            return EdgeContainer();
        }

        if (path.Length() < min_path_len_) {
            DEBUG("Path is too short");
            return EdgeContainer();
        }

        double path_coverage = provider_.EstimatePathCoverage(path);
        if (math::eq(path_coverage, -1.0) || math::le(path_coverage, 10.0)) {
            DEBUG("Path coverage can't be calculated of too low");
            return EdgeContainer();
        }
        DEBUG("Path coverage is " << path_coverage);

        for (auto e_d : edges) {
            if (path.Contains(g_.EdgeEnd(e_d.e_))) {
                DEBUG("Avoid to create loops");
                return EdgeContainer();
            }
        }
        return FindExtensionTroughRepeat(edges, path_coverage);
    }

private:

    void UpdateCanBeProcessed(VertexId v,
            std::queue<VertexId>& can_be_processed, double path_coverage) const {
        DEBUG("Updating can be processed");
        for (EdgeId e : g_.OutgoingEdges(v)) {
            VertexId neighbour_v = g_.EdgeEnd(e);
            if (g_.length(e) <= max_edge_length_in_repeat_ && CompatibleEdge(e, path_coverage)) {
                DEBUG("Adding vertex " << neighbour_v.int_id()
                                << "through edge " << g_.str(e));
                can_be_processed.push(neighbour_v);
            }
        }
    }

    GraphComponent<Graph> GetRepeatComponent(const VertexId start, double path_coverage) const {
        set<VertexId> vertices_of_component;
        vertices_of_component.insert(start);
        std::queue<VertexId> can_be_processed;
        UpdateCanBeProcessed(start, can_be_processed, path_coverage);
        while (!can_be_processed.empty()) {
            VertexId v = can_be_processed.front();
            can_be_processed.pop();
            if (vertices_of_component.count(v) != 0) {
                DEBUG("Component is too complex");
                return GraphComponent<Graph>::Empty(g_);
            }
            DEBUG("Adding vertex " << g_.str(v) << " to component set");
            vertices_of_component.insert(v);
            UpdateCanBeProcessed(v, can_be_processed, path_coverage);
        }

        return GraphComponent<Graph>::FromVertices(g_, vertices_of_component);
    }

    EdgeContainer FinalFilter(const EdgeContainer& edges,
            EdgeId edge_to_extend) const {
        EdgeContainer result;
        for (auto e_with_d : edges) {
            if (e_with_d.e_ == edge_to_extend) {
                result.push_back(e_with_d);
            }
        }
        return result;
    }

    bool CompatibleEdge(EdgeId e, double path_coverage) const {
        return math::ge(g_.coverage(e), path_coverage * delta_);
    }

    //returns lowest coverage among long compatible edges ahead of e
    //if std::numeric_limits<double>::max() -- no such edges were detected
    //if negative -- abort at once
    double AnalyzeExtension(EdgeId ext, double path_coverage) const {
        double answer = std::numeric_limits<double>::max();

        if (!CompatibleEdge(ext, path_coverage)) {
            DEBUG("Extension coverage is too low");
            return answer;
        }

        if (g_.length(ext) > max_edge_length_in_repeat_) {
            DEBUG("Long extension");
            return g_.coverage(ext);
        } 

        DEBUG("Short extension, launching repeat component analysis");
        GraphComponent<Graph> gc = GetRepeatComponent(g_.EdgeEnd(ext), path_coverage);
        if (gc.v_size() == 0) {
            DEBUG("Component search failed");
            return -1.;
        }

        for (auto e : gc.edges()) {
            if (g_.length(e) > max_edge_length_in_repeat_) {
                DEBUG("Repeat component contains long edges");
                return -1.;
            }
        }

        DEBUG("Checking long sinks");
        for (auto v : gc.exits()) {
            for (auto e : g_.OutgoingEdges(v)) {
                if (g_.length(e) > max_edge_length_in_repeat_ && 
                        CompatibleEdge(e, path_coverage) &&
                        math::ls(g_.coverage(e), answer)) {
                    DEBUG("Updating answer to coverage of edge " << g_.str(e));
                    answer = g_.coverage(e);
                }
            }
        }

        return answer;
    }

    EdgeContainer FindExtensionTroughRepeat(const EdgeContainer& edges, double path_coverage) const {
        static EdgeContainer EMPTY_CONTAINER;

        map<EdgeId, double> good_extension_to_ahead_cov;

        for (auto edge : edges) {
            DEBUG("Processing candidate extension " << g_.str(edge.e_));
            double analysis_res = AnalyzeExtension(edge.e_, path_coverage);

            if (analysis_res == std::numeric_limits<double>::max()) {
                DEBUG("Ignoring extension");
            } else if (math::ls(analysis_res, 0.)) {
                DEBUG("Troubles detected, abort mission");
                return EMPTY_CONTAINER;
            } else {
                good_extension_to_ahead_cov[edge.e_] = analysis_res;
                DEBUG("Extension mapped to ahead coverage of " << analysis_res);
            }
        }

        DEBUG("Number of good extensions is " << good_extension_to_ahead_cov.size());

        if (good_extension_to_ahead_cov.size() == 1) {
            auto extension_info = *good_extension_to_ahead_cov.begin();
            DEBUG("Single extension candidate " << g_.str(extension_info.first));
            if (math::le(extension_info.second, path_coverage / delta_)) {
                DEBUG("Extending");
                return FinalFilter(edges, extension_info.first);
            } else {
                DEBUG("Predicted ahead coverage is too high");
            }
        } else {
            DEBUG("Multiple extension candidates");
        }

        return EMPTY_CONTAINER;
    }
    
    CoverageAwareIdealInfoProvider provider_;
    const size_t max_edge_length_in_repeat_;
    const double delta_;
    const size_t min_path_len_;
    DECL_LOGGER("CoordCoverageExtensionChooser");
};

}
#endif /* EXTENSION_HPP_ */

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
#include "read_cloud_path_extend/read_cloud_dijkstras.hpp"
#include "weight_counter.hpp"
#include "pe_utils.hpp"
#include "assembly_graph/graph_support/scaff_supplementary.hpp"
#include "common/barcode_index/barcode_info_extractor.hpp"

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

class CompositeExtensionChooser: public ExtensionChooser {
    shared_ptr<ExtensionChooser> first_;
    shared_ptr<ExtensionChooser> second_;
public:
    CompositeExtensionChooser(const Graph& g, shared_ptr<ExtensionChooser> first, shared_ptr<ExtensionChooser> second) :
            ExtensionChooser(g), first_(first), second_(second) {}

    EdgeContainer Filter(const BidirectionalPath& path, const EdgeContainer& edges) const override {
        EdgeContainer answer;
        auto first_result = first_->Filter(path, edges);
        if (first_result.size() != 1) {
            auto second_result = second_->Filter(path, edges);
            return second_result;
        }
        return first_result;
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

class ReadCloudExtensionChooser : public ExtensionChooser {
protected:
    typedef shared_ptr<barcode_index::AbstractBarcodeIndexInfoExtractor> barcode_extractor_ptr_t;
    barcode_extractor_ptr_t barcode_extractor_ptr_;
    size_t fragment_len_;
    size_t distance_bound_;
    const ScaffoldingUniqueEdgeStorage& unique_storage_;
    PairedLibConnectionCondition paired_connection_condition_;
    friend class TenXExtensionChecker;

public:
    ReadCloudExtensionChooser(const conj_graph_pack &gp,
                              barcode_extractor_ptr_t extractor,
                              size_t fragment_len,
                              size_t distance_bound,
                              const ScaffoldingUniqueEdgeStorage &unique_storage,
                              PairedLibConnectionCondition paired_connection_condition) :
            ExtensionChooser(gp.g),
            barcode_extractor_ptr_(extractor),
            fragment_len_(fragment_len),
            distance_bound_(distance_bound),
            unique_storage_(unique_storage),
            paired_connection_condition_(paired_connection_condition) {}

    EdgeContainer Filter(const BidirectionalPath &path, const EdgeContainer& candidates) const override {
        EdgeContainer result;
        auto last_unique = FindLastUniqueInPath(path, unique_storage_);
        DEBUG("Searching for next unique edge.")
        result = FindNextUniqueEdge(last_unique.first, candidates);
        DEBUG("Searching finished.")

        DEBUG("next unique edges found, there are " << result.size() << " of them");
        return result;
    }

    EdgeContainer FindNextUniqueEdge(const EdgeId& last_unique, const EdgeContainer& initial_candidates) const {
        VERIFY(unique_storage_.IsUnique(last_unique));
        //find unique edges further in graph

        DEBUG("Searching for best")
        EdgeContainer best_candidates = GetBestCandidates(initial_candidates, last_unique);
        DEBUG("Finished searching")
        DEBUG("Best candidates: " << best_candidates.size());
        for (const auto& candidate: best_candidates) {
            DEBUG(candidate.e_.int_id());
            DEBUG("Conj: " << g_.conjugate(candidate.e_).int_id());
        }
        if (best_candidates.size() == 1) {
            DEBUG("Edge: " << best_candidates[0].e_);
            DEBUG("Distance: " << best_candidates[0].d_);
        }
        return best_candidates;
    }


    std::pair<EdgeId, int> FindLastUniqueInPath(const BidirectionalPath& path,
                                                       const ScaffoldingUniqueEdgeStorage& storage) const {
        for (int i =  (int)path.Size() - 1; i >= 0; --i) {
            if (storage.IsUnique(path.At(i))) {
                return std::make_pair(path.At(i), i);
            }
        }
        return std::make_pair(EdgeId(0), -1);
    }

protected:
    virtual EdgeContainer GetBestCandidates(const EdgeContainer& edges,
                                                        const EdgeId& decisive_edge) const = 0;

    //fixme very ineffective
    EdgeContainer FindClosestEdge(const vector<EdgeWithDistance>& edges) const {
        EdgeContainer closest_edges;
        auto edges_iter = edges.begin();
        size_t path_len_bound = cfg::get().ts_res.topsort_bound;
        do {
            auto edge = *edges_iter;
            VertexId start_vertex = g_.EdgeEnd(edge.e_);
            auto dijkstra = DijkstraHelper<Graph>::CreateBoundedDijkstra(g_, path_len_bound);
            dijkstra.Run(start_vertex);
            bool can_reach_everyone = true;
            for (auto& other_edge: edges) {
                if (other_edge.e_ == edge.e_)
                    continue;
                auto other_start = g_.EdgeStart(other_edge.e_);
                if (!dijkstra.DistanceCounted(other_start)) {
                    can_reach_everyone = false;
                    break;
                }
            }
            if (can_reach_everyone) {
                closest_edges.push_back(*edges_iter);
            }
            ++edges_iter;
        } while (edges_iter != edges.end());
        return closest_edges;
    }

private:
    DECL_LOGGER("ReadCloudExtensionChooser")
};

class TSLRExtensionChooser : public ReadCloudExtensionChooser {
    using ReadCloudExtensionChooser::barcode_extractor_ptr_t;
    using ReadCloudExtensionChooser::barcode_extractor_ptr_;
    using ReadCloudExtensionChooser::unique_storage_;
    using ReadCloudExtensionChooser::fragment_len_;
    double absolute_barcode_threshold_;

public:
    TSLRExtensionChooser(const conj_graph_pack& gp,
                         barcode_extractor_ptr_t extractor,
                         size_t fragment_len,
                         size_t distance_bound,
                         const ScaffoldingUniqueEdgeStorage& unique_storage,
                         const PairedLibConnectionCondition& paired_connection_condition,
                         double absolute_barcode_threshold) :
            ReadCloudExtensionChooser(gp, extractor, fragment_len, distance_bound,
                                      unique_storage, paired_connection_condition),
            absolute_barcode_threshold_(absolute_barcode_threshold) {}

private:
    double GetGapCoefficient(int gap) const {
        VERIFY(gap <= (int)fragment_len_)
        return static_cast<double>(fragment_len_ - gap) /
               static_cast<double>(fragment_len_);
    }

    EdgeContainer GetBestCandidates(const EdgeContainer& edges, const EdgeId& decisive_edge) const {
        //Find edges with barcode score greater than some threshold
        EdgeContainer best_candidates;
        std::copy_if(edges.begin(), edges.end(), std::back_inserter(best_candidates),
                     [this, &decisive_edge](const EdgeWithDistance& edge) {
                         return edge.e_ != decisive_edge and
                                this->barcode_extractor_ptr_->GetIntersectionSizeNormalizedBySecond(decisive_edge, edge.e_) >
                                absolute_barcode_threshold_ * GetGapCoefficient(edge.d_);
                     });
        return best_candidates;
    }
    DECL_LOGGER("TSLRExtensionChooser")
};

struct TenXExtensionChooserStatistics {
    size_t overall_;
    size_t no_candidates_;
    size_t single_candidate_;
    size_t no_barcodes_on_last_edge_;
    size_t too_many_input_candidates_;
    size_t no_candidates_after_initial_filter_;
    size_t initial_filter_helped_;
    size_t too_much_candidates_after_initial_;
    size_t no_candidates_after_middle_filter_;
    size_t middle_filter_helped_;
    size_t multiple_candidates_after_both_;
    size_t pair_of_conjugates_left_;
    size_t conjugate_resolved_;

    TenXExtensionChooserStatistics() :
            overall_(0),
            no_candidates_(0),
            single_candidate_(0),
            no_barcodes_on_last_edge_(0),
            too_many_input_candidates_(0),
            no_candidates_after_initial_filter_(0),
            initial_filter_helped_(0),
            too_much_candidates_after_initial_(0),
            no_candidates_after_middle_filter_(0),
            middle_filter_helped_(0),
            multiple_candidates_after_both_(0),
            pair_of_conjugates_left_(0),
            conjugate_resolved_(0)
    {}

    void PrintStats(const std::string& filename) {
        ofstream fout(filename);
        fout << "Overall: " << overall_ << endl;
        fout << "No candidates: " << no_candidates_ << endl;
        fout << "Single candidate: " << single_candidate_ << endl;
        fout << "No barcodes on last edge: " << no_barcodes_on_last_edge_ << endl;
        fout << "Too many input candidates: " << too_many_input_candidates_ << endl;
        fout << "No candidates after initial filter: " << no_candidates_after_initial_filter_ << endl;
        fout << "Initial filter helped: " << initial_filter_helped_ << endl;
        fout << "Too much after initial: " << too_much_candidates_after_initial_ << endl;
        fout << "No candidates after middle filter: " << no_candidates_after_middle_filter_ << endl;
        fout << "Multiple candidates after middle filter: " << multiple_candidates_after_both_ << endl;
        fout << "Pair of conjugates: " << pair_of_conjugates_left_ << endl;
        fout << "Middle filter helped: " << middle_filter_helped_ << endl;
        fout << "Conjugate resolved: " << conjugate_resolved_ << endl;
    }
};

class TenXExtensionChooser : public ReadCloudExtensionChooser {
    typedef ReadCloudExtensionChooser::barcode_extractor_ptr_t abstract_barcode_extractor_ptr_t;
    typedef barcode_index::FrameBarcodeIndexInfoExtractor frame_extractor_t;
    typedef debruijn_graph::config::debruijn_config::read_cloud_resolver::tenx_resolver tenx_configs_t;
    using ReadCloudExtensionChooser::unique_storage_;
    using ReadCloudExtensionChooser::fragment_len_;
    shared_ptr<frame_extractor_t> barcode_extractor_ptr_;
    tenx_configs_t tenx_configs_;
public:
    static TenXExtensionChooserStatistics stats_;


public:
    TenXExtensionChooser(const conj_graph_pack &gp, abstract_barcode_extractor_ptr_t extractor, size_t fragment_len,
                             size_t distance_bound, const ScaffoldingUniqueEdgeStorage &unique_storage,
                             const PairedLibConnectionCondition& paired_connection_condition,
                             const tenx_configs_t& tenx_configs) :
            ReadCloudExtensionChooser(gp, extractor, fragment_len, distance_bound,
                                      unique_storage, paired_connection_condition),
            barcode_extractor_ptr_(std::static_pointer_cast<frame_extractor_t>(extractor)),
            tenx_configs_(tenx_configs) {}


    static void PrintStats(const string& filename) {
        stats_.PrintStats(filename);
    }

private:
    EdgeContainer GetBestCandidates(const EdgeContainer& edges, const EdgeId& last_unique_edge) const {
        DEBUG("Last unique edge: " << last_unique_edge.int_id());
        DEBUG("Barcodes on last unique edge: " << barcode_extractor_ptr_->GetNumberOfBarcodes(last_unique_edge));
        DEBUG("Input candidates: " << edges.size());
        stats_.overall_++;
        if (edges.size() == 0) {
            stats_.no_candidates_++;
            return edges;
        } else if (edges.size() == 1) {
            stats_.single_candidate_++;
            return edges;
        }
        size_t input_edges_threshold = 2000;
        if (edges.size() > input_edges_threshold) {
            stats_.too_many_input_candidates_++;
            DEBUG("Too many input candidates");
            return edges;
        }
        size_t barcodes = barcode_extractor_ptr_->GetNumberOfBarcodes(last_unique_edge);
        if (barcodes == 0) {
            stats_.no_barcodes_on_last_edge_++;
            return edges;
        }
        //Find edges with barcode score greater than some threshold
        EdgeContainer initial_candidates = InitialSharedBarcodesFilter(edges, last_unique_edge,
                                                                       tenx_configs_.absolute_barcode_threshold,
                                                                       tenx_configs_.initial_abundancy_threshold,
                                                                       tenx_configs_.tail_threshold);
        DEBUG("Initial candidates: ");
        for(const auto& candidate: initial_candidates) {
            DEBUG(candidate.e_.int_id());
            DEBUG("Shared barcodes: " << barcode_extractor_ptr_->
                    CountSharedBarcodesWithFilter(last_unique_edge, candidate.e_,
                                                  tenx_configs_.initial_abundancy_threshold,
                                                  tenx_configs_.tail_threshold));
        }
        if (initial_candidates.size() == 0 ) {
            stats_.no_candidates_after_initial_filter_++;
            return initial_candidates;
        }
        if (initial_candidates.size() == 1) {
            stats_.initial_filter_helped_++;
            return initial_candidates;
        }
        if (initial_candidates.size() > tenx_configs_.max_initial_candidates) {
            stats_.too_much_candidates_after_initial_++;
            DEBUG("Too many initial candidates");
            return initial_candidates;
        }
        DEBUG("After initial check: " << initial_candidates.size());
        EdgeContainer next_candidates = ReadCloudClosestFilter(initial_candidates, last_unique_edge,
                                                               tenx_configs_.internal_gap_threshold,
                                                               tenx_configs_.middle_abundancy_threshold);
        DEBUG("After middle check: " << next_candidates.size());
        DEBUG("Middle candidates: ");
        for(const auto& candidate: next_candidates) {
            DEBUG(candidate.e_.int_id());
        }
        if (next_candidates.size() == 0) {
            stats_.no_candidates_after_middle_filter_++;
            DEBUG("No candidates after middle");
            return next_candidates;
        }
        if (next_candidates.size() == 1) {
            stats_.middle_filter_helped_++;
            DEBUG("Middle helped");
            DEBUG("Canditate: " << next_candidates.back().e_.int_id());
            return next_candidates;
        }
        EdgeContainer result = next_candidates;
        stats_.multiple_candidates_after_both_++;
        DEBUG("Multiple after middle");
        if (next_candidates.size() == 2) {
            EdgeId edge = next_candidates[0].e_;
            EdgeId other = next_candidates[1].e_;
            if (edge == g_.conjugate(other)) {
                DEBUG("Pair of conjugates left");
                stats_.pair_of_conjugates_left_++;
                //fixme magic constants
                result = ConjugateFilter(last_unique_edge, next_candidates[0], next_candidates[1],
                                         1000, 2000, 0.2);
                if (result.size() == 1) {
                    stats_.conjugate_resolved_++;
                }
            }
        }
        return result;
    }

    EdgeContainer InitialSharedBarcodesFilter(const EdgeContainer &candidates, const EdgeId &decisive_edge,
                                              size_t shared_threshold, size_t abundancy_threshold,
                                              size_t tail_threshold) const {
        EdgeContainer result;
        std::copy_if(candidates.begin(), candidates.end(), std::back_inserter(result),
                     [this, &decisive_edge, shared_threshold,
                             abundancy_threshold, tail_threshold](const EdgeWithDistance& edge) {
                         return edge.e_ != decisive_edge and
                                 this->barcode_extractor_ptr_->AreEnoughSharedBarcodesWithFilter(decisive_edge, edge.e_,
                                                                                                 shared_threshold,
                                                                                                 abundancy_threshold,
                                                                                                 tail_threshold);
                     });
        return result;
    }

    EdgeContainer TopologyClosestFilter(const EdgeContainer& candidates, const size_t distance_bound) const {
        EdgeContainer result;
        std::copy_if(candidates.begin(), candidates.end(), std::back_inserter(result),
                     [this, &candidates, distance_bound](const EdgeWithDistance& edge) {
                         return TopologyCheck(edge.e_, candidates, distance_bound);
                     });
        return result;
    }

    bool TopologyCheck(const EdgeId& candidate, const EdgeContainer& candidates, const size_t distance_bound) const {
        std::unordered_set<EdgeId> reached_unique_edges = GetReachedUniqueEdges(candidate, distance_bound);
        bool are_all_candidates_reached = AreCandidatesReached(candidate, candidates, reached_unique_edges);
        return are_all_candidates_reached;
    }

    unordered_set <EdgeId> GetReachedUniqueEdges(const EdgeId& edge, const size_t distance_bound) const {
        DijkstraHelper<debruijn_graph::Graph>::BoundedDijkstra dijkstra(
                DijkstraHelper<debruijn_graph::Graph>::CreateBoundedDijkstra(g_, distance_bound));
        dijkstra.Run(g_.EdgeEnd(edge));
        unordered_set <EdgeId> reached_unique_edges;
        for (auto v: dijkstra.ReachedVertices()) {
            for (auto connected: g_.OutgoingEdges(v)) {
                if (unique_storage_.IsUnique(connected) and
                    dijkstra.GetDistance(v) < distance_bound) {
                    reached_unique_edges.insert(connected);
                }
            }
        }
        return reached_unique_edges;
    }

    bool AreCandidatesReached(const EdgeId& candidate, const EdgeContainer& candidates, const unordered_set<EdgeId>& reached) const {
        for (const auto& ewd: candidates) {
            if (reached.find(ewd.e_) == reached.end() and reached.find(g_.conjugate(ewd.e_)) == reached.end()
                and ewd.e_ != candidate and ewd.e_ != g_.conjugate(candidate)) {
                return false;
            }
        }
        return true;
    }

    EdgeContainer ReadCloudClosestFilter(const EdgeContainer &candidates, const EdgeId &decisive_edge,
                                         size_t len_threshold, size_t abundancy_threshold) const {
        EdgeContainer result;
        std::copy_if(candidates.begin(), candidates.end(), std::back_inserter(result),
                     [this, &decisive_edge, &candidates, len_threshold, abundancy_threshold](const EdgeWithDistance& edge) {
                         return ReadCloudClosestCheck(decisive_edge, edge.e_, candidates, len_threshold,
                                                      abundancy_threshold);
                     });
        return result;
    }

    bool ReadCloudClosestCheck(const EdgeId &decisive_edge, const EdgeId &candidate,
                               const EdgeContainer &other_candidates, size_t len_threshold, size_t abundancy_threshold) const {
        DEBUG("Middle check for: " << candidate.int_id());
        for (const auto& other : other_candidates) {
            EdgeId other_edge = other.e_;
            if (other_edge != candidate and !IsBetween(candidate, decisive_edge, other_edge,
                                                       len_threshold, abundancy_threshold)) {
                return false;
                break;
            }
        }
        return true;
    }

    EdgeContainer ConjugateFilter(const EdgeId& decisive_edge, const EdgeWithDistance& edgewd,
                                  const EdgeWithDistance& conjugatewd, size_t gap_threshold_left,
                                  size_t gap_threshold_right, double fraction_threshold) const  {
        const EdgeId edge = edgewd.e_;
        const EdgeId conjugate = conjugatewd.e_;
        VERIFY(g_.conjugate(edge) == conjugate);
        EdgeContainer result;
        size_t edge_voters = 0;
        size_t conj_voters = 0;
        auto common_barcodes = barcode_extractor_ptr_->GetSharedBarcodes(decisive_edge, edge);
        for (const auto& barcode: common_barcodes) {
            size_t gap = barcode_extractor_ptr_->GetMinPos(edge, barcode);
            size_t conj_gap = barcode_extractor_ptr_->GetMinPos(conjugate, barcode);
            if (gap < gap_threshold_left and conj_gap > gap_threshold_right) {
                edge_voters++;
            }
            if (gap > gap_threshold_right and conj_gap < gap_threshold_left) {
                conj_voters++;
            }
        }
        size_t common_size = common_barcodes.size();
        double edge_fraction = static_cast<double>(edge_voters) / static_cast<double>(common_size);
        double conj_fraction = static_cast<double>(conj_voters) / static_cast<double>(common_size);
        if (edge_fraction - conj_fraction > fraction_threshold) {
            result.push_back(edgewd);
        }
        if (conj_fraction - edge_fraction > fraction_threshold) {
            result.push_back(conjugatewd);
        }
        return result;
    }

    bool IsBetween(const EdgeId& middle, const EdgeId& left, const EdgeId& right,
                   size_t len_threshold, size_t abundance_threshold) const {
        DEBUG("Checking against: " << right.int_id())
        auto side_barcodes = barcode_extractor_ptr_->GetSharedBarcodes(left, right);
        size_t middle_length = g_.length(middle);
        size_t current_length = 0;
        DEBUG("Side barcodes before filter: " << side_barcodes.size());
        size_t side_barcodes_after_filter = 0;
        for (const auto& barcode: side_barcodes) {
            //todo optimize after custom filter implementation
            size_t left_count = barcode_extractor_ptr_->GetNumberOfReads(left, barcode);
            size_t right_count = barcode_extractor_ptr_->GetNumberOfReads(right, barcode);
            if (!(barcode_extractor_ptr_->HasBarcode(middle, barcode))
                and left_count >= abundance_threshold
                and right_count >= abundance_threshold) {
                ++side_barcodes_after_filter;
                size_t right_length = barcode_extractor_ptr_->GetMinPos(right, barcode);
                size_t left_length = g_.length(left) - barcode_extractor_ptr_->GetMaxPos(left, barcode);
                current_length += (left_length + right_length + middle_length);
            }
        }
        DEBUG("Side barcodes after filter: " << side_barcodes_after_filter);
        size_t sum_length_threshold = len_threshold * side_barcodes_after_filter;
        DEBUG("Current length: " << current_length);
        return current_length <= sum_length_threshold;
    }

    DECL_LOGGER("10XExtensionChooser")
};

class ReadCloudGapExtensionChooser : public ExtensionChooser {
protected:
    typedef shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr_t;
    barcode_extractor_ptr_t barcode_extractor_ptr_;
    EdgeId target_edge_;
    const ScaffoldingUniqueEdgeStorage& unique_storage_;
    struct ReadCloudGapExtensionChooserConfigs {
        const size_t shared_threshold_;
        const size_t read_threshold;
        const size_t gap_threshold_;
        const size_t scan_bound_;

        ReadCloudGapExtensionChooserConfigs(const size_t shared_threshold_, const size_t count_threshold_,
                                            const size_t gap_threshold_, const size_t scan_bound_) :
                                                                           shared_threshold_(shared_threshold_),
                                                                           read_threshold(count_threshold_),
                                                                           gap_threshold_(gap_threshold_),
                                                                           scan_bound_(scan_bound_) {}
    };
    const ReadCloudGapExtensionChooserConfigs configs_;

public:
    ReadCloudGapExtensionChooser(const Graph &g, const barcode_extractor_ptr_t &barcode_extractor_ptr_,
                                 const EdgeId &target_edge_, const ScaffoldingUniqueEdgeStorage& unique_storage_)
            : ExtensionChooser(g), barcode_extractor_ptr_(barcode_extractor_ptr_),
              target_edge_(target_edge_), unique_storage_(unique_storage_), configs_(0, 1, 5000, 500) {}

    virtual EdgeContainer Filter(const BidirectionalPath& path, const EdgeContainer& edges) const override {
        DEBUG("Filter started");
        EdgeContainer result;
        const EdgeId last_unique = FindLastUniqueInPath(path);
        DEBUG("Trying to close gap");
        DEBUG("Last edge: " << path.Back().int_id());
        DEBUG("Last unique: " << last_unique.int_id());
        DEBUG("Target edge: " << target_edge_.int_id());
        DEBUG(edges.size() << " candidates");
        if (last_unique.int_id() == 0 or target_edge_.int_id() == 0) {
            return result;
        }

        auto is_edge_supported_by_clouds = [this, &last_unique](const EdgeWithDistance& edge) {
            return IsEdgeSupportedByClouds(edge.e_, last_unique, this->target_edge_, this->configs_.shared_threshold_,
                                           this->configs_.read_threshold, this->configs_.gap_threshold_);
        };
        std::copy_if(edges.begin(), edges.end(), std::back_inserter(result), is_edge_supported_by_clouds);
        if (result.size() > 0) {
            DEBUG("Found supported edges");
            return result;
        }

        auto are_supported_by_clouds_reachable = [this, &last_unique](const EdgeWithDistance& edge) {
            return AreSupportedByCloudsReachable(edge.e_, last_unique, this->target_edge_,
                                                 this->configs_.shared_threshold_, this->configs_.read_threshold,
                                                 this->configs_.gap_threshold_, this->configs_.scan_bound_);
        };
        std::copy_if(edges.begin(), edges.end(), std::back_inserter(result), are_supported_by_clouds_reachable);
        if (result.size() > 0) {
            DEBUG("Supported edges are reachable from some candidates");
        } else {
            DEBUG("Candidates were not found.");
        }
        DEBUG("Filter finished")
        return result;
    }
private:
    bool IsEdgeSupportedByClouds(const EdgeId &edge, const EdgeId &last_edge, const EdgeId &target_edge,
                                 const size_t shared_threshold, const size_t count_threshold,
                                 const size_t gap_threshold) const {
        return barcode_extractor_ptr_->AreEnoughSharedBarcodesWithFilter(last_edge, edge, shared_threshold,
                                                                         count_threshold, gap_threshold) and
               barcode_extractor_ptr_->AreEnoughSharedBarcodesWithFilter(edge, target_edge, shared_threshold,
                                                                         count_threshold, gap_threshold);
    }

    bool AreSupportedByCloudsReachable(const EdgeId &edge, const EdgeId &last_edge, const EdgeId &target_edge,
                                       const size_t shared_threshold, const size_t count_threshold,
                                       const size_t gap_threshold, const size_t scan_bound) const {
        DijkstraHelper<debruijn_graph::Graph>::BoundedDijkstra dijkstra(
                DijkstraHelper<debruijn_graph::Graph>::CreateBoundedDijkstra(g_, scan_bound));
        dijkstra.Run(g_.EdgeEnd(edge));
        DEBUG("Running dij from candidate " << edge.int_id());
        bool found_supported = false;
        //todo use custom dijkstra that stops on supported edges to improve performance
        for (auto v: dijkstra.ReachedVertices()) {
            if (dijkstra.GetDistance(v) < scan_bound) {
                for (auto connected: g_.OutgoingEdges(v)) {
                    DEBUG("Checking edge " << connected.int_id());
                    if (IsEdgeSupportedByClouds(connected, last_edge, target_edge,
                                                shared_threshold, count_threshold, gap_threshold)) {
                        found_supported = true;
                        DEBUG("supported");
                        break;
                    }
                }
            }
            if (found_supported) break;
        }
        return found_supported;
    }

    EdgeId FindLastUniqueInPath(const BidirectionalPath& path) const {
        for (int i =  static_cast<int>(path.Size()) - 1; i >= 0; --i) {
            if (unique_storage_.IsUnique(path.At(i))) {
                return path.At(i);
            }
        }
        return EdgeId(0);
    }

    DECL_LOGGER("ReadCloudGapExtensionChooser");
};
}
#endif /* EXTENSION_HPP_ */

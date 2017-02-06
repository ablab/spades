//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef EXTENSION_HPP_
#define EXTENSION_HPP_

#include "pe_utils.hpp"
#include "weight_counter.hpp"

#include "alignment/rna/ss_coverage.hpp"
#include "assembly_graph/components/graph_component.hpp"
#include "assembly_graph/graph_support/scaff_supplementary.hpp"
#include "barcode_index/barcode_mapper.hpp"
#include "bounded_dijkstra.hpp"
#include "modules/alignment/rna/ss_coverage.hpp"

#include <cfloat>
#include <fstream>
#include <iostream>
#include <map>
#include <queue>

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
        if (path.Empty())
            return 0;

        int lastEdge = (int) path.Size() - 1;
        do {
            lastEdge = ExcludeTrivial(path, edges, lastEdge);
            bool bulge = true;

            if (lastEdge >= 0) {
                VertexId v = g_.EdgeEnd(path[lastEdge]);
                VertexId u = g_.EdgeStart(path[lastEdge]);

                for (EdgeId e : g_.IncomingEdges(v)) {
                    if (g_.EdgeStart(e) != u) {
                        bulge = false;
                        break;
                    }
                }

                if (!bulge)
                    break;

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
    PreserveSimplePathsAnalyzer(const Graph &g)
            : PathAnalyzer(g) { }

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
        if (path.Empty())
            return 0;

        int lastEdge = (int) path.Size() - 1;
        bool has_bulge = false;
        do {
            lastEdge = PathAnalyzer::ExcludeTrivial(path, edges, lastEdge);

            if (lastEdge >= 0) {
                VertexId v = g_.EdgeEnd(path[lastEdge]);
                VertexId u = g_.EdgeStart(path[lastEdge]);
                has_bulge = true;

                for (EdgeId e : g_.IncomingEdges(v)) {
                    if (g_.EdgeStart(e) != u) {
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
    std::shared_ptr<WeightCounter> wc_;
    //FIXME memory leak?!
    std::vector<ExtensionChooserListener *> listeners_;

    double weight_threshold_;

public:
    ExtensionChooser(const Graph &g, std::shared_ptr<WeightCounter> wc = nullptr, double weight_threshold = -1.):
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

    bool WeightCounterBased() const noexcept {
        return wc_ != nullptr;
    }

    std::shared_ptr<WeightCounter> wc() const {
        return wc_;
    }

protected:
    double IdealInfo(EdgeId e1, EdgeId e2, size_t dist, bool additive = false) const {
        return wc_->PairedLibrary().IdealPairedInfo(e1, e2, (int) dist, additive);
    }

    bool HasIdealInfo(EdgeId e1, EdgeId e2, size_t dist) const {
        return math::gr(IdealInfo(e1, e2, dist), 0.);
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
    std::shared_ptr<ExtensionChooser> first_;
    std::shared_ptr<ExtensionChooser> second_;

public:
    JointExtensionChooser(const Graph& g,
                          std::shared_ptr<ExtensionChooser> first,
                          std::shared_ptr<ExtensionChooser> second): ExtensionChooser(g),
        first_(first), second_(second) {
    }

    EdgeContainer Filter(const BidirectionalPath& path, const EdgeContainer& edges) const override {
        EdgeContainer answer;
        auto r1 = first_->Filter(path, edges);
        auto r2 = second_->Filter(path, edges);
        for (const auto& ewd1 : r1) {
            for (const auto& ewd2 : r2) {
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
    //> 1
    double coverage_margin_;
    //> 1
    double max_coverage_variation_;

    double min_upper_coverage_;

public:
    SimpleCoverageExtensionChooser(const SSCoverageStorage& coverage_storage, const Graph& g,
                                   double coverage_margin, double max_coverage_variation, double min_upper_coverage) :
        ExtensionChooser(g), coverage_storage_(coverage_storage),
        coverage_margin_(coverage_margin),
        max_coverage_variation_(max_coverage_variation),
        min_upper_coverage_(min_upper_coverage) {
        VERIFY(math::gr(coverage_margin_, 1.0));
        VERIFY(math::gr(max_coverage_variation_, 1.0));
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

        bool reverse_cov =  math::ls(coverage_storage_.GetCoverage(path_edge_at_split), coverage_storage_.GetCoverage(path_edge_at_split, true));
        return Filter(path, edges, index, reverse_cov);
    }

private:
    EdgeContainer Filter(const BidirectionalPath& path, const EdgeContainer& edges, size_t split_index,  bool reverse) const {
        DEBUG("Coverage extension chooser");
        VERIFY(split_index > 0);
        VERIFY(edges.size() == 2);

        DEBUG("Split found at " << split_index);
        EdgeId path_edge_at_split = path[split_index - 1];
        EdgeId other_edge_at_split = GetOtherEdgeAtSplit(g_.EdgeEnd(path_edge_at_split), path_edge_at_split);

        EdgeId candidate1 = edges.front().e_;
        EdgeId candidate2 = edges.back().e_;
        VERIFY(other_edge_at_split != EdgeId());

        auto cov_path = coverage_storage_.GetCoverage(path_edge_at_split, reverse);
        auto cov_other = coverage_storage_.GetCoverage(other_edge_at_split, reverse);
        auto cov_e1 = coverage_storage_.GetCoverage(candidate1, reverse);
        auto cov_e2 = coverage_storage_.GetCoverage(candidate2, reverse);

        if (IsCoverageSimilar(cov_path, cov_other, coverage_margin_) || IsCoverageSimilar(cov_e1, cov_e2, coverage_margin_)) {
            DEBUG("Margin is too low: path = " << cov_path << ", other = " << cov_other << ", ex1 = " << cov_e1 << ", ex2 = " << cov_e2);
            return EdgeContainer();
        }

        double high_path = std::max(cov_path, cov_other);
        double low_path = std::min(cov_path, cov_other);
        double high_ex = std::max(cov_e1, cov_e2);
        double low_ex = std::min(cov_e1, cov_e2);

        if (!IsEnoughCoverage(low_path, high_path) || !IsEnoughCoverage(low_ex, high_ex)) {
            DEBUG("Coverage is too low: path = " << high_path << ", other = " << high_ex);
            return EdgeContainer();
        }

        EdgeContainer result;
        if (math::gr(cov_path, cov_other)) {
            DEBUG("Path coverage is high, edge " << g_.int_id(path_edge_at_split) << ", path cov = " << cov_path << ", other " << cov_other);
            if (IsCoverageSimilar(high_path, high_ex, max_coverage_variation_))
                result.emplace_back(math::gr(cov_e1, cov_e2) ? candidate1 : candidate2, 0);
        } else {
            DEBUG("Path coverage is low, edge " << g_.int_id(path_edge_at_split) << ", path cov = " << cov_path << ", other " << cov_other);
            if (IsCoverageSimilar(low_path, low_ex, max_coverage_variation_))
                result.emplace_back(math::ls(cov_e1, cov_e2) ? candidate1 : candidate2, 0);
        }

        VERIFY(result.size() <= 1);
        return result;
    }

    bool IsEnoughCoverage(double low, double high) const {
        return math::eq(low, 0.0) || math::ge(high, min_upper_coverage_);
    }

    bool IsCoverageSimilar(double cov1, double cov2, double threshold) const {
        if (math::eq(cov2, 0.0) || math::eq(cov1, 0.0)) {
            return false;
        }

        double diff = math::gr(cov1, cov2) ? cov1 / cov2 : cov2 / cov1;
        return math::le(diff, threshold);
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
            weights.emplace(weight, *iter);
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
    ExcludingExtensionChooser(const Graph &g, std::shared_ptr<WeightCounter> wc, PathAnalyzer analyzer,
                              double weight_threshold, double priority) :
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
        std::map<size_t, unsigned> edge_2_extension_cnt;
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

    SimpleExtensionChooser(const Graph &g, std::shared_ptr<WeightCounter> wc, double weight_threshold, double priority) :
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
                               std::shared_ptr<WeightCounter> wc,
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

    RNAExtensionChooser(const Graph &g, std::shared_ptr<WeightCounter> wc, double weight_threshold, double priority) :
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
    LongEdgeExtensionChooser(const Graph& g, std::shared_ptr<WeightCounter> wc, double weight_threshold, double priority) :
        ExcludingExtensionChooser(g, wc, PathAnalyzer(g), weight_threshold, priority) {
    }
};

class ScaffoldingExtensionChooser : public ExtensionChooser {
    typedef ExtensionChooser base;
    typedef std::vector<std::pair<int, double>> Histogram;
    double raw_weight_threshold_;
    double cl_weight_threshold_;
    const double is_scatter_coeff_ = 3.0;
    double relative_cutoff_;

    void AddInfoFromEdge(const std::vector<int> &distances, const std::vector<double> &weights,
                         Histogram &histogram, size_t len_to_path_end) const {
        for (size_t l = 0; l < distances.size(); ++l) {
            //todo commented out condition seems unnecessary and should be library dependent! do we need "max(0" there?
            if (/*distances[l] > max(0, (int) len_to_path_end - int(1000)) && */math::ge(weights[l], raw_weight_threshold_)) {
                histogram.emplace_back(distances[l] - (int) len_to_path_end, weights[l]);
            }
        }
    }

    int CountMean(const Histogram &histogram) const {
        double dist = 0.0;
        double sum = 0.0;
        for (size_t i = 0; i < histogram.size(); ++i) {
            dist += histogram[i].first * histogram[i].second;
            sum += histogram[i].second;
        }
        dist /= sum;
        return (int) round(dist);
    }

    void GetDistances(EdgeId e1, EdgeId e2, std::vector<int> &dist, std::vector<double> &w) const {
        wc_->PairedLibrary().CountDistances(e1, e2, dist, w);
    }

    void CountAvrgDists(const BidirectionalPath &path, EdgeId e, Histogram &histogram) const {
        for (size_t j = 0; j < path.Size(); ++j) {
            std::vector<int> distances;
            std::vector<double> weights;
            GetDistances(path.At(j), e, distances, weights);
            if (distances.size() > 0) {
                AddInfoFromEdge(distances, weights, histogram, path.LengthAt(j));
            }
        }
    }

    void FindBestFittedEdgesForClustered(const BidirectionalPath &path,
                                         const std::set<EdgeId>& edges,
                                         EdgeContainer& result) const {
        for (EdgeId e : edges) {
            Histogram histogram;
            DEBUG("Analyzing edge " << g_.int_id(e))
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
            if (math::gr(relative_cutoff_, 0.0)) {
                //double avg_cov = (g_.coverage(path.Back()) + g_.coverage(e)) / 2.0;
                double iw = CountIdealInfo(path, e, gap);
                if (math::ls(sum / iw, relative_cutoff_)) {
                    continue;
                }
            }

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

    std::set<EdgeId> FindCandidates(const BidirectionalPath &path) const {
        std::set<EdgeId> jumping_edges;
        const auto& lib = wc_->PairedLibrary();
        //todo lib (and FindJumpEdges) knows its var so it can be counted there
        int is_scatter = int(math::round(lib.GetIsVar() * is_scatter_coeff_));
        for (int i = (int) path.Size() - 1; i >= 0 && path.LengthAt(i) - g_.length(path.At(i)) <= lib.GetISMax(); --i) {
            std::set<EdgeId> jump_edges_i;
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
    ScaffoldingExtensionChooser(const Graph& g, std::shared_ptr<WeightCounter> wc,
                                double cl_weight_threshold,
                                double is_scatter_coeff,
                                double relative_cutoff = 0.0) :
        ExtensionChooser(g, wc), raw_weight_threshold_(0.0),
        cl_weight_threshold_(cl_weight_threshold),
        is_scatter_coeff_(is_scatter_coeff), relative_cutoff_(relative_cutoff) {
    }

    EdgeContainer Filter(const BidirectionalPath& path, const EdgeContainer& edges) const override {
        DEBUG("Extension chooser filter, threshold = " << cl_weight_threshold_)
        if (edges.empty()) {
            return edges;
        }
        auto candidates = FindCandidates(path);
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

public:
    double CountIdealInfo(const BidirectionalPath& path, EdgeId edge, size_t gap) const {
        double sum = 0.0;
        for (int i = (int) path.Size() - 1; i >= 0; --i) {
            double w = IdealInfo(path[i], edge, (int) path.LengthAt(i) + gap);
            sum += w;
        }
        double avg_cov = (g_.coverage(path.Back()) + g_.coverage(edge)) / 2.0;
        double correction_coeff = avg_cov / ((double(wc_->PairedLibrary().GetRL()) - double(g_.k())) * 2);
        return sum * correction_coeff;
    }

private:
    DECL_LOGGER("ScaffoldingExtensionChooser");
};

template <class T>
inline bool EdgeWithWeightCompareReverse(const std::pair<T, double>& p1,
                                         const std::pair<T, double>& p2) {
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

    bool CheckInconsistence(const BidirectionalPath& path1, size_t pos1,
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

    std::pair<double, double> GetSubPathsWeights(const BidirectionalPath& cand1,
                                                 const BidirectionalPath& cand2,
                                                 const BidirectionalPathSet& cov_paths) const {
        double weight1 = 0.0;
        double weight2 = 0.0;
        for (const BidirectionalPath *path : cov_paths) {
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
       for (EdgeId e : g_.edges()) {
           total_len += g_.length(e);
           if (g_.length(e) >= max_repeat_length_) {
               sum_cov += g_.coverage(e) * (double) g_.length(e);
               sum_len += g_.length(e);
           }
       }
       if (sum_len * 4 < total_len) return;
       sum_cov /= (double)sum_len;
       DEBUG("average coverage of long edges: " << sum_cov) ;
       for (EdgeId e : g_.edges()) {
           if (g_.length(e) > 500 && (double)g_.coverage(e) < 1.2 * sum_cov) {
               if (unique_edges_.count(e))
                   continue;

               unique_edges_.insert(e);
               unique_edges_.insert(g_.conjugate(e));
               DEBUG("Added coverage based unique edge " << g_.int_id(e) << " len "<< g_.length(e) << " " << g_.coverage(e));
           }
       }
   }


    void FindAllUniqueEdges() {
        DEBUG("Looking for unique edges");
        for (EdgeId e : g_.edges()) {
            if (!UniqueEdge(e))
                continue;

            unique_edges_.insert(e);
            unique_edges_.insert(g_.conjugate(e));
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

class LongReadsRNAExtensionChooser : public ExtensionChooser {
public:
    LongReadsRNAExtensionChooser(const Graph& g,
                                 const GraphCoverageMap& read_paths_cov_map,
                                 double filtering_threshold,
                                 size_t min_significant_overlap)
        : ExtensionChooser(g),
          filtering_threshold_(filtering_threshold),
          min_significant_overlap_(min_significant_overlap),
          cov_map_(read_paths_cov_map)
    {
        DEBUG("Created LongReadsRNAExtensionChooser with params: filtering_threshold_ = " << filtering_threshold_ << ", min_significant_overlap_ = " << min_significant_overlap_);
    }

    /* Choose extension as correct only if we have reads that traverse a unique edge from the path and this extension.
     * Edge is unique if all reads mapped to this edge are consistent.
     * Two reads are consistent if they can form one path in the graph.
     */
    EdgeContainer Filter(const BidirectionalPath& path,
                         const EdgeContainer& edges) const override {
        if (edges.empty()) {
            return edges;
        }
        DEBUG("We are in Filter of LongReadsRNAExtensionChooser");
        path.PrintDEBUG();

        std::map<EdgeId, double> weights_candidates;
        DEBUG("Candidates")
        for (auto it = edges.begin(); it != edges.end(); ++it) {
            DEBUG(g_.int_id(it->e_));
            weights_candidates.emplace(it->e_, 0.0);
        }

        std::set<EdgeId> filtered_candidates;
        auto support_paths = cov_map_.GetCoveringPaths(path.Back());
        DEBUG("Found " << support_paths.size() << " supporting paths");
        for (auto it = support_paths.begin(); it != support_paths.end(); ++it) {
            auto supporting_path = *it;
            auto positions = supporting_path->FindAll(path.Back());

            for (size_t i = 0; i < positions.size(); ++i) {
                if ((int) positions[i] < (int) supporting_path->Size() - 1
                    && EqualBegins(path, path.Size() - 1, *supporting_path, positions[i], false)) {
                    DEBUG("Supporting path matches, Checking unique path_back for " << supporting_path->GetId());

                    if (UniqueBackPath(*supporting_path, positions[i])) {
                        DEBUG("Success");

                        EdgeId next = supporting_path->At(positions[i] + 1);
                        weights_candidates[next] += supporting_path->GetWeight();
                        filtered_candidates.insert(next);
                    }
                }
            }
        }
        DEBUG("Supported candidates");
        for (auto iter = weights_candidates.begin(); iter != weights_candidates.end(); ++iter) {
            DEBUG("Candidate " << g_.int_id(iter->first) << " weight " << iter->second);
        }
        std::vector<std::pair<EdgeId, double> > sort_res = MapToSortVector(weights_candidates);
        if (sort_res.size() < 1 || sort_res[0].second < filtering_threshold_) {
            filtered_candidates.clear();
        }
        EdgeContainer result;
        for (auto it = edges.begin(); it != edges.end(); ++it) {
            if (filtered_candidates.find(it->e_) != filtered_candidates.end()) {
                result.push_back(*it);
                DEBUG("Inserting to results: " << g_.int_id(it->e_));
            }
        }
        if (result.size() == 0) {
            DEBUG("Long reads don't help =(");
        }
        return result;
    }

private:

    bool UniqueBackPath(const BidirectionalPath& path, size_t pos) const {
        int int_pos = (int) pos;
        while (int_pos >= 0) {
            if (g_.length(path.At(int_pos)) >= min_significant_overlap_)
                return true;
            int_pos--;
        }
        return false;
    }

    std::vector<std::pair<EdgeId, double> > MapToSortVector(const std::map<EdgeId, double>& map) const {
        std::vector<std::pair<EdgeId, double> > result(map.begin(), map.end());
        std::sort(result.begin(), result.end(), EdgeWithWeightCompareReverse<EdgeId>);
        return result;
    }

    double filtering_threshold_;
    size_t min_significant_overlap_;
    const GraphCoverageMap& cov_map_;

    DECL_LOGGER("LongReadsRNAExtensionChooser");
};

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
    EdgeContainer Filter(const BidirectionalPath &path, const EdgeContainer &edges) const override {
        if (edges.empty()) {
            return edges;
        }
        DEBUG("We are in Filter of LongReadsExtensionChooser");
        path.PrintDEBUG();
        std::map<EdgeId, double> weights_cands;
        for (const auto &edge : edges) {
            weights_cands.emplace(edge.e_, 0.0);
        }
        std::set<EdgeId> filtered_cands;
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
        auto sort_res = MapToSortVector(weights_cands);
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

    std::vector<std::pair<EdgeId, double>> MapToSortVector(const std::map<EdgeId, double>& map) const {
        std::vector<std::pair<EdgeId, double>> result(map.begin(), map.end());
        std::sort(result.begin(), result.end(), EdgeWithWeightCompareReverse<EdgeId>);
        return result;
    }

    double filtering_threshold_;
    double weight_priority_threshold_;
    size_t min_significant_overlap_;
    const GraphCoverageMap& cov_map_;
    LongReadsUniqueEdgeAnalyzer unique_edge_analyzer_;

    DECL_LOGGER("LongReadsExtensionChooser");
};

class TrustedContigsExtensionChooser : public ExtensionChooser {
public:
    TrustedContigsExtensionChooser(const Graph& g,
                                    const GraphCoverageMap& read_paths_cov_map,
                                    double filtering_threshold,
                                    double weight_priority_threshold,
                                    double unique_edge_priority_threshold,
                                    size_t min_significant_overlap,
                                    size_t max_repeat_length,
                                    bool uneven_depth,
                                    bool use_low_quality_matching = false)
            : ExtensionChooser(g)
            , filtering_threshold_(filtering_threshold)
            , weight_priority_threshold_(weight_priority_threshold)
            , min_significant_overlap_(min_significant_overlap)
            , cov_map_(read_paths_cov_map)
            , unique_edge_analyzer_(g, cov_map_, filtering_threshold,
                                    unique_edge_priority_threshold,
                                    max_repeat_length, uneven_depth)
            , use_low_quality_matching_(use_low_quality_matching)
    {}

    /// @returns the possible next edge of the path
    EdgeContainer Filter(const BidirectionalPath &path, const EdgeContainer &) const override {
        DEBUG("We are in Filter of TrustedContigsExtensionChooser");
        path.PrintDEBUG();

        std::map<EdgeWithDistance, double> weights_cands;
        auto filtered_cands = GetHighQualityCandidats(path, weights_cands);
        if (use_low_quality_matching_ && filtered_cands.empty())
            filtered_cands = GetLowQualityCandidats(path, weights_cands);

        DEBUG("Candidates:");
        for (auto const & iter : weights_cands)
            DEBUG("Candidate " << g_.int_id(iter.first.e_) << " weight " << iter.second);

        auto sort_res = MapToSortVector(weights_cands);
        DEBUG("sort res " << sort_res.size() << " tr " << weight_priority_threshold_);
        if (sort_res.empty() || sort_res[0].second < filtering_threshold_) {
            filtered_cands.clear();
        } else if (sort_res.size() > 1) {
            if (sort_res[0].second > weight_priority_threshold_ * sort_res[1].second) {
                filtered_cands.clear();
                filtered_cands.insert(sort_res[0].first);
            } else {
                for (size_t i = 0; i < sort_res.size(); ++i) {
                    if (sort_res[i].second * weight_priority_threshold_ < sort_res[0].second)
                        filtered_cands.erase(sort_res[i].first);
                }
            }
        }
        EdgeContainer result(filtered_cands.begin(), filtered_cands.end());
        if (result.size() != 1)
            DEBUG("Trusted contigs don't help =(");
        return result;
    }

private:

    /// @return (there are covering paths that contains 'path[start_pos]')
    std::pair<bool, std::set<EdgeWithDistance>> GetCandidates(const BidirectionalPath &path,
                                             std::map<EdgeWithDistance, double> &weights_cands,
                                             size_t start_pos,
                                             const std::function<std::pair<bool, size_t>(const BidirectionalPath&, size_t)> &comparator) const
    {
        VERIFY(start_pos < path.Size());
        std::set<EdgeWithDistance> filtered_cands;
        auto support_paths = cov_map_.GetCoveringPaths(path[start_pos]);
        DEBUG("Found " << support_paths.size() << " covering paths!!!");
        for (auto const & it : support_paths) {
            for (auto pos : it->FindAll(path[start_pos])) {
                if (pos + 1 < it->Size()) {
                    auto tmp = comparator(*it, pos);
                    auto& is_good_path = tmp.first;
                    auto& matched_len = tmp.second;
                    if (is_good_path) {
                        auto gap = it->GapAt(pos + 1);
                        EdgeWithDistance next = {it->At(pos + 1), gap.gap, std::move(gap.gap_seq)};
                        weights_cands[next] += static_cast<double>(matched_len)*it->GetWeight();
                        filtered_cands.insert(next);
                    }
                }
            }
        }
        return {!support_paths.empty(), std::move(filtered_cands)};
    }

    std::set<EdgeWithDistance> GetHighQualityCandidats(const BidirectionalPath &path, std::map<EdgeWithDistance, double> &weights_cands) const {
        auto start_pos = path.Size() - 1;
        auto comparator = [&path, start_pos, th = this] (const BidirectionalPath &coverage_path, size_t pos) {
            auto is_good_path = (FirstNotEqualPosition(path, start_pos, coverage_path, pos, false) == -1ul);
            auto matched_len = start_pos + 1;
            auto privilege_scalar = (th->HasUniqueEdge(path, 0, matched_len) ? 3 : 1);
            return std::pair<bool, size_t>(is_good_path, matched_len * privilege_scalar);
        };
        return GetCandidates(path, weights_cands, start_pos, comparator).second;
    }

    std::set<EdgeWithDistance> GetLowQualityCandidats(const BidirectionalPath &path, std::map<EdgeWithDistance, double> &weights_cands) const {
        DEBUG("Fallback mode");
        auto get_comparator = [&path, th = this](size_t start_pos){
            /// compares the two path prefixes from 'start_pos', might skip several nonunique edges
            return [&path, start_pos, th] (const BidirectionalPath &coverage_path, size_t pos) {
                auto pos1 = static_cast<int>(start_pos);
                auto pos2 = static_cast<int>(pos);
                auto skipped_non_unique_edges = 0;
                while (pos1 >= 0 && pos2 >= 0) {
                    if (path[pos1] == coverage_path[pos2]) {
                        skipped_non_unique_edges = 0;
                        --pos2;
                    } else {
                        ++skipped_non_unique_edges;
                        if (th->IsUniqueEdge(path[pos1]) || skipped_non_unique_edges > 0)
                            break;
                    }
                    --pos1;
                }
                auto matched_len = pos-pos2;
                auto is_good_path = th->HasUniqueEdge(coverage_path, pos2+1, matched_len);
                return std::pair<bool, size_t>(is_good_path, matched_len);
            };
        };
        bool used_unique_edge = false;
        for (int i = 0; i < std::min((int)path.Size(), 2) && !used_unique_edge; ++i) {
            DEBUG("iteration: " << i);
            auto pos = path.Size() - 1 - i;
            used_unique_edge = IsUniqueEdge(path[pos]);
            auto tmp = GetCandidates(path, weights_cands, pos, get_comparator(pos));
            auto& coverage_paths_are_found = tmp.first;
            auto& filtered_cands = tmp.second;
            if (!filtered_cands.empty() || coverage_paths_are_found)
                return std::move(filtered_cands);
        }
        DEBUG("Fallback mode does not help");
        return {};
    }

    bool IsUniqueEdge(EdgeId edge) const {
        return unique_edge_analyzer_.IsUnique(edge) && g_.length(edge) >= min_significant_overlap_;
    }

    bool HasUniqueEdge(const BidirectionalPath& path, size_t from, size_t len) const {
        for (size_t i = 0; i < len; ++i) {
            auto edge = path.At(from + i);
            if (IsUniqueEdge(edge))
                return true;
        }
        return false;
    }

    std::vector<std::pair<EdgeWithDistance, double>> MapToSortVector(const std::map<EdgeWithDistance, double>& map) const {
        std::vector<std::pair<EdgeWithDistance, double>> result(map.begin(), map.end());
        std::sort(result.begin(), result.end(), EdgeWithWeightCompareReverse<EdgeWithDistance>);
        return result;
    }

    double filtering_threshold_;
    double weight_priority_threshold_;
    size_t min_significant_overlap_;
    const GraphCoverageMap& cov_map_;
    LongReadsUniqueEdgeAnalyzer unique_edge_analyzer_;
    bool use_low_quality_matching_;

    DECL_LOGGER("TrustedContigsExtensionChooser");
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

        for (const auto& e_d : edges) {
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

    omnigraph::GraphComponent<Graph> GetRepeatComponent(const VertexId start, double path_coverage) const {
        std::set<VertexId> vertices_of_component;
        vertices_of_component.insert(start);
        std::queue<VertexId> can_be_processed;
        UpdateCanBeProcessed(start, can_be_processed, path_coverage);
        while (!can_be_processed.empty()) {
            VertexId v = can_be_processed.front();
            can_be_processed.pop();
            if (vertices_of_component.count(v) != 0) {
                DEBUG("Component is too complex");
                return omnigraph::GraphComponent<Graph>::Empty(g_);
            }
            DEBUG("Adding vertex " << g_.str(v) << " to component set");
            vertices_of_component.insert(v);
            UpdateCanBeProcessed(v, can_be_processed, path_coverage);
        }

        return omnigraph::GraphComponent<Graph>::FromVertices(g_, vertices_of_component);
    }

    EdgeContainer FinalFilter(const EdgeContainer& edges,
            EdgeId edge_to_extend) const {
        EdgeContainer result;
        for (const auto& e_with_d : edges) {
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
        auto gc = GetRepeatComponent(g_.EdgeEnd(ext), path_coverage);
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

        std::map<EdgeId, double> good_extension_to_ahead_cov;

        for (const auto& edge : edges) {
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
    ScaffoldingUniqueEdgeStorage unique_storage_;
    friend class TenXExtensionChecker;

public:
    ReadCloudExtensionChooser(const conj_graph_pack& gp,
                              barcode_extractor_ptr_t extractor,
                              size_t fragment_len,
                              size_t distance_bound,
                              const ScaffoldingUniqueEdgeStorage& unique_storage) :
            ExtensionChooser(gp.g),
            barcode_extractor_ptr_(extractor),
            fragment_len_(fragment_len),
            distance_bound_(distance_bound),
            unique_storage_(unique_storage) {
    }

    //todo remove code duplication with ExtensionChooser2015
    EdgeContainer Filter(const BidirectionalPath &path, const EdgeContainer&) const override {
        EdgeContainer result;

        pair<EdgeId, int> last_unique = FindLastUniqueInPath(path, unique_storage_);

        if (last_unique.second == -1) {
            return result;
        }
        EdgeContainer candidates = GetInitialCandidates(last_unique.first, path);

        DEBUG("At edge " << path.Back().int_id());
        DEBUG("Decisive edge " << last_unique.first.int_id());
        DEBUG("Decisive edge barcodes: " << barcode_extractor_ptr_->GetTailBarcodeNumber(last_unique.first));

        DEBUG("Searching for next unique edge.")
        result = FindNextUniqueEdge(last_unique.first, candidates);
        DEBUG("Searching finished.")

        DEBUG("next unique edges found, there are " << result.size() << " of them");
//Backward check. We connected edges iff they are best continuation to each other.
//        if (result.size() == 1) {
//            //We should reduce gap size with length of the edges that came after last unique.
//            result[0].d_ -= int (path.LengthAt(last_unique.second) - g_.length(last_unique.first));
//            DEBUG("For edge " << g_.int_id(last_unique.first) << " unique next edge "<< result[0].e_ <<" found, doing backwards check ");
//            EdgeContainer backwards_check = FindNextUniqueEdge(g_.conjugate(result[0].e_));
//            if ((backwards_check.size() != 1) || (g_.conjugate(backwards_check[0].e_) != last_unique.first)) {
//                DEBUG("Backward check failed")
//                result.clear();
//            }
//        }
        return result;
    }

    EdgeContainer FindNextUniqueEdge(const EdgeId& last_unique, const EdgeContainer& initial_candidates) const {
        VERIFY(unique_storage_.IsUnique(last_unique));
        //find unique edges further in graph

        DEBUG("Searching for best")
        EdgeContainer best_candidates = GetBestCandidates(initial_candidates, last_unique);
        DEBUG("Finished searching")
        DEBUG("Best candidates: " << best_candidates.size());

        //Try to find topologically closest edge to resolve loops
//        if (best_candidates.size() > 1) {
//            FilterByTopSort(best_candidates, result);
//        }

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
    virtual EdgeContainer GetInitialCandidates(EdgeId last_edge, const BidirectionalPath& path) const {
        ExtensionChooser::EdgeContainer candidates;
        vector<EdgeId> initial_candidates;
        auto put_checker = BarcodePutChecker<Graph>(g_, last_edge, unique_storage_, initial_candidates);
        auto dij = BarcodeDijkstra<Graph>::CreateBarcodeBoundedDijkstra(g_, distance_bound_, put_checker);
        dij.Run(g_.EdgeEnd(last_edge));
        DEBUG("Initial candidates: " << initial_candidates.size());
        DEBUG("Dijkstra finished");
        candidates.reserve(initial_candidates.size());
        for (const auto& edge : initial_candidates) {
            VERIFY(unique_storage_.IsUnique(last_edge));
            if (edge != last_edge and edge != g_.conjugate(last_edge) and
                path.FindFirst(edge) == -1 and path.FindFirst(g_.conjugate(edge)) == -1) {
                candidates.push_back(EdgeWithDistance(edge, dij.GetDistance(g_.EdgeStart(edge))));
            }
        }
        return candidates;
    }

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

    void FilterByTopSort(const vector<EdgeWithDistance> &best_candidates, EdgeContainer &result) const {
        auto closest_edges =  FindClosestEdge(best_candidates);
        if (closest_edges.size() != 1) {
            DEBUG("Unable to find single topmin edge.");
            for (auto& edge : best_candidates) {
                result.push_back(edge);
            }
        }
        else {
            DEBUG("Found topmin edge");
            result.push_back(closest_edges.back());
        }
    }

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
                         double absolute_barcode_threshold) :
            ReadCloudExtensionChooser(gp, extractor, fragment_len, distance_bound, unique_storage),
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
    using ReadCloudExtensionChooser::unique_storage_;
    using ReadCloudExtensionChooser::fragment_len_;
    shared_ptr<frame_extractor_t> barcode_extractor_ptr_;
    size_t absolute_barcode_threshold_;
    size_t tail_threshold_;
    size_t max_initial_candidates_;
    size_t internal_gap_threshold_;
    size_t initial_abundancy_threshold_;
    size_t middle_abundancy_threshold_;
    friend class TenXExtensionChecker;
public:
    static TenXExtensionChooserStatistics stats_;


public:
    TenXExtensionChooser(const conj_graph_pack &gp, abstract_barcode_extractor_ptr_t extractor, size_t fragment_len,
                             size_t distance_bound, const ScaffoldingUniqueEdgeStorage &unique_storage,
                             size_t absolute_barcode_threshold, size_t tail_threshold, size_t max_initial_candidates,
                             size_t internal_gap_threshold, size_t initial_abundancy_threshold_,
                             size_t middle_abundancy_threshold_) :
            ReadCloudExtensionChooser(gp, extractor, fragment_len, distance_bound, unique_storage),
            barcode_extractor_ptr_(std::static_pointer_cast<frame_extractor_t>(extractor)),
            absolute_barcode_threshold_(absolute_barcode_threshold),
            tail_threshold_(tail_threshold),
            max_initial_candidates_(max_initial_candidates),
            internal_gap_threshold_(internal_gap_threshold),
            initial_abundancy_threshold_(initial_abundancy_threshold_),
            middle_abundancy_threshold_(middle_abundancy_threshold_) {}


    static void PrintStats(const string& filename) {
        stats_.PrintStats(filename);
    }

private:
    EdgeContainer GetBestCandidates(const EdgeContainer& edges, const EdgeId& decisive_edge) const {
        DEBUG(decisive_edge.int_id());
        stats_.overall_++;
        if (edges.size() == 0) {
            stats_.no_candidates_++;
            return edges;
        } else if (edges.size() == 1) {
            stats_.single_candidate_++;
            return edges;
        }
        size_t barcodes = barcode_extractor_ptr_->GetTailBarcodeNumber(decisive_edge);
        if (barcodes == 0) {
            stats_.no_barcodes_on_last_edge_++;
            return edges;
        }
        //Find edges with barcode score greater than some threshold
        EdgeContainer initial_candidates = InitialFilter(edges, decisive_edge, absolute_barcode_threshold_,
                                                         initial_abundancy_threshold_, tail_threshold_);
        DEBUG("Initial candidates: ");
        for(const auto& candidate: initial_candidates) {
            DEBUG(candidate.e_.int_id());
        }
        if (initial_candidates.size() == 0 ) {
            stats_.no_candidates_after_initial_filter_++;
            return initial_candidates;
        }
        if (initial_candidates.size() == 1) {
            stats_.initial_filter_helped_++;
            return initial_candidates;
        }
        if (initial_candidates.size() > max_initial_candidates_) {
            stats_.too_much_candidates_after_initial_++;
            return initial_candidates;
        }
        DEBUG("After initial check: " << initial_candidates.size());
        EdgeContainer next_candidates = MiddleFilter(initial_candidates, decisive_edge,
                                                     internal_gap_threshold_, middle_abundancy_threshold_);
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
            return next_candidates;
        }

        EdgeContainer result = next_candidates;

        stats_.multiple_candidates_after_both_++;
        DEBUG("Multiple after middle");
        if (next_candidates.size() == 2) {
            EdgeId edge = next_candidates[0].e_;
            EdgeId other = next_candidates[1].e_;
            if (edge == g_.conjugate(other)) {
                stats_.pair_of_conjugates_left_++;
                //fixme magic constants
                result = ConjugateFilter(decisive_edge, next_candidates[0], next_candidates[1],
                                         1000, 2000, 0.2);
                if (result.size() == 1) {
                    stats_.conjugate_resolved_++;
                }
            }
        }
        return result;
    }

    EdgeContainer InitialFilter(const EdgeContainer& candidates, const EdgeId& decisive_edge,
                                size_t shared_threshold, size_t abundancy_threshold,
                                size_t tail_threshold) const {
        EdgeContainer result;
        std::copy_if(candidates.begin(), candidates.end(), std::back_inserter(result),
                     [this, &decisive_edge, shared_threshold,
                             abundancy_threshold, tail_threshold](const EdgeWithDistance& edge) {
                         return edge.e_ != decisive_edge and
                                 this->barcode_extractor_ptr_->AreEnoughSharedBarcodes(decisive_edge, edge.e_,
                                                                                       shared_threshold,
                                                                                       abundancy_threshold,
                                                                                       tail_threshold);
                     });
        return result;
    }

    EdgeContainer MiddleFilter(const EdgeContainer& candidates, const EdgeId& decisive_edge,
                               size_t len_threshold, size_t abundancy_threshold) const {
        EdgeContainer result;
        std::copy_if(candidates.begin(), candidates.end(), std::back_inserter(result),
                     [this, &decisive_edge, &candidates, len_threshold, abundancy_threshold](const EdgeWithDistance& edge) {
                         return MiddleCheck(decisive_edge, edge.e_, candidates, len_threshold, abundancy_threshold);
                     });
        return result;
    }

    bool MiddleCheck(const EdgeId& decisive_edge, const EdgeId& candidate,
                     const EdgeContainer& other_candidates, size_t len_threshold, size_t abundancy_threshold) const {
        bool result = true;
        for (const auto& other : other_candidates) {
            EdgeId other_edge = other.e_;
            if (other_edge != candidate and !IsBetween(candidate, decisive_edge, other_edge,
                                                       len_threshold, abundancy_threshold)) {
                result = false;
                break;
            }
        }
        return result;
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
        auto common_barcodes = barcode_extractor_ptr_->GetIntersection(decisive_edge, edge);
        for (const auto barcode: common_barcodes) {
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
                   size_t len_threshold, size_t abundancy_threshold) const {
        auto side_barcodes = barcode_extractor_ptr_->GetIntersection(left, right);
        size_t middle_length = g_.length(middle);
        size_t sum_length_threshold = len_threshold * side_barcodes.size();
        size_t current_length = 0;
        DEBUG("Side barcodes: " << side_barcodes.size());
        VERIFY(side_barcodes.size() > absolute_barcode_threshold_);
        for (const auto& barcode: side_barcodes) {
            //todo optimize after custom filter implementation
            size_t left_count = barcode_extractor_ptr_->GetInfo(left, barcode).GetCount();
            size_t right_count = barcode_extractor_ptr_->GetInfo(right, barcode).GetCount();
            if (!(barcode_extractor_ptr_->has_barcode(middle, barcode))
                and left_count >= abundancy_threshold
                and right_count >= abundancy_threshold) {
                size_t right_length = barcode_extractor_ptr_->GetMinPos(right, barcode);
                size_t left_length = g_.length(left) - barcode_extractor_ptr_->GetMaxPos(left, barcode);
                current_length += (left_length + right_length + middle_length);
            }
        }
        DEBUG("Current length: " << current_length);
        return current_length <= sum_length_threshold;
    }

    DECL_LOGGER("10XExtensionChooser")
};
}
#endif /* EXTENSION_HPP_ */

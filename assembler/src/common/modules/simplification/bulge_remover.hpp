//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * bulge_remover.hpp
 *
 *  Created on: Apr 13, 2011
 *      Author: sergey
 */

#pragma once

#include "assembly_graph/graph_support/parallel_processing.hpp"
#include "assembly_graph/graph_support/basic_edge_conditions.hpp"
#include "assembly_graph/graph_support/graph_processing_algorithm.hpp"
#include "assembly_graph/paths/path_processor.hpp"
#include "assembly_graph/graph_support/comparators.hpp"
#include "assembly_graph/components/graph_component.hpp"
#include "sequence/sequence_tools.hpp"
#include "utils/standard_base.hpp"
#include <cmath>
#include <stack>
#include "math/xmath.h"

namespace omnigraph {

template<class Graph>
struct SimplePathCondition {
    typedef typename Graph::EdgeId EdgeId;
    const Graph& g_;

    SimplePathCondition(const Graph& g) :
            g_(g) {

    }

    bool operator()(EdgeId edge, const vector<EdgeId>& path) const {
        if (edge == g_.conjugate(edge))
            return false;
        for (size_t i = 0; i < path.size(); ++i)
            if (edge == path[i] || edge == g_.conjugate(path[i]))
                return false;
        for (size_t i = 0; i < path.size(); ++i) {
            if (path[i] == g_.conjugate(path[i])) {
                return false;
            }
            for (size_t j = i + 1; j < path.size(); ++j)
                if (path[i] == path[j] || path[i] == g_.conjugate(path[j]))
                    return false;
        }
        return true;
    }
};

template<class Graph>
bool TrivialCondition(typename Graph::EdgeId,
        const vector<typename Graph::EdgeId>& path) {
    for (size_t i = 0; i < path.size(); ++i)
        for (size_t j = i + 1; j < path.size(); ++j)
            if (path[i] == path[j])
                return false;
    return true;
}

template<class Graph>
class MostCoveredSimpleAlternativePathChooser: public PathProcessor<Graph>::Callback {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    const Graph& g_;
    EdgeId forbidden_edge_;

    double max_coverage_;
    vector<EdgeId> most_covered_path_;

public:

    MostCoveredSimpleAlternativePathChooser(const Graph& g, EdgeId edge) :
            g_(g), forbidden_edge_(edge), max_coverage_(-1.0) {

    }

    void HandleReversedPath(const vector<EdgeId>& reversed_path) override {
        vector<EdgeId> path = this->ReversePath(reversed_path);
        double path_cov = AvgCoverage(g_, path);
        for (size_t i = 0; i < path.size(); i++) {
            if (path[i] == forbidden_edge_)
                return;
        }
        if (path_cov > max_coverage_ && SimplePathCondition<Graph>(g_)(forbidden_edge_, path)) {
            max_coverage_ = path_cov;
            most_covered_path_ = path;
        }
    }

    double max_coverage() {
        return max_coverage_;
    }

    const vector<EdgeId>& most_covered_path() {
        return most_covered_path_;
    }
};

inline size_t CountMaxDifference(size_t absolute_diff, size_t length, double relative_diff) {
    return std::max((size_t) std::floor(relative_diff * (double) length), absolute_diff);
}

template<class Graph>
class BulgeGluer {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef std::function<void(EdgeId edge, const vector<EdgeId>& path)> BulgeCallbackF;
    Graph& g_;
    BulgeCallbackF opt_callback_;
    std::function<void(EdgeId)> removal_handler_;

    void InnerProcessBulge(EdgeId edge, const vector<EdgeId>& path) {

        EnsureEndsPositionAligner aligner(CumulativeLength(g_, path),
                g_.length(edge));
        double prefix_length = 0.;
        vector<size_t> bulge_prefix_lengths;

        for (EdgeId e : path) {
            prefix_length += (double) g_.length(e);
            bulge_prefix_lengths.push_back(aligner.GetPosition((size_t) prefix_length));
        }

        EdgeId edge_to_split = edge;
        size_t prev_length = 0;

        TRACE("Process bulge " << path.size() << " edges");

        VERIFY(bulge_prefix_lengths.back() == g_.length(edge));

        for (size_t i = 0; i < path.size(); ++i) {
            if (bulge_prefix_lengths[i] > prev_length) {
                if (bulge_prefix_lengths[i] - prev_length
                        != g_.length(edge_to_split)) {

                    TRACE("SplitEdge " << g_.str(edge_to_split));
                    TRACE(
                            "Start: " << g_.str(g_.EdgeStart(edge_to_split)));
                    TRACE(
                            "Start: " << g_.str(g_.EdgeEnd(edge_to_split)));

                    pair<EdgeId, EdgeId> split_result = g_.SplitEdge(
                            edge_to_split,
                            bulge_prefix_lengths[i] - prev_length);

                    edge_to_split = split_result.second;

                    TRACE("GlueEdges " << g_.str(split_result.first));
                    g_.GlueEdges(split_result.first, path[i]);

                } else {
                    TRACE("GlueEdges " << g_.str(edge_to_split));
                    g_.GlueEdges(edge_to_split, path[i]);
                }
            }
            prev_length = bulge_prefix_lengths[i];
        }
    }

public:

    BulgeGluer(Graph& g, BulgeCallbackF opt_callback = 0,
               std::function<void(EdgeId)> removal_handler = 0) :
               g_(g),
               opt_callback_(opt_callback),
               removal_handler_(removal_handler) {

    }

    void operator()(EdgeId edge, const vector<EdgeId>& path) {
        if (opt_callback_)
            opt_callback_(edge, path);

        if (removal_handler_)
            removal_handler_(edge);

        VertexId start = g_.EdgeStart(edge);
        VertexId end = g_.EdgeEnd(edge);

        TRACE("Projecting edge " << g_.str(edge));
        InnerProcessBulge(edge, path);

        TRACE("Compressing start vertex " << g_.str(start));
        g_.CompressVertex(start);

        TRACE("Compressing end vertex " << g_.str(end));
        g_.CompressVertex(end);
    }

};

template<class Graph>
class AlternativesAnalyzer {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    const Graph& g_;
    double max_coverage_;
    size_t max_length_;
    double max_relative_coverage_;
    size_t max_delta_;
    double max_relative_delta_;
    size_t max_edge_cnt_;

    static vector<EdgeId> EmptyPath() {
        static vector<EdgeId> vec = {};
        return vec;
    }

    /**
     * Checks if alternative path is simple (doesn't contain conjugate edges, edge e or conjugate(e))
     * and its average coverage * max_relative_coverage_ is greater than g.coverage(e)
     */
    bool BulgeCondition(EdgeId e, const vector<EdgeId>& path,
            double path_coverage) const {
        return math::ge(path_coverage * max_relative_coverage_,
                g_.coverage(e)) && SimplePathCondition<Graph>(g_)(e, path);
    }

public:
    AlternativesAnalyzer(const Graph& g, double max_coverage, size_t max_length,
                         double max_relative_coverage, size_t max_delta,
                         double max_relative_delta, size_t max_edge_cnt) :
                         g_(g),
                         max_coverage_(max_coverage),
                         max_length_(max_length),
                         max_relative_coverage_(max_relative_coverage),
                         max_delta_(max_delta),
                         max_relative_delta_(max_relative_delta),
                         max_edge_cnt_(max_edge_cnt) {
        DEBUG("Created alternatives analyzer max_length=" << max_length
        << " max_coverage=" << max_coverage
        << " max_relative_coverage=" << max_relative_coverage
        << " max_delta=" << max_delta
        << " max_relative_delta=" << max_relative_delta);
    }

    vector<EdgeId> operator()(EdgeId e) const {
        if (g_.length(e) > max_length_ || math::gr(g_.coverage(e), max_coverage_)) {
            return EmptyPath();
        }

        size_t kplus_one_mer_coverage = (size_t) math::round((double) g_.length(e) * g_.coverage(e));
        TRACE("Processing edge " << g_.str(e) << " and coverage " << kplus_one_mer_coverage);

        size_t delta = CountMaxDifference(max_delta_, g_.length(e), max_relative_delta_);

        MostCoveredSimpleAlternativePathChooser<Graph> path_chooser(g_, e);

        VertexId start = g_.EdgeStart(e);
        TRACE("Start " << g_.str(start));
        VertexId end = g_.EdgeEnd(e);
        TRACE("End " << g_.str(end));

        ProcessPaths(g_, (g_.length(e) > delta) ? g_.length(e) - delta : 0,
                g_.length(e) + delta, start, end, path_chooser, max_edge_cnt_);

        const vector<EdgeId>& path = path_chooser.most_covered_path();
        if (!path.empty()) {
            VERIFY(g_.EdgeStart(path[0]) == start);
            VERIFY(g_.EdgeEnd(path.back()) == end);
        }

        double path_coverage = path_chooser.max_coverage();
        if (math::gr(path_coverage, 0.)) {
            TRACE("Best path with coverage " << path_coverage << " is " << PrintPath(g_, path));

            if (BulgeCondition(e, path, path_coverage)) {
                TRACE("Satisfied condition");
                return path;
            } else {
                TRACE("Didn't satisfy condition");
                return EmptyPath();
            }
        } else {
            TRACE("Didn't find alternative");
            return EmptyPath();
        }
    }

    double max_coverage() const {
        return max_coverage_;
    }

    size_t max_length() const {
        return max_length_;
    }

private:
    DECL_LOGGER("AlternativesAnalyzer");
};

template<class Graph>
func::TypedPredicate<typename Graph::EdgeId>
NecessaryBulgeCondition(const Graph& g, size_t max_length, double max_coverage) {
    return AddAlternativesPresenceCondition(g,
                                            func::And(LengthUpperBound<Graph>(g, max_length),
                                                     CoverageUpperBound<Graph>(g, max_coverage)));
}

/**
 * This class removes simple bulges from given graph with the following algorithm: it iterates through all edges of
 * the graph and for each edge checks if this edge is likely to be a simple bulge
 * if edge is judged to be one it is removed.
 */
template<class Graph>
class BulgeRemover: public PersistentProcessingAlgorithm<Graph,
                                                        typename Graph::EdgeId,
                                                        CoverageComparator<Graph>> {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef PersistentProcessingAlgorithm<Graph, EdgeId, CoverageComparator<Graph>> base;

protected:

    /*virtual*/
    bool Process(EdgeId e) {
        TRACE("Considering edge " << this->g().str(e)
                      << " of length " << this->g().length(e)
                      << " and avg coverage " << this->g().coverage(e));

        if (!HasAlternatives(this->g(), e)) {
            TRACE("Not possible bulge edge");
            return false;
        }

        vector<EdgeId> alternative = alternatives_analyzer_(e);
        if (!alternative.empty()) {
            gluer_(e, alternative);
            return true;
        }
        return false;
    }

public:

    typedef std::function<void(EdgeId edge, const vector<EdgeId>& path)> BulgeCallbackF;

    BulgeRemover(Graph& g, const std::shared_ptr<InterestingElementFinder<Graph, EdgeId>>& interesting_finder,
            const AlternativesAnalyzer<Graph>& alternatives_analyzer,
            BulgeCallbackF opt_callback = 0,
            std::function<void(EdgeId)> removal_handler = 0,
            bool track_changes = true) :
            base(g,
                 interesting_finder,
                 /*canonical_only*/true,
                 CoverageComparator<Graph>(g),
                 track_changes),
            alternatives_analyzer_(alternatives_analyzer),
            gluer_(g, opt_callback, removal_handler) {
    }

private:
    AlternativesAnalyzer<Graph> alternatives_analyzer_;
    BulgeGluer<Graph> gluer_;
private:
    DECL_LOGGER("BulgeRemover")
};

template<class Graph>
class ParallelBulgeRemover : public PersistentAlgorithmBase<Graph> {
private:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef std::shared_ptr<InterestingElementFinder<Graph, EdgeId>> CandidateFinderPtr;
    typedef SmartSetIterator<Graph, EdgeId, CoverageComparator<Graph>> SmartEdgeSet;

    size_t buff_size_;
    double buff_cov_diff_;
    double buff_cov_rel_diff_;
    AlternativesAnalyzer<Graph> alternatives_analyzer_;
    BulgeGluer<Graph> gluer_;
    CandidateFinderPtr interesting_edge_finder_;
    //todo remove
    bool tracking_;

    size_t curr_iteration_;

    SmartEdgeSet it_;

    static vector<EdgeId> EmptyPath() {
        static vector<EdgeId> vec = {};
        return vec;
    }

    struct BulgeInfo : private boost::noncopyable {
        size_t id;
        EdgeId e;
        std::vector<EdgeId> alternative;

        BulgeInfo() :
            id(-1ul) {
        }

        //passing by value is not a mistake!
        BulgeInfo(size_t id_, EdgeId e_, std::vector<EdgeId> alternative_) :
            id(id_), e(e_), alternative(std::move(alternative_)) {

        }

        BulgeInfo(BulgeInfo&& that) {
            *this = std::move(that);
        }

        BulgeInfo& operator= (BulgeInfo&& that) {
            id = that.id;
            e = that.e;
            alternative = std::move(that.alternative);
            return *this;
        }

//        BulgeInfo(size_t id_, EdgeId e_, std::vector<EdgeId>&& alternative_) :
//            id(id_), e(e_), alternative(std::move(alternative_)) {
//
//        }
//
        bool operator< (const BulgeInfo& that) const {
//            VERIFY_MSG(id != that.id, "Ooops " << id);
            return id < that.id;
        }

        std::string str(const Graph& g) const {
            std::stringstream ss;
            ss << "BulgeInfo " << id
                    << " e: " << g.str(e)
                    << " path: " << PrintPath(g, alternative);
            return ss.str();
        }

    };

    bool CheckInteracting(const BulgeInfo& info, const std::unordered_set<EdgeId>& involved_edges) const {
        if (involved_edges.count(info.e))
            return true;
        for (EdgeId e : info.alternative)
            if (involved_edges.count(e))
                return true;
        return false;
    }

    void AccountEdge(EdgeId e, std::unordered_set<EdgeId>& involved_edges) const {
        TRACE("Pushing edge " << this->g().str(e));
        involved_edges.insert(e);
        EdgeId conj = this->g().conjugate(e);
        TRACE("Pushing edge " << this->g().str(conj));
        involved_edges.insert(conj);
    }

    void AccountEdges(const BulgeInfo& info, std::unordered_set<EdgeId>& involved_edges) const {
        AccountEdge(info.e, involved_edges);
        for (EdgeId e : info.alternative) {
            AccountEdge(e, involved_edges);
        }
    }

    //false if time to stop
    bool FillEdgeBuffer(vector<EdgeId>& buffer, func::TypedPredicate<EdgeId> proceed_condition) {
        VERIFY(buffer.empty());
        DEBUG("Filling edge buffer of size " << buff_size_);
        perf_counter perf;
        double low_cov = 0.;
        double cov_diff = 0.;
        while (!it_.IsEnd() && buffer.size() < buff_size_) {
            EdgeId e = *it_;
            TRACE("Current edge " << this->g().str(e));
            if (!proceed_condition(e)) {
                TRACE("Stop condition was reached.");
                //need to release last element of the iterator to make it replaceable by new elements
                it_.ReleaseCurrent();
                return false;
            }

            double cov = this->g().coverage(e);
            if (buffer.empty()) {
                low_cov = cov;
                cov_diff = max(buff_cov_diff_, buff_cov_rel_diff_ * low_cov);
            } else {
                if (math::gr(cov, low_cov + cov_diff)) {
                    //need to release last element of the iterator to make it replaceable by new elements
                    it_.ReleaseCurrent();
                    return true;
                }
            }
            TRACE("Potential bulge edge");
            buffer.push_back(e);
            ++it_;
        }

        DEBUG("Filled in " << perf.time() << " seconds");
        if (buffer.size() == buff_size_) {
            TRACE("Buffer filled");
            return true;
        } else {
            TRACE("No more edges in iterator");
            return false;
        }
    }

    std::vector<std::vector<BulgeInfo>> FindBulges(const std::vector<EdgeId>& edge_buffer) const {
        DEBUG("Looking for bulges (in parallel). Edge buffer size " << edge_buffer.size());
        perf_counter perf;
        std::vector<std::vector<BulgeInfo>> bulge_buffers(omp_get_max_threads());
        size_t n = edge_buffer.size();
        //order is in agreement with coverage
        #pragma omp parallel for schedule(guided)
        for (size_t i = 0; i < n; ++i) {
            EdgeId e = edge_buffer[i];
            auto alternative = alternatives_analyzer_(e);
            if (!alternative.empty()) {
                bulge_buffers[omp_get_thread_num()].push_back(BulgeInfo(i, e, std::move(alternative)));
            }
        }
        DEBUG("Bulges found in " << perf.time() << " seconds");
        return bulge_buffers;
    }

    std::vector<BulgeInfo> MergeBuffers(std::vector<std::vector<BulgeInfo>>&& buffers) const {
        DEBUG("Merging bulge buffers");
        perf_counter perf;

        std::vector<BulgeInfo> merged_bulges;
        for (auto& bulge_buffer : buffers) {
            std::copy(std::make_move_iterator(bulge_buffer.begin()),
                      std::make_move_iterator(bulge_buffer.end()),
                      std::back_inserter(merged_bulges));
        }

        DEBUG("Sorting");
        //order is in agreement with coverage
        std::sort(merged_bulges.begin(), merged_bulges.end());
        DEBUG("Total bulges " << merged_bulges.size());
        DEBUG("Buffers merged in " << perf.time() << " seconds");
        return merged_bulges;
    }

    SmartEdgeSet RetainIndependentBulges(std::vector<BulgeInfo>& bulges) const {
        DEBUG("Looking for independent bulges");
        size_t total_cnt = bulges.size();
        perf_counter perf;

        std::vector<BulgeInfo> filtered;
        filtered.reserve(bulges.size());
        //fixme switch to involved vertices to bring fully parallel glueing closer
        std::unordered_set<EdgeId> involved_edges;
        SmartEdgeSet interacting_edges(this->g(), false, CoverageComparator<Graph>(this->g()));

        for (BulgeInfo& info : bulges) {
            TRACE("Analyzing interactions of " << info.str(this->g()));
            if (CheckInteracting(info, involved_edges)) {
                TRACE("Interacting");
                interacting_edges.push(info.e);
            } else {
                TRACE("Independent");
                AccountEdges(info, involved_edges);
                filtered.push_back(std::move(info));
            }
        }
        bulges = std::move(filtered);

        DEBUG("Independent bulges identified in " << perf.time() << " seconds");
        DEBUG("Independent cnt " << bulges.size());
        DEBUG("Interacting cnt " << interacting_edges.size());
        VERIFY(bulges.size() + interacting_edges.size() == total_cnt);

        return interacting_edges;
    }

    size_t ProcessBulges(const std::vector<BulgeInfo>& independent_bulges, SmartEdgeSet&& interacting_edges) {
        DEBUG("Processing bulges");
        perf_counter perf;

        size_t triggered = 0;

        for (const BulgeInfo& info : independent_bulges) {
            TRACE("Processing bulge " << info.str(this->g()));
            triggered++;
            gluer_(info.e, info.alternative);
        }

        DEBUG("Independent bulges glued in " << perf.time() << " seconds");
        perf.reset();

        DEBUG("Processing remaining interacting bulges " << interacting_edges.size());
        //usual br strategy
        for (; !interacting_edges.IsEnd(); ++interacting_edges) {
            EdgeId e = *interacting_edges;
            TRACE("Processing edge " << this->g().str(e));
            std::vector<EdgeId> alternative = alternatives_analyzer_(e);
            if (!alternative.empty()) {
                gluer_(e, alternative);
                triggered++;
            }
        }
        DEBUG("Interacting edges processed in " << perf.time() << " seconds");
        return triggered;
    }

public:

    typedef std::function<void(EdgeId edge, const vector<EdgeId>& path)> BulgeCallbackF;

    ParallelBulgeRemover(Graph& g, const CandidateFinderPtr& interesting_edge_finder,
                         size_t buff_size, double buff_cov_diff,
                         double buff_cov_rel_diff, const AlternativesAnalyzer<Graph>& alternatives_analyzer,
                         BulgeCallbackF opt_callback = 0,
                         std::function<void(EdgeId)> removal_handler = 0,
                         bool track_changes = true) :
                         PersistentAlgorithmBase<Graph>(g),
                         buff_size_(buff_size),
                         buff_cov_diff_(buff_cov_diff),
                         buff_cov_rel_diff_(buff_cov_rel_diff),
                         alternatives_analyzer_(alternatives_analyzer),
                         gluer_(g, opt_callback, removal_handler),
                         interesting_edge_finder_(interesting_edge_finder),
                         tracking_(track_changes),
                         curr_iteration_(0),
                         it_(g, true, CoverageComparator<Graph>(g), true) {
        VERIFY(buff_size_ > 0);
        it_.Detach();
    }

    size_t Run(bool force_primary_launch = false) override {
        bool primary_launch = force_primary_launch ? true : curr_iteration_ == 0;
        //todo remove if not needed;
        //potentially can vary coverage threshold in coordination with ec threshold
        auto proceed_condition = func::AlwaysTrue<EdgeId>();

        if (!it_.IsAttached()) {
            it_.Attach();
        }
        if (primary_launch) {
            it_.clear();
            TRACE("Primary launch.");
            TRACE("Start search for interesting edges");
            interesting_edge_finder_->Run(this->g(), [&](EdgeId e) {it_.push(e);});
            TRACE(it_.size() << " interesting edges to process");
        } else {
            VERIFY(tracking_);
            TRACE(it_.size() << " edges to process");
        }

        size_t triggered = 0;
        bool proceed = true;
        while (proceed) {
            std::vector<EdgeId> edge_buffer;
            edge_buffer.reserve(buff_size_);
            proceed = FillEdgeBuffer(edge_buffer, proceed_condition);

            std::vector<BulgeInfo> bulges = MergeBuffers(FindBulges(edge_buffer));

            auto interacting_edges = RetainIndependentBulges(bulges);

            size_t inner_triggered = ProcessBulges(bulges, std::move(interacting_edges));
            proceed |= (inner_triggered > 0);
            triggered += inner_triggered;
        }

        TRACE("Finished processing. Triggered = " << triggered);
        if (!tracking_)
            it_.Detach();

        curr_iteration_++;

        return triggered;
    }

private:
    DECL_LOGGER("ParallelBulgeRemover")
};

}

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

#include <cmath>
#include <stack>
#include "standard_base.hpp"
#include "omni_utils.hpp"
#include "graph_component.hpp"
#include "xmath.h"
#include "sequence/sequence_tools.hpp"
#include "path_processor.hpp"
#include "graph_processing_algorithm.hpp"

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

	virtual void HandleReversedPath(const vector<EdgeId>& reversed_path) {
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

        //fixme remove after checking results
        bool flag = false;
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
                    flag = true;
                    g_.GlueEdges(split_result.first, path[i]);

                } else {
                    TRACE("GlueEdges " << g_.str(edge_to_split));
                    flag = true;
                    g_.GlueEdges(edge_to_split, path[i]);
                }
            }
            prev_length = bulge_prefix_lengths[i];
        }
        VERIFY(flag);
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
     * and its average coverage is greater than max_relative_coverage_ * g.coverage(e)
     */
    bool BulgeCondition(EdgeId e, const vector<EdgeId>& path,
            double path_coverage) const {
        return math::ge(path_coverage * max_relative_coverage_,
                g_.coverage(e)) && SimplePathCondition<Graph>(g_)(e, path);
    }

public:
    AlternativesAnalyzer(const Graph& g, double max_relative_coverage, size_t max_delta,
                         double max_relative_delta, size_t max_edge_cnt) :
                         g_(g),
                         max_relative_coverage_(max_relative_coverage),
                         max_delta_(max_delta),
                         max_relative_delta_(max_relative_delta),
                         max_edge_cnt_(max_edge_cnt) {
    }

    vector<EdgeId> operator()(EdgeId edge) const {
        size_t kplus_one_mer_coverage = (size_t) math::round((double) g_.length(edge) * g_.coverage(edge));
        TRACE("Processing edge " << g_.str(edge) << " and coverage " << kplus_one_mer_coverage);

        size_t delta = CountMaxDifference(max_delta_, g_.length(edge), max_relative_delta_);

        MostCoveredSimpleAlternativePathChooser<Graph> path_chooser(g_, edge);

        VertexId start = g_.EdgeStart(edge);
        TRACE("Start " << g_.str(start));
        VertexId end = g_.EdgeEnd(edge);
        TRACE("End " << g_.str(end));

        ProcessPaths(g_, (g_.length(edge) > delta) ? g_.length(edge) - delta : 0,
                g_.length(edge) + delta, start, end, path_chooser, max_edge_cnt_);

        const vector<EdgeId>& path = path_chooser.most_covered_path();
        if(path.size() != 0) {
            VERIFY(graph_.EdgeStart(path[0]) == start);
            VERIFY(graph_.EdgeEnd(path.back()) == end);
        }

        double path_coverage = path_chooser.max_coverage();
        if (math::gr(path_coverage, 0.)) {
            TRACE("Best path with coverage " << path_coverage << " is " << PrintPath(g_, path));

            if (BulgeCondition(edge, path, path_coverage)) {
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
};

/**
 * This class removes simple bulges from given graph with the following algorithm: it iterates through all edges of
 * the graph and for each edge checks if this edge is likely to be a simple bulge
 * if edge is judged to be one it is removed.
 */
template<class Graph>
class BulgeRemover: public EdgeProcessingAlgorithm<Graph> {
    typedef EdgeProcessingAlgorithm<Graph> base;
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	bool PossibleBulgeEdge(EdgeId e) const {
	  return (this->g().length(e) <= max_length_ && this->g().coverage(e) < max_coverage_ &&
	          this->g().OutgoingEdgeCount(this->g().EdgeStart(e)) > 1 &&
	          this->g().IncomingEdgeCount(this->g().EdgeEnd(e)) > 1);
	}

protected:

	/*virtual*/
    bool ProcessEdge(EdgeId e) {
        if (this->g().conjugate(e) < e) {
            TRACE("Noncanonical edge");
            VERIFY(false);
            return false;
        }

        TRACE("Considering edge " << this->g().str(e)
                      << " of length " << this->g().length(e)
                      << " and avg coverage " << this->g().coverage(e));
        TRACE("Is possible bulge " << PossibleBulgeEdge(e));

        if (!PossibleBulgeEdge(e)) {
            return false;
        }

        vector<EdgeId> alternative = alternatives_analyzer_(e);
        if (!alternative.empty()) {
            gluer_(e, alternative);
            return true;
        } else {
            return false;
        }
    }

public:

	typedef std::function<void(EdgeId edge, const vector<EdgeId>& path)> BulgeCallbackF;
    
	BulgeRemover(Graph& g, size_t max_length, double max_coverage,
			double max_relative_coverage, size_t max_delta,
			double max_relative_delta,
			size_t max_edge_cnt,
			BulgeCallbackF opt_callback = 0,
			std::function<void(EdgeId)> removal_handler = 0) :
			base(g, true),
			max_length_(max_length),
			max_coverage_(max_coverage),
			alternatives_analyzer_(g, max_relative_coverage, max_delta, 
                                    max_relative_delta, max_edge_cnt),
			gluer_(g, opt_callback, removal_handler) {
                DEBUG("Launching br max_length=" << max_length 
                << " max_coverage=" << max_coverage 
                << " max_relative_coverage=" << max_relative_coverage
                << " max_delta=" << max_delta 
                << " max_relative_delta=" << max_relative_delta
                << " max_number_edges=" << max_edge_cnt);
	}

private:
	size_t max_length_;
	double max_coverage_;
	AlternativesAnalyzer<Graph> alternatives_analyzer_;
	BulgeGluer<Graph> gluer_;
private:
	DECL_LOGGER("BulgeRemover")
};

template<class Graph>
class ParallelBulgeRemover {
    typedef EdgeProcessingAlgorithm<Graph> base;
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

//    typedef std::pair<EdgeId, vector<EdgeId>> BulgeInfo;

    struct BulgeInfo : private boost::noncopyable {
        size_t id;
        EdgeId e;
        std::vector<EdgeId> alternative;

        BulgeInfo() :
            id(-1ul) {
        }

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
            VERIFY(id != that.id);
            return id < that.id;
        }

    };

    bool PossibleBulgeEdge(EdgeId e) const {
      return (g_.length(e) <= max_length_ && g_.coverage(e) < max_coverage_ &&
              g_.OutgoingEdgeCount(g_.EdgeStart(e)) > 1 &&
              g_.IncomingEdgeCount(g_.EdgeEnd(e)) > 1);
    }

    bool CheckInteracting(const BulgeInfo& info, const std::unordered_set<EdgeId>& involved_edges) const {
        if (involved_edges.count(info.e))
            return true;
        for (EdgeId e : info.alternative)
            if (involved_edges.count(e))
                return true;
        return false;
    }

    void AccountEdges(const BulgeInfo& info, std::unordered_set<EdgeId>& involved_edges) const {
        involved_edges.insert(info.e);
        involved_edges.insert(g_.conjugate(info.e));
        for (EdgeId e : info.alternative) {
            involved_edges.insert(e);
            involved_edges.insert(g_.conjugate(e));
        }
    }

public:

    typedef std::function<void(EdgeId edge, const vector<EdgeId>& path)> BulgeCallbackF;
    typedef SmartSetIterator<Graph, EdgeId, CoverageComparator<Graph>> SmartEdgeSet;

    ParallelBulgeRemover(Graph& g, size_t chunk_size,
                         size_t max_length, double max_coverage,
                         double max_relative_coverage, size_t max_delta,
                         double max_relative_delta,
                         size_t max_edge_cnt,
                         BulgeCallbackF opt_callback = 0,
                         std::function<void(EdgeId)> removal_handler = 0) :
                         g_(g),
                         chunk_size_(chunk_size),
                         max_length_(max_length),
                         max_coverage_(max_coverage),
                         alternatives_analyzer_(g, max_relative_coverage, max_delta, 
                                   max_relative_delta, max_edge_cnt),
                         gluer_(g, opt_callback, removal_handler),
                         usual_br_(g, max_length, max_coverage, max_relative_coverage,
                                   max_delta, max_relative_delta, max_edge_cnt, 
                                   opt_callback, removal_handler) {
                             DEBUG("Launching br max_length=" << max_length
                             << " max_coverage=" << max_coverage
                             << " max_relative_coverage=" << max_relative_coverage
                             << " max_delta=" << max_delta
                             << " max_relative_delta=" << max_relative_delta);
                             VERIFY(chunk_size_ > 0);
                 }

    bool ProcessEdge(EdgeId e) {
        if (g_.conjugate(e) < e) {
            TRACE("Noncanonical edge");
            VERIFY(false);
            return false;
        }

        TRACE("Considering edge " << g_.str(e)
                      << " of length " << g_.length(e)
                      << " and avg coverage " << g_.coverage(e));
        TRACE("Is possible bulge " << PossibleBulgeEdge(e));

        if (!PossibleBulgeEdge(e)) {
            return false;
        }

        vector<EdgeId> alternative = alternatives_analyzer_(e);
        if (!alternative.empty()) {
            gluer_(e, alternative);
            return true;
        } else {
            return false;
        }
    }

    //false if time to stop
    template<class SmartEdgeIt>
    bool FillEdgeBuffer(SmartEdgeIt& it, vector<EdgeId>& buffer) const {
        INFO("Filling edge buffer");
        perf_counter perf;
        VERIFY(buffer.empty());
        auto proceed_condition = make_shared<CoverageUpperBound<Graph>>(g_, max_coverage_);

        for (; !it.IsEnd(); ++it) {
            EdgeId e = *it;
            TRACE("Current edge " << g_.str(e));
            if (!proceed_condition->Check(e)) {
                TRACE("Stop condition was reached.");
                //need to release last element of the iterator to make it replacable by new elements
                it.ReleaseCurrent();
                return false;
            }

            if (PossibleBulgeEdge(e)) {
                TRACE("Potential bulge edge");
                buffer.push_back(e);

                if (buffer.size() == chunk_size_) {
                    TRACE("Buffer filled");
                    //need to release last element of the iterator to make it replacable by new elements
                    it.ReleaseCurrent();
                    return true;
                }
            } else {
                TRACE("Can not be bulge");
            }
        }
        INFO("Buffer filled in " << perf.time() << " seconds");
        TRACE("No more edges in iterator");
        return false;
    }
    
    std::vector<std::vector<BulgeInfo>> FindBulges(const std::vector<EdgeId> edge_buffer) const {
    	INFO("Looking for bulges (in parallel)");
        perf_counter perf;
        std::vector<std::vector<BulgeInfo>> bulge_buffers(omp_get_max_threads());
        size_t n = edge_buffer.size();
        #pragma omp parallel for schedule(guided)
        for (size_t i = 0; i < n; ++i) {
            EdgeId e = edge_buffer[i];
            auto alternative = alternatives_analyzer_(e);
            if (!alternative.empty()) {
                bulge_buffers[omp_get_thread_num()].push_back(BulgeInfo(i, e, std::move(alternative)));
            }
        }
        INFO("Buffers found in " << perf.time() << " seconds");
        return bulge_buffers;
    }

    std::vector<BulgeInfo> MergeBuffers(std::vector<std::vector<BulgeInfo>>&& buffers) const {
        INFO("Merging buffers");
        perf_counter perf;

        std::vector<BulgeInfo> merged_bulges;
        for (auto& bulge_buffer : buffers) {
            std::copy(std::make_move_iterator(bulge_buffer.begin()),
                      std::make_move_iterator(bulge_buffer.end()),
                      std::back_inserter(merged_bulges));
        }

        std::sort(merged_bulges.begin(), merged_bulges.end());
        INFO("Buffers merged in " << perf.time() << " seconds");
        return merged_bulges;
    }

    SmartEdgeSet RetainIndependentBulges(std::vector<BulgeInfo>& bulges) const {
        INFO("Looking for independent bulges");
        perf_counter perf;

        std::vector<BulgeInfo> filtered;
        filtered.reserve(bulges.size());
        //fixme switch to involved vertices for fully parallel glueing
        std::unordered_set<EdgeId> involved_edges;
        SmartEdgeSet interacting_edges(g_, CoverageComparator<Graph>(g_));

        for (BulgeInfo& info : bulges) {
            if (CheckInteracting(info, involved_edges)) {
                interacting_edges.push(info.e);
            } else {
                AccountEdges(info, involved_edges);
                filtered.push_back(std::move(info));
            }
        }
        bulges = std::move(filtered);

        INFO("Independent bulges identified in " << perf.time() << " seconds");
        INFO("Independent cnt " << bulges.size());
        INFO("Interacting cnt " << interacting_edges.size());
        return interacting_edges;
    }

    bool ProcessBulges(const std::vector<BulgeInfo>& bulges, SmartEdgeSet& interacting_edges) {
        INFO("Processing bulges");
        perf_counter perf;

    	bool triggered = false;

        for (const BulgeInfo& info : bulges) {
        	triggered = true;
            gluer_(info.e, info.alternative);
        }
        INFO("Independent bulges glued in " << perf.time() << " seconds");
        perf.reset();

        triggered |= usual_br_.RunFromIterator(interacting_edges,
                              make_shared<CoverageUpperBound<Graph>>(g_, max_coverage_));

        INFO("Interacting edges processed in " << perf.time() << " seconds");
        return triggered;
    }

    template<class SmartEdgeIt>
    bool RunFromIterator(SmartEdgeIt& it) {
        VERIFY(it.canonical_only());
        TRACE("Start processing");

        perf_counter perf;

        bool triggered = false;

        bool proceed = true;
        while (proceed) {
            std::vector<EdgeId> edge_buffer;
            edge_buffer.reserve(chunk_size_);
            proceed = FillEdgeBuffer(it, edge_buffer);

            std::vector<BulgeInfo> bulges = MergeBuffers(FindBulges(edge_buffer));

            auto interacting_edges = RetainIndependentBulges(bulges);

            triggered |= ProcessBulges(bulges, interacting_edges);
        }

        TRACE("Finished processing. Triggered = " << triggered);
        return triggered;
    }

    bool Run() {
        auto it = g_.SmartEdgeBegin(CoverageComparator<Graph>(g_), true);
        return RunFromIterator(it);
    }

private:
    Graph& g_;
    size_t chunk_size_;
    size_t max_length_;
    double max_coverage_;
    AlternativesAnalyzer<Graph> alternatives_analyzer_;
    BulgeGluer<Graph> gluer_;
    BulgeRemover<Graph> usual_br_;

private:
    DECL_LOGGER("ParallelBulgeRemover")
};

}

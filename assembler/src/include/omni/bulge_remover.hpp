//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

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
class MostCoveredAlternativePathChooser: public PathProcessor<Graph>::Callback {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	Graph& g_;
	EdgeId forbidden_edge_;
	double max_coverage_;
	vector<EdgeId> most_covered_path_;

public:

	MostCoveredAlternativePathChooser(Graph& g, EdgeId edge) :
			g_(g), forbidden_edge_(edge), max_coverage_(-1.0) {

	}

	virtual void HandleReversedPath(const vector<EdgeId>& reversed_path) {
		vector<EdgeId> path = this->ReversePath(reversed_path);
		double path_cov = AvgCoverage(g_, path);
		for (size_t i = 0; i < path.size(); i++) {
			if (path[i] == forbidden_edge_)
				return;
		}
		if (path_cov > max_coverage_) {
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

/**
 * This class removes simple bulges from given graph with the following algorithm: it iterates through all edges of
 * the graph and for each edge checks if this edge is likely to be a simple bulge
 * if edge is judged to be one it is removed.
 */
template<class Graph>
class BulgeRemover: public EdgeProcessingAlgorithm<Graph, CoverageComparator<Graph>> {
    typedef EdgeProcessingAlgorithm<Graph, CoverageComparator<Graph>> base;
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	bool PossibleBulgeEdge(EdgeId e) const {
	  return (graph_.length(e) <= max_length_ && graph_.coverage(e) < max_coverage_ &&
	          graph_.OutgoingEdgeCount(graph_.EdgeStart(e)) > 1 &&
	          graph_.IncomingEdgeCount(graph_.EdgeEnd(e)) > 1);
	}

	/**
	 * Checks if alternative path is simple (doesn't contain conjugate edges, edge e or conjugate(e))
	 * and its average coverage is greater than max_relative_coverage_ * g.coverage(e)
	 */
	bool BulgeCondition(EdgeId e, const vector<EdgeId>& path,
			double path_coverage) {
		return math::ge(path_coverage * max_relative_coverage_,
				graph_.coverage(e)) && SimplePathCondition<Graph>(graph_)(e, path);
	}

	void ProcessBulge(EdgeId edge, const vector<EdgeId>& path) {
		if (opt_callback_)
			opt_callback_(edge, path);

		if (removal_handler_)
			removal_handler_(edge);

		VertexId start = graph_.EdgeStart(edge);
		VertexId end = graph_.EdgeEnd(edge);

		TRACE("Projecting edge " << graph_.str(edge));
		InnerProcessBulge(edge, path);

		TRACE("Compressing start vertex " << graph_.str(start));
		graph_.CompressVertex(start);

		TRACE("Compressing end vertex " << graph_.str(end));
		graph_.CompressVertex(end);
	}

	void InnerProcessBulge(EdgeId edge, const vector<EdgeId>& path) {

		EnsureEndsPositionAligner aligner(CumulativeLength(graph_, path),
				graph_.length(edge));
		double prefix_length = 0.;
		vector<size_t> bulge_prefix_lengths;

		for (EdgeId e : path) {
			prefix_length += (double) graph_.length(e);
			bulge_prefix_lengths.push_back(aligner.GetPosition((size_t) prefix_length));
		}

		EdgeId edge_to_split = edge;
		size_t prev_length = 0;

		TRACE("Process bulge " << path.size() << " edges");

		//fixme remove after checking results
		bool flag = false;
        VERIFY(bulge_prefix_lengths.back() == graph_.length(edge));

		for (size_t i = 0; i < path.size(); ++i) {
			if (bulge_prefix_lengths[i] > prev_length) {
				if (bulge_prefix_lengths[i] - prev_length
						!= graph_.length(edge_to_split)) {

					TRACE("SplitEdge " << graph_.str(edge_to_split));
					TRACE(
							"Start: " << graph_.str(graph_.EdgeStart(edge_to_split)));
					TRACE(
							"Start: " << graph_.str(graph_.EdgeEnd(edge_to_split)));

					pair<EdgeId, EdgeId> split_result = graph_.SplitEdge(
							edge_to_split,
							bulge_prefix_lengths[i] - prev_length);

					edge_to_split = split_result.second;

					TRACE("GlueEdges " << graph_.str(split_result.first));
					flag = true;
					graph_.GlueEdges(split_result.first, path[i]);

				} else {
					TRACE("GlueEdges " << graph_.str(edge_to_split));
					flag = true;
					graph_.GlueEdges(edge_to_split, path[i]);
				}
			}
			prev_length = bulge_prefix_lengths[i];
		}
		VERIFY(flag);
	}

protected:
    /*virtual*/
    bool ProcessEdge(EdgeId edge) {
        TRACE("Considering edge " << graph_.str(edge) << " of length " << graph_.length(edge) << " and avg coverage " << graph_.coverage(edge));
        TRACE("Is possible bulge " << PossibleBulgeEdge(edge));

        if (!PossibleBulgeEdge(edge)) {
            return false;
        }

        size_t kplus_one_mer_coverage = (size_t) math::round((double) graph_.length(edge) * graph_.coverage(edge));
        TRACE("Processing edge " << graph_.str(edge) << " and coverage " << kplus_one_mer_coverage);

        size_t delta = CountMaxDifference(max_delta_, graph_.length(edge), max_relative_delta_);

        MostCoveredAlternativePathChooser<Graph> path_chooser(graph_, edge);

        VertexId start = graph_.EdgeStart(edge);
        TRACE("Start " << graph_.str(start));

        VertexId end = graph_.EdgeEnd(edge);
        TRACE("End " << graph_.str(end));

        PathProcessor<Graph> path_finder(graph_, (graph_.length(edge) > delta) ? graph_.length(edge) - delta : 0, graph_.length(edge) + delta, start, end, path_chooser);

        path_finder.Process();

        const vector<EdgeId>& path = path_chooser.most_covered_path();

        double path_coverage = path_chooser.max_coverage();
        TRACE("Best path with coverage " << path_coverage << " is " << PrintPath<Graph>(graph_, path));

        if (BulgeCondition(edge, path, path_coverage)) {
            TRACE("Satisfied condition");

            ProcessBulge(edge, path);
            return true;
        } else {
            TRACE("Didn't satisfy condition");
            return false;
        }
    }

public:

	typedef boost::function<void(EdgeId edge, const vector<EdgeId>& path)> BulgeCallbackF;
    
	BulgeRemover(Graph& graph, size_t max_length, double max_coverage,
			double max_relative_coverage, size_t max_delta,
			double max_relative_delta,
			BulgeCallbackF opt_callback = 0,
			boost::function<void(EdgeId)> removal_handler = 0) :
			base(graph,
			     CoverageComparator<Graph>(graph),
			     make_shared<CoverageUpperBound<Graph>>(graph, max_coverage)),
			graph_(graph),
			max_length_(max_length),
			max_coverage_(max_coverage),
			max_relative_coverage_(max_relative_coverage),
			max_delta_(max_delta),
			max_relative_delta_(max_relative_delta),
			opt_callback_(opt_callback),
			removal_handler_(removal_handler) {
	}

//  Old version. If it was math::gr then it would be equivalent to new one.
//	bool RemoveBulges() {
//		bool changed = false;
//		CoverageComparator<Graph> comparator(graph_);
//		for (auto iterator = graph_.SmartEdgeBegin(comparator);
//				!iterator.IsEnd(); ++iterator) {
//			EdgeId e = *iterator;
//			if (math::ge(graph_.coverage(e), max_coverage_))
//				break;
//			changed |= ProcessNext(e);
//		}
//		return changed;
//	}

private:
	Graph& graph_;
	size_t max_length_;
	double max_coverage_;
	double max_relative_coverage_;
	size_t max_delta_;
	double max_relative_delta_;
	BulgeCallbackF opt_callback_;
	boost::function<void(EdgeId)> removal_handler_;

private:
	DECL_LOGGER("BulgeRemover")
};

}

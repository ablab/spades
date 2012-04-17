//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * bulge_remover.hpp
 *
 *  Created on: Apr 13, 2011
 *      Author: sergey
 */

#ifndef BULGE_REMOVER_HPP_
#define BULGE_REMOVER_HPP_

#include <cmath>
#include <stack>
#include "omni_utils.hpp"
#include "xmath.h"
#include "sequence/sequence_tools.hpp"

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

//template<class Graph>
//struct TrivialCondition {
//	typedef typename Graph::EdgeId EdgeId;
//
//	bool operator()(EdgeId edge, const vector<EdgeId>& path) const {
//		for (size_t i = 0; i < path.size(); ++i)
//			for (size_t j = i + 1; j < path.size(); ++j)
//				if (path[i] == path[j])
//					return false;
//		return true;
//	}
//};

template<class Graph>
bool TrivialCondition(typename Graph::EdgeId edge,
		const vector<typename Graph::EdgeId>& path) {
	typedef typename Graph::EdgeId EdgeId;
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

	double PathAvgCoverage(const vector<EdgeId>& path) {
		double unnormalized_coverage = 0;
		size_t path_length = 0;
		for (size_t i = 0; i < path.size(); ++i) {
			EdgeId edge = path[i];
			size_t length = g_.length(edge);
			path_length += length;
			unnormalized_coverage += g_.coverage(edge) * length;
		}
		return unnormalized_coverage / path_length;
	}

public:

	MostCoveredAlternativePathChooser(Graph& g, EdgeId edge) :
			g_(g), forbidden_edge_(edge), max_coverage_(-1.0) {

	}

	virtual void HandlePath(const vector<EdgeId>& path) {
		double path_cov = PathAvgCoverage(path);
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

/**
 * This class removes simple bulges from given graph with the following algorithm: it iterates through all edges of
 * the graph and for each edge checks if this edge is likely to be a simple bulge
 * if edge is judged to be one it is removed.
 */
template<class Graph>
class BulgeRemover {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
public:
	typedef boost::function<bool(EdgeId edge, const vector<EdgeId>& path)> BulgeCallbackF;

	/**
	 * Create BulgeRemover with specified parameters.
	 */
	BulgeRemover(Graph& g, size_t max_length, double max_coverage,
			double max_relative_coverage, double max_delta,
			double max_relative_delta, BulgeCallbackF bulge_condition,
			BulgeCallbackF opt_callback = 0,
			boost::function<void(EdgeId)> removal_handler = 0) :
			g_(g), max_length_(max_length), max_coverage_(max_coverage), max_relative_coverage_(
					max_relative_coverage), max_delta_(max_delta), max_relative_delta_(
					max_relative_delta), bulge_condition_(bulge_condition), opt_callback_(
					opt_callback), removal_handler_(removal_handler) {
	}

	void RemoveBulges();

private:
	Graph& g_;
	size_t max_length_;
	double max_coverage_;
	double max_relative_coverage_;
	double max_delta_;
	double max_relative_delta_;
	BulgeCallbackF bulge_condition_;
	BulgeCallbackF opt_callback_;
	boost::function<void(EdgeId)> removal_handler_;

	bool PossibleBulgeEdge(EdgeId e);

	size_t PathLength(const vector<EdgeId>& path);

	/**
	 * Returns most covered path from start to the end such that its length doesn't exceed length_left.
	 * Returns pair of empty vector and -1 if no such path could be found.
	 * Edges are returned in reverse order!
	 */
	//	pair<vector<EdgeId> , int> BestPath(VertexId start,
	//			VertexId end, int length_left);
	/**
	 * Checks if alternative path is simple (doesn't contain conjugate edges, edge e or conjugate(e))
	 * and its average coverage is greater than max_relative_coverage_ * g.coverage(e)
	 */
	bool BulgeCondition(EdgeId e, const vector<EdgeId>& path,
			double path_coverage) {
		return math::ge(path_coverage * max_relative_coverage_, g_.coverage(e))
				&& bulge_condition_(e, path);
		//		return path_and_coverage.second > max_relative_coverage * g.kplus_one_mer_coverage(edge);
	}

	void ProcessBulge(EdgeId edge, const vector<EdgeId>& path) {
//		UniformPositionAligner aligner(PathLength(path) + 1, g_.length(edge) + 1);
		EnsureEndsPositionAligner aligner(PathLength(path), g_.length(edge));
		double prefix_length = 0.;
		vector<size_t> bulge_prefix_lengths;
		for (auto it = path.begin(); it != path.end(); ++it) {
			prefix_length += g_.length(*it);
			bulge_prefix_lengths.push_back(aligner.GetPosition(prefix_length));
		}
		EdgeId edge_to_split = edge;
		size_t prev_length = 0;
		TRACE("Process bulge " << path.size() << " edges");
		for (size_t i = 0; i < path.size(); ++i) {
			if (bulge_prefix_lengths[i] > prev_length) {
				if (bulge_prefix_lengths[i] - prev_length
						!= g_.length(edge_to_split)) {
					pair<EdgeId, EdgeId> split_result = g_.SplitEdge(
							edge_to_split,
							bulge_prefix_lengths[i] - prev_length);
					edge_to_split = split_result.second;
					g_.GlueEdges(split_result.first, path[i]);
				} else {
					g_.GlueEdges(edge_to_split, path[i]);
				}
			}
			prev_length = bulge_prefix_lengths[i];
		}
	}

private:
	DECL_LOGGER("BulgeRemover")
};

template<class Graph>
bool BulgeRemover<Graph>::PossibleBulgeEdge(EdgeId e) {
	return g_.length(e) <= max_length_ && g_.coverage(e) < max_coverage_;
}

template<class Graph>
size_t BulgeRemover<Graph>::PathLength(const vector<EdgeId>& path) {
	size_t length = 0;
	for (size_t i = 0; i < path.size(); ++i) {
		length += g_.length(path[i]);
	}
	return length;
}

template<class Graph>
void BulgeRemover<Graph>::RemoveBulges() {
	TRACE("Bulge remove process started");

	size_t it_count = 0;

	DEBUG("RemoveBulges function started");

	CoverageComparator<Graph> comparator(g_);
	for (auto iterator = g_.SmartEdgeBegin(comparator); !iterator.IsEnd();
			++iterator) {

		EdgeId edge = *iterator;

		VERBOSE_POWER_T(++it_count, 1000, "th iteration of bulge processing");

		TRACE(
				"Considering edge " << g_.int_id(edge) << " of length " << g_.length(edge)
						<< " and avg coverage " << g_.coverage(edge));
		TRACE("Is possible bulge " << PossibleBulgeEdge(edge));

		if (PossibleBulgeEdge(edge)) {

			size_t kplus_one_mer_coverage = math::round(
					g_.length(edge) * g_.coverage(edge));

			TRACE("Processing edge " << g_.int_id(edge) << " and coverage "
							<< kplus_one_mer_coverage);

			VertexId start = g_.EdgeStart(edge);
			TRACE("Start " << g_.int_id(start));

			VertexId end = g_.EdgeEnd(edge);
			TRACE("End " << g_.int_id(end));
			size_t delta = std::floor(
					std::max(max_relative_delta_ * g_.length(edge),
							max_delta_));
			MostCoveredAlternativePathChooser<Graph> path_chooser(g_, edge);

			PathProcessor<Graph> path_finder(g_,
					(g_.length(edge) > delta) ? g_.length(edge) - delta : 0,
					g_.length(edge) + delta, start, end, path_chooser);

			path_finder.Process();

			const vector<EdgeId>& path = path_chooser.most_covered_path();
			double path_coverage = path_chooser.max_coverage();

			TRACE("Best path with coverage " << path_coverage << " is " << PrintPath<Graph>(g_, path));

			//if edge was returned, this condition will fail
			if (BulgeCondition(edge, path, path_coverage)) {

				TRACE("Satisfied condition");

				if (opt_callback_)
					opt_callback_(edge, path);

				if (removal_handler_)
					removal_handler_(edge);

				TRACE("Projecting edge " << g_.int_id(edge));
				ProcessBulge(edge, path);

				TRACE("Compressing start vertex " << g_.int_id(start))
				g_.CompressVertex(start);

				TRACE("Compressing end vertex " << g_.int_id(end))
				g_.CompressVertex(end);

			} else {
				TRACE("Didn't satisfy condition");
			}
		}
		TRACE("-----------------------------------");
	}
}

}
#endif /* BULGE_REMOVER_HPP_ */

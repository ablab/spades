//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * tip_clipper.hpp
 *
 *  Created on: Mar 25, 2011
 *      Author: sergey
 */

#pragma once

#include <set>

#include "omni_utils.hpp"
#include "xmath.h"
#include "basic_edge_conditions.hpp"

namespace omnigraph {

template<class Graph>
class RelativeCoverageTipCondition: public EdgeCondition<Graph> {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef EdgeCondition<Graph> base;

	const double max_relative_coverage_;

	template<class IteratorType>
	double MaxCompetitorCoverage(EdgeId tip, IteratorType begin, IteratorType end) const {
		double result = 0;
		for (auto it = begin; it != end; ++it) {
			if (*it != tip)
				result = std::max(result, this->g().coverage(*it));
		}
		return result;
	}

	double MaxCompetitorCoverage(EdgeId tip) const {
		const Graph &g = this->g();
		VertexId start = g.EdgeStart(tip), end = g.EdgeEnd(tip);
		auto out = g.OutgoingEdges(start);
		auto in = g.IncomingEdges(end);
		return std::max(
						MaxCompetitorCoverage(tip, out.begin(),	out.end()),
						MaxCompetitorCoverage(tip, in.begin(), in.end()));
//		return std::max(
//				MaxCompetitorCoverage(tip, g.out_begin(start),
//						g.out_end(start)),
//				MaxCompetitorCoverage(tip, g.in_begin(end), g.in_end(end)));
	}

public:

	RelativeCoverageTipCondition(const Graph& g, double max_relative_coverage) :
			base(g), max_relative_coverage_(max_relative_coverage) {
	}

	bool Check(EdgeId e) const {
		//+1 is a trick to deal with edges of 0 coverage from iterative run
		double max_coverage = MaxCompetitorCoverage(e) + 1;
		return math::le(this->g().coverage(e),
				max_relative_coverage_ * max_coverage);
	}
};

//todo refactor
/**
 * This class removes tips from given graph with the following algorithm: it iterates through all edges of
 * the graph(in order defined by certain comparator) and for each edge checks if this edge is likely to be
 * a tip and if edge is judged to be one it is removed.
 * todo should extend EdgeProcessingAlgorithm
 */
template<class Graph>
class TipClipper: private boost::noncopyable {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef function<bool(EdgeId)> edge_condition_t;

	Graph &graph_;
	shared_ptr<Predicate<EdgeId>> additional_condition_;
	size_t removed_;

public:
	const size_t max_tip_length_;

private:
	boost::function<void(EdgeId)> removal_handler_;

protected:

	const Graph& graph() const {
		return graph_;
	}

	Graph& graph() {
		return graph_;
	}

	/**
	 * This method checks if given vertex topologically looks like end of tip
	 * @param v vertex to be checked
	 * @return true if vertex judged to be tip and false otherwise.
	 */
	bool IsTip(VertexId v) const {
		return graph_.IncomingEdgeCount(v) + graph_.OutgoingEdgeCount(v) == 1;
	}

	/**
	 * This method checks if given edge topologically looks like a tip.
	 * @param edge edge vertex to be checked
	 * @return true if edge judged to be tip and false otherwise.
	 */
	bool IsTip(EdgeId edge) const {
		return graph_.length(edge) <= max_tip_length_
				&& (IsTip(graph_.EdgeEnd(edge)) || IsTip(graph_.EdgeStart(edge)))
				&& (graph_.OutgoingEdgeCount(graph_.EdgeStart(edge))
						+ graph_.IncomingEdgeCount(graph_.EdgeEnd(edge)) > 2);
	}

	void CompressSplitVertex(VertexId splitVertex) {
		if (graph_.CanCompressVertex(splitVertex)) {
			graph_.CompressVertex(splitVertex);
		}
	}

	bool DeleteTipVertex(VertexId vertex) {
		if (graph_.IsDeadEnd(vertex) && graph_.IsDeadStart(vertex)) {
			graph_.DeleteVertex(vertex);
			return true;
		}
		return false;
	}

	void ProcessVertex(VertexId v) {
		if (!DeleteTipVertex(v)) {
			CompressSplitVertex(v);
		}
	}

	virtual void RemoveTip(EdgeId tip) {
		VertexId start = graph_.EdgeStart(tip);
		VertexId end = graph_.EdgeEnd(tip);
		if (removal_handler_) {
			removal_handler_(tip);
		}
		graph_.DeleteEdge(tip);
		ProcessVertex(start);
		ProcessVertex(end);
	}

	virtual bool TryToRemoveTip(EdgeId tip) {
		RemoveTip(tip);
		TRACE("Edge removed");
		return true;
	}

	/**
	 * Method clips tips of the graph.
	 */
	bool ProcessNext(const EdgeId& tip) {
		TRACE("Checking edge for being tip " << this->graph().str(tip));
		if (this->IsTip(tip)) {
			TRACE("Edge " << this->graph().str(tip) << " judged to look like a tip topologically");
			if (additional_condition_->Check(tip)) {
				TRACE("Edge " << this->graph().str(tip) << " judged to be a tip");

				if (this->TryToRemoveTip(tip)) {
					removed_++;
				} else {
					return false;
				}

			} else {
				TRACE("Edge " << this->graph().str(tip) << " judged NOT to be tip");
			}
		} else {
			TRACE("Edge " << this->graph().str(tip) << " judged NOT to look like tip topologically");
		}
		return true;
	}

public:

	/**
	 * Create TipClipper with specified parameters. Those parameters could probably be replaced later with
	 * certain generic checker class.
	 */
	TipClipper(Graph &graph, size_t max_tip_length,
			const shared_ptr<Predicate<EdgeId>>& additional_condition,
			boost::function<void(EdgeId)> removal_handler = 0,
			boost::function<double(EdgeId)> qual_f = 0) :
			graph_(graph), additional_condition_(additional_condition), removed_(
					0), max_tip_length_(max_tip_length), removal_handler_(
					removal_handler) {

	}

	bool ClipTips() {
		LengthComparator<Graph> comparator(graph_);
		for (auto iterator = graph_.SmartEdgeBegin(comparator); !iterator.IsEnd();
				++iterator) {
			EdgeId e = *iterator;
			if (graph_.length(e) > max_tip_length_)
				break;
			this->ProcessNext(e);
		}
		return removed_ > 0;
	}

private:
	DECL_LOGGER("AbstractTipClipper")
};

template<class Graph>
class DefaultTipClipper: public TipClipper<Graph> {

	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef TipClipper<Graph> base;

public:

	DefaultTipClipper(Graph &graph, size_t max_tip_length, size_t max_coverage,
			double max_relative_coverage,
			boost::function<void(EdgeId)> removal_handler = 0,
			boost::function<double(EdgeId)> qual_f = 0) :
			base(graph, max_tip_length,
					And<EdgeId>(
							make_shared<CoverageUpperBound<Graph>>(graph,
									max_coverage),
							make_shared<RelativeCoverageTipCondition<Graph>>(
									graph, max_relative_coverage)),
					removal_handler, qual_f) {
	}

private:
	DECL_LOGGER("DefaultTipClipper")
};

template<class Graph>
class TopologyTipClipper: public TipClipper<Graph> {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef TipClipper<Graph> base;

public:
	TopologyTipClipper(Graph &graph, size_t max_tip_length,
			size_t uniqueness_length, size_t plausibility_length,
			boost::function<void(EdgeId)> removal_handler = 0) :
			base(graph, max_tip_length,
					make_shared<DefaultUniquenessPlausabilityCondition<Graph>>(graph, uniqueness_length,
							plausibility_length), removal_handler) {
	}

private:
	DECL_LOGGER("TopologyTipClipper")
};

} // namespace omnigraph

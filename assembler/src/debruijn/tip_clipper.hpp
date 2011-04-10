/*
 * tip_clipper.hpp
 *
 *  Created on: Mar 25, 2011
 *      Author: sergey
 */

#ifndef TIP_CLIPPER_HPP_
#define TIP_CLIPPER_HPP_

#include <set>
#include "edge_graph.hpp"
#include "set"
#include "utils.hpp"

#define DEFAULT_COVERAGE_BOUND 5
#define DEFAULT_RELATIVE_COVERAGE_BOUND 2.

namespace edge_graph {

using de_bruijn::PriorityQueue;

struct TipComparator {
private:
	const EdgeGraph &graph_;
public:
	TipComparator(EdgeGraph &graph) :
		graph_(graph) {
	}

	bool operator()(Edge* edge1, Edge* edge2) const {
		return graph_.EdgeNucls(edge1).size() < graph_.EdgeNucls(edge2).size();
	}
};

template<typename Comparator>
class TipClipper {
private:
	EdgeGraph &graph_;
	PriorityQueue<Edge *, Comparator> tipQueue_;
	const size_t maxTipLength_;
	const size_t coverageBound_;
	const double relativeCoverageBound_;

	bool isTip(Vertex *v) {
		if (!graph_.CheckUniqueIncomingEdge(v) || !graph_.IsDeadEnd(v))
			return false;
		Edge *edge = graph_.GetUniqueIncomingEdge(v);
		return graph_.length(edge) <= maxTipLength_;
	}

	bool isTip(Edge *edge) {
		return isTip(graph_.EdgeEnd(edge));
	}

	void FindTips() {
		for (EdgeGraph::VertexIterator it = graph_.begin(); it
				!= graph_.begin(); ++it) {
			if (isTip(*it)) {
				tipQueue_.offer(graph_.GetUniqueIncomingEdge(*it));
			}
		}
	}

	size_t maxCompetotorCoverage(Vertex *splitVertex, Edge *tip) {
		assert(!graph_.CheckUniqueOutgiongEdge(splitVertex));
		if (graph_.CheckUniqueOutgiongEdge(splitVertex)) {
			assert(false);//such situation should never occur
		}
		const vector<Edge *> competitors = graph_.OutgoingEdges(splitVertex);
		size_t result = 0;
		for (vector<Edge *>::const_iterator it = competitors.begin(); it
				!= competitors.end(); ++it) {
			if (*it != tip)
				result = max(result, graph_.coverage(*it));
		}
		return result;
	}

	bool tipShouldBeRemoved(Edge *tip) {
		if (graph_.length(tip) > maxTipLength_ || graph_.coverage(tip)
				> coverageBound_)
			return false;
		Vertex *splitVertex = graph_.EdgeStart(tip);
		if (graph_.CheckUniqueOutgiongEdge(splitVertex))
			return false;
		size_t maxCoverage = maxCompetotorCoverage(splitVertex, tip);
		return graph_.coverage(tip) <= relativeCoverageBound_ * maxCoverage;
	}

	void compressSplitVertex(Vertex *splitVertex) {
		if (graph_.CanCompressVertex(splitVertex)) {
			Edge *edge1 = graph_.GetUniqueOutgoingEdge(splitVertex);
			Edge *edge2 = graph_.GetUniqueOutgoingEdge(
					graph_.Complement(splitVertex));
			if (isTip(edge1) || isTip(edge2)) {
				graph_.CompressVertex(splitVertex);
			}
		}
	}

	//	void compressSplitVertex(Vertex *splitVertex) {
	//		if (graph_.CanCompressVertex(splitVertex)) {
	//			graph_.CompressVertex(splitVertex);
	//		}
	//	}

	void removeTip(Edge *tip) {
		Vertex *splitVertex = graph_.EdgeStart(tip);
		Vertex *tipVertex = graph_.EdgeEnd(tip);
		graph_.DeleteEdge(tip);
		graph_.DeleteVertex(tipVertex);
		compressSplitVertex(splitVertex);
	}

	void RemoveTips() {
		while (!tipQueue_.empty()) {
			Edge * tip = tipQueue_.poll();
			if (tipShouldBeRemoved(tip)) {
				removeTip(tip);
			}
		}
	}

public:
	TipClipper(EdgeGraph &graph, Comparator comparator, size_t maxTipLength,
			size_t coverageBound,
			double relativeCoverageBound = DEFAULT_RELATIVE_COVERAGE_BOUND) :
		graph_(graph), tipQueue_(comparator), maxTipLength_(maxTipLength),
				coverageBound_(coverageBound),
				relativeCoverageBound_(coverageBound) {
	}

	TipClipper(EdgeGraph &graph, Comparator comparator) :
		graph_(graph), tipQueue_(comparator), maxTipLength_(2 * graph.k() - 1),
				coverageBound_(DEFAULT_COVERAGE_BOUND),
				relativeCoverageBound_(DEFAULT_RELATIVE_COVERAGE_BOUND) {
	}

	void ClipTips() {
		FindTips();
		RemoveTips();
		graph_.CompressAllVertices();
	}

};

}

#endif /* TIP_CLIPPER_HPP_ */

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
#define DEFAULT_MAX_TIP_LENGTH 50

namespace edge_graph {

using de_bruijn::PriorityQueue;

struct TipComparator {
private:
	EdgeGraph *graph_;
public:
	TipComparator() {
		assert(false);
	}

	TipComparator(EdgeGraph &graph) :
		graph_(&graph) {
	}

	bool operator()(Edge* edge1, Edge* edge2) const {
		return graph_->EdgeNucls(edge1).size() < graph_->EdgeNucls(edge2).size();
	}
};

template<typename Comparator>
class TipClipper {
private:
//	EdgeGraph *graph_;
//	PriorityQueue<Edge *, Comparator> tipQueue_;
	const size_t maxTipLength_;
	const size_t coverageBound_;
	const double relativeCoverageBound_;
	Comparator comparator_;

	bool isTip(EdgeGraph &graph, Vertex *v) {
		if (!graph.CheckUniqueIncomingEdge(v) || !graph.IsDeadEnd(v))
			return false;
		Edge *edge = graph.GetUniqueIncomingEdge(v);
		return graph.length(edge) <= maxTipLength_;
	}

	bool isTip(EdgeGraph &graph, Edge *edge) {
		return isTip(graph, graph.EdgeEnd(edge));
	}

//	void FindTips() {
//		for (EdgeGraph::VertexIterator it = graph_.begin(); it
//				!= graph_.begin(); ++it) {
//			if (isTip(*it)) {
//				tipQueue_.offer(graph_.GetUniqueIncomingEdge(*it));
//			}
//		}
//	}

	size_t maxCompetitorCoverage(EdgeGraph &graph, Vertex *splitVertex, Edge *tip) {
		assert(!graph.CheckUniqueOutgiongEdge(splitVertex));
		if (graph.CheckUniqueOutgiongEdge(splitVertex)) {
			assert(false);//such situation should never occur
		}
		const vector<Edge *> competitors = graph.OutgoingEdges(splitVertex);
		size_t result = 0;
		for (vector<Edge *>::const_iterator it = competitors.begin(); it
				!= competitors.end(); ++it) {
			if (*it != tip)
				result = max(result, graph.coverage(*it));
		}
		return result;
	}

	bool tipShouldBeRemoved(EdgeGraph &graph, Edge *tip) {
		if (graph.length(tip) > maxTipLength_ || graph.coverage(tip)
				> coverageBound_)
			return false;
		Vertex *splitVertex = graph.EdgeStart(tip);
		if (graph.CheckUniqueOutgiongEdge(splitVertex))
			return false;
		size_t maxCoverage = maxCompetitorCoverage(graph, splitVertex, tip);
		return graph.coverage(tip) <= relativeCoverageBound_ * maxCoverage;
	}

	void compressSplitVertex(EdgeGraph &graph, Vertex *splitVertex) {
		if (graph.CanCompressVertex(splitVertex)) {
			Edge *edge1 = graph.GetUniqueOutgoingEdge(splitVertex);
			Edge *edge2 = graph.GetUniqueOutgoingEdge(
					graph.Complement(splitVertex));
			if (isTip(graph, edge1) || isTip(graph, edge2)) {
				graph.CompressVertex(splitVertex);
			}
		}
	}

	//	void compressSplitVertex(Vertex *splitVertex) {
	//		if (graph_.CanCompressVertex(splitVertex)) {
	//			graph_.CompressVertex(splitVertex);
	//		}
	//	}

	void removeTip(EdgeGraph &graph, Edge *tip) {
		Vertex *splitVertex = graph.EdgeStart(tip);
		Vertex *tipVertex = graph.EdgeEnd(tip);
		graph.DeleteEdge(tip);
		graph.DeleteVertex(tipVertex);
		compressSplitVertex(graph, splitVertex);
	}

//	void RemoveTips() {
//		while (!tipQueue_.empty()) {
//			Edge * tip = tipQueue_.poll();
//			if (tipShouldBeRemoved(tip)) {
//				removeTip(tip);
//			}
//		}
//	}

public:
	TipClipper(Comparator comparator, size_t maxTipLength,
			size_t coverageBound,
			double relativeCoverageBound = DEFAULT_RELATIVE_COVERAGE_BOUND) : comparator_(comparator), maxTipLength_(maxTipLength),
				coverageBound_(coverageBound),
				relativeCoverageBound_(coverageBound) {
	}

	TipClipper(Comparator comparator) : comparator_(comparator), maxTipLength_(DEFAULT_MAX_TIP_LENGTH),
				coverageBound_(DEFAULT_COVERAGE_BOUND),
				relativeCoverageBound_(DEFAULT_RELATIVE_COVERAGE_BOUND) {
	}

	void ClipTips(EdgeGraph &graph) {
		PriorityQueue<Edge *, Comparator> tipQueue(comparator_);
		de_bruijn::SmartEdgeIterator<EdgeGraph, Comparator> iterator(graph);
		de_bruijn::SmartEdgeIterator<EdgeGraph, Comparator> end;
		while(end != iterator) {
			Edge * tip = *iterator;
			if (tipShouldBeRemoved(graph, tip)) {
				removeTip(graph, tip);
			}
			++iterator;
		}
//		FindTips(tipQueue);
//		RemoveTips();
		graph.CompressAllVertices();
	}

};

}

#endif /* TIP_CLIPPER_HPP_ */

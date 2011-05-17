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
#include "utils.hpp"
#include "omni_utils.hpp"
#include "omni_tools.hpp"

#define DEFAULT_COVERAGE_BOUND 1000
#define DEFAULT_RELATIVE_COVERAGE_BOUND 2.
#define DEFAULT_MAX_TIP_LENGTH 50

namespace edge_graph {

using omnigraph::Compresser;

template<class Graph>
struct TipComparator {
private:
	Graph *graph_;
public:
	TipComparator() {
		assert(false);
	}

	TipComparator(Graph &graph) :
		graph_(&graph) {
	}

	bool operator()(Edge* edge1, Edge* edge2) const {
		if (graph_->EdgeNucls(edge1).size() == graph_->EdgeNucls(edge2).size()) {
			return edge1 < edge2;
		}
		return graph_->EdgeNucls(edge1).size()
				< graph_->EdgeNucls(edge2).size();
	}
};

template<class Graph, typename Comparator>
class TipClipper {
private:
	//	Graph *graph_;
	//	PriorityQueue<Edge *, Comparator> tipQueue_;
	Comparator comparator_;
	const size_t maxTipLength_;
	const size_t coverageBound_;
	const double relativeCoverageBound_;

	bool isTip(Graph &graph, Vertex *v) {
		if (!graph.CheckUniqueIncomingEdge(v) || !graph.IsDeadEnd(v))
			return false;
		Edge *edge = graph.GetUniqueIncomingEdge(v);
		return graph.length(edge) <= maxTipLength_;
	}

	bool isTip(Graph &graph, Edge *edge) {
		return isTip(graph, graph.EdgeEnd(edge));
	}

	//	void FindTips() {
	//		for (Graph::VertexIterator it = graph_.begin(); it
	//				!= graph_.begin(); ++it) {
	//			if (isTip(*it)) {
	//				tipQueue_.offer(graph_.GetUniqueIncomingEdge(*it));
	//			}
	//		}
	//	}

	double maxCompetitorCoverage(Graph &graph, Vertex *splitVertex,
			Edge *tip) {
		assert(!graph.CheckUniqueOutgiongEdge(splitVertex));
		if (graph.CheckUniqueOutgiongEdge(splitVertex)) {
			assert(false);//such situation should never occur
		}
		const vector<Edge *> competitors = graph.OutgoingEdges(splitVertex);
		double result = 0;
		for (vector<Edge *>::const_iterator it = competitors.begin(); it
				!= competitors.end(); ++it) {
			if (*it != tip)
				result = max(result, graph.coverage(*it));
		}
		return result;
	}

	bool tipShouldBeRemoved(Graph &graph, Edge *tip) {
		if (graph.length(tip) > maxTipLength_ || graph.coverage(tip)
				> coverageBound_)
			return false;
		Vertex *splitVertex = graph.EdgeStart(tip);
		if (graph.CheckUniqueOutgiongEdge(splitVertex))
			return false;
		double maxCoverage = maxCompetitorCoverage(graph, splitVertex, tip);
		return graph.coverage(tip) <= relativeCoverageBound_ * maxCoverage;
	}

	void compressSplitVertex(Graph &graph, Vertex *splitVertex) {
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
	//			graph_.CompressVertex(s	plitVertex);
	//		}
	//	}

	void removeTip(Graph &graph, Edge *tip) {
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
			double relativeCoverageBound = DEFAULT_RELATIVE_COVERAGE_BOUND) :
		comparator_(comparator), maxTipLength_(maxTipLength),
				coverageBound_(coverageBound),
				relativeCoverageBound_(coverageBound) {
	}

	TipClipper(Comparator comparator) :
		comparator_(comparator), maxTipLength_(DEFAULT_MAX_TIP_LENGTH),
				coverageBound_(DEFAULT_COVERAGE_BOUND),
				relativeCoverageBound_(DEFAULT_RELATIVE_COVERAGE_BOUND) {
	}

	void ClipTips(Graph &graph) {
		SmartEdgeIterator<Graph, Comparator> iterator =
				graph.SmartEdgeBegin(comparator_);
		SmartEdgeIterator<Graph, Comparator> end =
				graph.SmartEdgeEnd(comparator_);
		while (end != iterator) {
			EdgeId tip = *iterator;
			if (isTip(graph, tip)) {
				bool tmp = tipShouldBeRemoved(graph, tip);
				if (tmp) {
					removeTip(graph, tip);
				}
			}
			++iterator;
		}
		Compresser<Graph> compresser(graph);
		compresser.CompressAllVertices();
	}

};

}

#endif /* TIP_CLIPPER_HPP_ */

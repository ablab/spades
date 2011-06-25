/*
 * tip_clipper.hpp
 *
 *  Created on: Mar 25, 2011
 *      Author: sergey
 */

#ifndef TIP_CLIPPER_HPP_
#define TIP_CLIPPER_HPP_

#include <set>
//#include "edge_graph.hpp"
//#include "utils.hpp"
#include "omni_utils.hpp"
#include "omni_tools.hpp"
//
//#define DEFAULT_COVERAGE_BOUND 1000
//#define DEFAULT_RELATIVE_COVERAGE_BOUND 2.0
//#define DEFAULT_MAX_TIP_LENGTH 50

namespace omnigraph {

/**
 * This class defines which edge is more likely to be tip. In this case we just assume shorter edges
 * are more likely tips then longer ones.
 */
template<class Graph>
struct TipComparator {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	Graph *graph_;
public:
	/**
	 * TipComparator should never be created with default constructor but it is necessary on order for
	 * code to compile.
	 */
	//	TipComparator() {
	//		assert(false);
	//	}

	/**
	 * Construct TipComparator for given graph
	 * @param graph graph for which comparator is created
	 */
	TipComparator(Graph &graph) :
		graph_(&graph) {
	}

	/**
	 * Standard comparator function as used in collections.
	 */
	bool operator()(EdgeId edge1, EdgeId edge2) const {
		if (graph_->EdgeNucls(edge1).size() == graph_->EdgeNucls(edge2).size()) {
			return edge1 < edge2;
		}
		return graph_->EdgeNucls(edge1).size()
				< graph_->EdgeNucls(edge2).size();
	}
};

/**
 * This class removes tips from given graph with the following algorithm: it iterates through all edges of
 * the graph(in order defined by certain comparator) and for each edge checks if this edge is likely to be
 * a tip and if edge is judged to be one it is removed.
 */
template<class Graph, typename Comparator>
class TipClipper {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	const size_t max_tip_length_;
	const size_t max_coverage_;
	const double max_relative_coverage_;

	Graph &graph_;
	Comparator comparator_;

	/**
	 * This method checks if given vertex topologically looks like end of tip
	 * @param v vertex to be checked
	 * @return true if vertex judged to be tip and false otherwise.
	 */
	bool IsTip(VertexId v) {
		return graph_.IncomingEdgeCount(v) + graph_.OutgoingEdgeCount(v) == 1;
	}

	/**
	 * This method checks if given edge topologically looks like a tip.
	 * @param edge edge vertex to be checked
	 * @return true if edge judged to be tip and false otherwise.
	 */
	bool IsTip(EdgeId edge) {
		return graph_.length(edge) <= max_tip_length_ && (IsTip(
				graph_.EdgeEnd(edge)) || IsTip(graph_.EdgeStart(edge)));
	}

	//	void FindTips() {
	//		for (Graph::VertexIterator it = graph_.begin(); it
	//				!= graph_.begin(); ++it) {
	//			if (isTip(*it)) {
	//				tipQueue_.offer(graph_.GetUniqueIncomingEdge(*it));
	//			}
	//		}
	//	}

	double MaxCompetitorCoverage(EdgeId tip, vector<EdgeId> competitors) {
		double result = 0;
		for (auto it = competitors.begin(); it != competitors.end(); ++it) {
			if (*it != tip)
				result = max(result, graph_.coverage(*it));
		}
		return result;
	}

	double MaxCompetitorCoverage(EdgeId tip) {
		return max(
				MaxCompetitorCoverage(tip,
						graph_.OutgoingEdges(graph_.EdgeStart(tip))),
				MaxCompetitorCoverage(tip,
						graph_.IncomingEdges(graph_.EdgeEnd(tip))));
	}

	/**
	 * This method checks if given edge is a tip and thus should be removed
	 * @param tip edge to check
	 */
	bool TipShouldBeRemoved(EdgeId tip) {
		if (graph_.length(tip) > max_tip_length_ || graph_.coverage(tip)
				> max_coverage_)
			return false;
		//		VertexId splitVertex = graph_.EdgeStart(tip);
		if (graph_.OutgoingEdgeCount(graph_.EdgeStart(tip))
				+ graph_.IncomingEdgeCount(graph_.EdgeEnd(tip)) == 2)
			return false;
		double max_coverage = MaxCompetitorCoverage(tip);
		return graph_.coverage(tip) <= max_relative_coverage_ * max_coverage;
	}

	void CompressSplitVertex(VertexId splitVertex) {
		if (graph_.CanCompressVertex(splitVertex)) {
			EdgeId edge1 = graph_.GetUniqueOutgoingEdge(splitVertex);
			EdgeId edge2 = graph_.GetUniqueIncomingEdge(splitVertex);
			if (IsTip(edge1) || IsTip(edge2)) {
//				if (graph_.CanCompressVertex(splitVertex)
				graph_.CompressVertex(splitVertex);
			}
		}
	}

	void DeleteTipVertex(VertexId vertex) {
		if (graph_.IsDeadEnd(vertex) && graph_.IsDeadStart(vertex)) {
			TRACE("DELETE IT")
			graph_.DeleteVertex(vertex);
			TRACE("DELETE OK")
		}
	}

	//	void compressSplitVertex(Vertex *splitVertex) {
	//		if (graph_.CanCompressVertex(splitVertex)) {
	//			graph_.CompressVertex(s	plitVertex);
	//		}
	//	}


	void removeTip(EdgeId tip) {
		TRACE("remove Tip 1")
		VertexId start = graph_.EdgeStart(tip);
		TRACE("remove Tip 2")
		VertexId end = graph_.EdgeEnd(tip);
		TRACE("remove Tip 3")
		graph_.DeleteEdge(tip);
		TRACE("remove Tip 4")
		if (graph_.CanCompressVertex(start))
				CompressSplitVertex(start);
		else
			DeleteTipVertex(start);
		if (graph_.CanCompressVertex(end))
				CompressSplitVertex(end);
		else
			DeleteTipVertex(end);
		TRACE("remove Tip 8")
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

	/**
	 * Create TipClipper with specified parameters. Those parameters could probably be replaced later with
	 * certain generic checker class.
	 */
	TipClipper(Graph &graph, Comparator comparator, size_t max_tip_length,
			size_t max_coverage, double max_relative_coverage) :
		max_tip_length_(max_tip_length), max_coverage_(max_coverage),
				max_relative_coverage_(max_relative_coverage), graph_(graph),
				comparator_(comparator) {
	}

	/**
	 * Method clips tips of the graph.
	 */
	void ClipTips() {
		TRACE("Tip clipping started");
		for (auto iterator = graph_.SmartEdgeBegin(comparator_); !iterator.IsEnd(); ) {
			EdgeId tip = *iterator;
			TRACE("Checking edge for being tip " << tip);
			if (IsTip(tip)) {
				TRACE("Edge " << tip << " judged to look like tip topologically");
				if (TipShouldBeRemoved(tip)) {
					TRACE("Edge " << tip << " judged to be tip");
					removeTip(tip);
					TRACE("Edge " << tip << " removed as tip");
				} else {
					TRACE("Edge " << tip << " judged NOT to be tip");
				}
			} else {
				TRACE("Edge " << tip << " judged NOT to look like tip topologically");
			}
			TRACE("Try to find next edge");
			++iterator;
			TRACE("Use next edge");
		}
		TRACE("Tip clipping finished");
		Compressor<Graph> compresser(graph_);
		compresser.CompressAllVertices();
	}

private:
	DECL_LOGGER("TipClipper")
};

}

#endif /* TIP_CLIPPER_HPP_ */

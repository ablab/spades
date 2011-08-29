/*
 * erroneous_connection_remover.hpp
 *
 *  Created on: May 31, 2011
 *      Author: sergey
 */

#ifndef ERRONEOUS_CONNECTION_REMOVER_HPP_
#define ERRONEOUS_CONNECTION_REMOVER_HPP_

namespace omnigraph {

#include "omni_tools.hpp"
#include "xmath.h"

template<class Graph>
struct CoverageComparator {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	const Graph& graph_;
public:
	CoverageComparator(const Graph &graph) :
		graph_(graph) {
	}

	/**
	 * Standard comparator function as used in collections.
	 */
	bool operator()(EdgeId edge1, EdgeId edge2) const {
		if (math::eq(graph_.coverage(edge1), graph_.coverage(edge2))) {
			return edge1 < edge2;
		}
		return math::ls(graph_.coverage(edge1), graph_.coverage(edge2));
	}
};

template <class Graph>
class LowCoverageEdgeRemover {
	size_t max_length_;
	double max_coverage_;
public:
	LowCoverageEdgeRemover(size_t max_length, double max_coverage)
	: max_length_(max_length), max_coverage_(max_coverage) {

	}

	void RemoveEdges(Graph& g) {
		for (auto it = g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
			typename Graph::EdgeId e = *it;
			if (g.length(e) < max_length_ && g.coverage(e) < max_coverage_) {
				g.DeleteEdge(e);
			}
		}
		omnigraph::Compressor<Graph> compressor(g);
		compressor.CompressAllVertices();
		omnigraph::Cleaner<Graph> cleaner(g);
		cleaner.Clean();
	}

};

template <class Graph>
class IterativeLowCoverageEdgeRemover {
	size_t max_length_;
	double max_coverage_;
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

public:
	IterativeLowCoverageEdgeRemover(size_t max_length, double max_coverage)
	: max_length_(max_length), max_coverage_(max_coverage) {

	}

	void RemoveEdges(Graph& g) {
		CoverageComparator<Graph> comparator(g);
		for (auto it = g.SmartEdgeBegin(comparator); !it.IsEnd(); ++it) {
			typename Graph::EdgeId e = *it;
			if (math::gr(g.coverage(e), max_coverage_)) {
				return;
			}
			if (g.length(e) < max_length_) {
				VertexId start = g.EdgeStart(e);
				VertexId end = g.EdgeEnd(e);
				g.DeleteEdge(e);
				g.CompressVertex(start);
				g.CompressVertex(end);
			}
		}
		omnigraph::Cleaner<Graph> cleaner(g);
		cleaner.Clean();
	}

};

}

#endif /* ERRONEOUS_CONNECTION_REMOVER_HPP_ */

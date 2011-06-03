/*
 * erroneous_connection_remover.hpp
 *
 *  Created on: May 31, 2011
 *      Author: sergey
 */

#ifndef ERRONEOUS_CONNECTION_REMOVER_HPP_
#define ERRONEOUS_CONNECTION_REMOVER_HPP_

namespace debruijn_graph {

#include "omni_tools.hpp"

template <class Graph>
class LowCoverageEdgeRemover {
	size_t max_length_;
	double max_coverage_;
public:
	LowCoverageEdgeRemover(size_t max_length, double max_coverage)
	: max_length_(max_length), max_coverage_(max_coverage) {

	}

	void RemoveEdges(Graph& g) {
		for (auto it = g.SmartEdgeBegin(); !it.isEnd(); ++it) {
			typename Graph::EdgeId e = *it;
			if (g.length(e) < max_length_ && g.coverage(e) < max_coverage_) {
				g.DeleteEdge(e);
			}
		}
		omnigraph::Compressor<Graph> compressor(g);
		compressor.CompressAllVertices();
	}

};

}

#endif /* ERRONEOUS_CONNECTION_REMOVER_HPP_ */

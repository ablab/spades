/*
 * erroneous_connection_remover.hpp
 *
 *  Created on: May 31, 2011
 *      Author: sergey
 */

#ifndef ERRONEOUS_CONNECTION_REMOVER_HPP_
#define ERRONEOUS_CONNECTION_REMOVER_HPP_

namespace debruijn_graph {

template <class Graph>
class LowCoverageEdgeRemover {
	size_t length_density_;
	double coverage_density_;
public:
	LowCoverageEdgeRemover(size_t length_density, double coverage_density)
	: length_density_(length_density), coverage_density_(coverage_density) {

	}

	void RemoveEdges(Graph& g) {
		for (auto it = g.SmartEdgeBegin(); !it.isEnd(); ++it) {
			typename Graph::EdgeId e = *it;
			if (g.length(e) < length_density_ && g.coverage(e) < coverage_density_) {
				g.DeleteEdge(e);
			}
		}
	}

};

}

#endif /* ERRONEOUS_CONNECTION_REMOVER_HPP_ */

/*
 * simplification_quality_counter.hpp
 *
 *  Created on: May 31, 2011
 *      Author: sergey
 */

#ifndef SIMPLIFICATION_QUALITY_COUNTER_HPP_
#define SIMPLIFICATION_QUALITY_COUNTER_HPP_

#include "utils.hpp"

namespace debruijn_graph {

template<size_t k, class Graph>
class SimplificationQualityCounter {
	typedef Graph::EdgeId EdgeId;
public:
	size_t CountBlackEdges(Graph& g, EdgeIndex<k + 1, Graph>& index, const string& genome) {
		size_t black_count = 0;
		debruijn_graph::SimpleSequenceMapper<k, Graph> sequence_mapper(g, index);
		Path<EdgeId> path = sequence_mapper.MapSequence(Sequence(genome));
		const vector<EdgeId> path_edges = path.sequence();
		set<EdgeId> colored_edges(path_edges.first(), path_edges.second());
		for (auto it = g.SmartEdgeBegin(); it != g.SmartEdgeEnd(); ++it) {
			if (colored_edges.count(*it) == 0) {
				black_count++;
			}
		}
		return black_count;
	}
};

}

#endif /* SIMPLIFICATION_QUALITY_COUNTER_HPP_ */

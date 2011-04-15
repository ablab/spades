/*
 * bulge_remover.hpp
 *
 *  Created on: Apr 13, 2011
 *      Author: sergey
 */

#ifndef BULGE_REMOVER_HPP_
#define BULGE_REMOVER_HPP_

namespace de_bruijn {
template<class Graph>
class BulgeRemover {
public:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	void RemoveBulges(Graph& g);
private:
	bool PossibleBulgeEdge(const Graph& g, EdgeId e);

	//returns reversed path!!!
	pair<vector<EdgeId> , int> BestPath(const Graph& g, VertexId start,
			VertexId end, size_t length_left);
};

template<class Graph>
bool BulgeRemover<Graph>::PossibleBulgeEdge(const Graph& g, EdgeId e) {
	//todo
	size_t max_length = 15;
	return g.length(e) <= max_length;
}

template<class Graph>
void BulgeRemover<Graph>::RemoveBulges(Graph& g) {
	typedef de_bruijn::SmartEdgeIterator<Graph> EdgeIter;

	//todo
	size_t delta = 3;
	double coverage_gap = 2.;

	assert(coverage_gap > 1.);

	for (EdgeIter iterator = g.SmartEdgeBegin(), end = g.SmartEdgeEnd(); end != iterator; ++iterator) {
		EdgeId edge = *iterator;
		if (PossibleBulgeEdge(g, edge)) {
			VertexId start = g.EdgeStart(edge);
			VertexId end = g.EdgeEnd(edge);
			pair<vector<EdgeId> , int> path_and_coverage = BestPath(g, start,
					end, g.length(edge) + delta);

			//if edge was returned, this condition will fail
			//it will not! coverage_gap > 1. And also g.coverage(edge) is average coverage while coverage_gap is coverage sum!
			if (path_and_coverage.second  < coverage_gap * g.coverage(edge)) {
				g.DeleteEdge(edge);
				g.CompressVertex(start);
				g.CompressVertex(end);
			}
		}
	}
}

template<class Graph>
pair<vector<typename Graph::EdgeId> , int> BulgeRemover<Graph>::BestPath(
		const Graph& g, typename Graph::VertexId start,
		typename Graph::VertexId end, size_t length_left) {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	if (length_left < 0) {
		return make_pair(vector<EdgeId> (0), -1);
	}
	if (start == end) {
		return make_pair(vector<EdgeId> (), 0);
	}
	vector<EdgeId> outgoing_edges = g.OutgoingEdges(start);
	int best_path_coverage = -1;
	vector<EdgeId> best_path(0);
	for (size_t i = 0; i < outgoing_edges.size(); ++i) {
		EdgeId edge = outgoing_edges[i];
		size_t coverage = g.coverage(edge);
		pair<vector<EdgeId> , int> path_and_coverage = BestPath(g,
				g.EdgeEnd(edge), end, length_left - g.length(edge));
		if (path_and_coverage.second >= 0 && path_and_coverage.second + (int) coverage > best_path_coverage) {
			best_path_coverage = path_and_coverage.second + coverage;
			best_path = path_and_coverage.first;
			best_path.push_back(edge);
		}
	}
	return make_pair(best_path, best_path_coverage);
}

}
#endif /* BULGE_REMOVER_HPP_ */

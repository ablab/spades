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
			VertexId end, int length_left);
};

template<class Graph>
bool BulgeRemover<Graph>::PossibleBulgeEdge(const Graph& g, EdgeId e) {
	//todo
	size_t max_length = 30;
	return g.length(e) <= max_length;
}

template<class Graph>
const string PrintPath(Graph& g, const vector<typename Graph::EdgeId>& edges) {
	string delim = "";
	stringstream ss;
	for (size_t i = 0; i < edges.size(); ++i) {
		ss << delim << g.EdgeNucls(edges[i]).str() << " (" + g.length(edges[i]) << ")";
		delim = " -> ";
	}
	return ss.str();
}

template<class Graph>
void BulgeRemover<Graph>::RemoveBulges(Graph& g) {
	typedef de_bruijn::SmartEdgeIterator<Graph> EdgeIter;

	//todo
	size_t delta = 3;
	double coverage_gap = 1.1;

	assert(coverage_gap > 1.);

	DEBUG("Bulge remove process started");

	for (EdgeIter iterator = g.SmartEdgeBegin(), end = g.SmartEdgeEnd(); end != iterator; ++iterator) {
		EdgeId edge = *iterator;
		DEBUG("Considering edge of length " << g.length(edge) << " and avg coverage " << g.coverage(edge));
		DEBUG("Is possible bulge " << PossibleBulgeEdge(g, edge));
		if (PossibleBulgeEdge(g, edge)) {
			DEBUG("Processing edge " << g.EdgeNucls(edge) << " and coverage " << g.kplus_one_mer_coverage(edge));

			VertexId start = g.EdgeStart(edge);
			DEBUG("Start " << g.VertexNucls(start));

			VertexId end = g.EdgeEnd(edge);
			DEBUG("End " << g.VertexNucls(end));

			DEBUG("Looking for best path between start and end shorter than " << g.length(edge) + delta);

			pair<vector<EdgeId> , int> path_and_coverage = BestPath(g, start,
					end, g.length(edge) + delta);
			DEBUG("Best path with coverage " << path_and_coverage.second << " is " << PrintPath<Graph>(g, path_and_coverage.first));


			//if edge was returned, this condition will fail
			DEBUG("Satisfies condition " << (path_and_coverage.second  > coverage_gap * g.kplus_one_mer_coverage(edge)));
			if (path_and_coverage.second  > coverage_gap * g.kplus_one_mer_coverage(edge)) {
				DEBUG("Deleting edge " << g.EdgeNucls(edge) << " and compressing ends");
				g.DeleteEdge(edge);
				g.CompressVertex(start);
				g.CompressVertex(end);
			}
		}
		DEBUG("-----------------------------------");
	}
}

template<class Graph>
pair<vector<typename Graph::EdgeId> , int> BulgeRemover<Graph>::BestPath(
		const Graph& g, typename Graph::VertexId start,
		typename Graph::VertexId end, int length_left) {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	DEBUG("Recursive call for vertex " << g.VertexNucls(start) << " path length left " << length_left);
	if (length_left < 0) {
		DEBUG("Length threshold was exceeded, returning no path");
		return make_pair(vector<EdgeId> (0), -1);
	}
	if (start == end) {
		DEBUG("Path found, backtracking");
		return make_pair(vector<EdgeId> (), 0);
	}
	vector<EdgeId> outgoing_edges = g.OutgoingEdges(start);
	int best_path_coverage = -1;
	vector<EdgeId> best_path(0);
	DEBUG("Iterating through outgoing edges, finding best path");
	for (size_t i = 0; i < outgoing_edges.size(); ++i) {
		EdgeId edge = outgoing_edges[i];
		DEBUG("Going along edge " << g.EdgeNucls(edge));
		size_t kplus_one_mer_coverage = g.kplus_one_mer_coverage(edge);
		pair<vector<EdgeId> , int> path_and_coverage = BestPath(g,
				g.EdgeEnd(edge), end, length_left - g.length(edge));
		if (path_and_coverage.second >= 0 && path_and_coverage.second + (int) kplus_one_mer_coverage > best_path_coverage) {
			best_path_coverage = path_and_coverage.second + kplus_one_mer_coverage;
			best_path = path_and_coverage.first;
			best_path.push_back(edge);
		}
	}
	DEBUG("Best path from vertex "<< g.VertexNucls(start) << " is " << PrintPath(g, best_path) << " with coverage " << best_path_coverage);
	return make_pair(best_path, best_path_coverage);
}

}
#endif /* BULGE_REMOVER_HPP_ */

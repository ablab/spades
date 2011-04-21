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

	BulgeRemover(size_t max_length = 5 * K, double coverage_threshold = 1000.,
			double coverage_gap = 1.1, double delta_threshold = 4, double relative_delta = 0.1) :
		max_length_(max_length), coverage_threshold_(coverage_threshold),
				coverage_gap_(coverage_gap), delta_threshold_(delta_threshold), relative_delta_(relative_delta) {
		assert(coverage_gap_ > 1.);
	}

private:

	size_t max_length_;
	double coverage_threshold_;
	double coverage_gap_;
	double delta_threshold_;
	double relative_delta_;

	bool PossibleBulgeEdge(const Graph& g, EdgeId e);

	bool PathCondition(const Graph& g, EdgeId edge, const vector<EdgeId>& path);

	bool BulgeCondition(const Graph& g, EdgeId edge, const vector<EdgeId>& path) {
		return PathAvgCoverage(g, path) > coverage_gap_ * g.coverage(edge) && PathCondition(g, edge, path);
		//		return path_and_coverage.second > coverage_gap * g.kplus_one_mer_coverage(edge);
	}

	double PathAvgCoverage(const Graph& g, const vector<EdgeId> path);

	size_t PathLength(const Graph& g, const vector<EdgeId> path);

	//returns reversed path!!!
	pair<vector<EdgeId> , int> BestPath(const Graph& g, VertexId start,
			VertexId end, int length_left);
};

template<class Graph>
bool BulgeRemover<Graph>::PossibleBulgeEdge(const Graph& g, EdgeId e) {
	return g.length(e) <= max_length_ && g.coverage(e) < coverage_threshold_;
}

template<class Graph>
size_t BulgeRemover<Graph>::PathLength(const Graph& g, const vector<EdgeId> path) {
	size_t length = 0;
	for (size_t i = 0; i < path.size(); ++i) {
		length += g.length(path[i]);
	}
	return length;
}

template<class Graph>
bool BulgeRemover<Graph>::PathCondition(const Graph& g, EdgeId edge, const vector<EdgeId>& path) {
	for (size_t i = 0; i < path.size(); ++i)
		if (edge == path[i] || edge == g.Complement(path[i]))
			return false;
	for (size_t i = 0; i < path.size(); ++i)
		for (size_t j = i + 1; j < path.size(); ++j)
			if (path[i] == path[j] || path[i] == g.Complement(path[j]))
				return false;
	return true;
}

template<class Graph>
double BulgeRemover<Graph>::PathAvgCoverage(const Graph& g, const vector<EdgeId> path) {
	double unnormalized_coverage = 0;
	size_t path_length = 0;
	for (size_t i = 0; i < path.size(); ++i) {
		EdgeId edge = path[i];
		size_t length = g.length(edge);
		path_length += length;
		unnormalized_coverage += g.coverage(edge) * length;
	}
	return unnormalized_coverage / path_length;
}

//todo move to some common place
template<class Graph>
const string PrintPath(Graph& g, const vector<typename Graph::EdgeId>& edges) {
	string delim = "";
	stringstream ss;
	for (size_t i = 0; i < edges.size(); ++i) {
		ss << delim << g.EdgeNucls(edges[i]).str() << " ("
				<< g.length(edges[i]) << ")";
		delim = " -> ";
	}
	return ss.str();
}

template<class Graph>
void BulgeRemover<Graph>::RemoveBulges(Graph& g) {
	typedef de_bruijn::SmartEdgeIterator<Graph> EdgeIter;

	DEBUG("Bulge remove process started");

	for (EdgeIter iterator = g.SmartEdgeBegin(), end = g.SmartEdgeEnd(); end
			!= iterator; ++iterator) {
		EdgeId edge = *iterator;
		DEBUG(
				"Considering edge of length " << g.length(edge)
						<< " and avg coverage " << g.coverage(edge));
		DEBUG("Is possible bulge " << PossibleBulgeEdge(g, edge));
		if (PossibleBulgeEdge(g, edge)) {
			DEBUG("Processing edge " << g.EdgeNucls(edge) << " and coverage "
							<< g.KPlusOneMerCoverage(edge));

			VertexId start = g.EdgeStart(edge);
			DEBUG("Start " << g.VertexNucls(start));

			VertexId end = g.EdgeEnd(edge);
			DEBUG("End " << g.VertexNucls(end));

			size_t length_threshold = g.length(edge) + std::max(relative_delta_ * g.length(edge), delta_threshold_);
			DEBUG("Looking for best path between start and end shorter than " << length_threshold);

			pair<vector<EdgeId> , int> path_and_coverage = BestPath(g, start,
					end, length_threshold);
			const vector<EdgeId>& path = path_and_coverage.first;
			double path_coverage = (double)path_and_coverage.second / PathLength(g, path);
			DEBUG("Best path with coverage " << path_coverage
							<< " is " << PrintPath<Graph> (g, path));

			//if edge was returned, this condition will fail
			if (BulgeCondition(g, edge, path)) {
				DEBUG("Satisfied condition");
				DEBUG("Deleting edge " << g.EdgeNucls(edge) << " and compressing ends");
				g.DeleteEdge(edge);
				g.CompressVertex(start);
				g.CompressVertex(end);
			} else {
				DEBUG("Didn't satisfy condition");
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

	DEBUG(
			"Recursive call for vertex " << g.VertexNucls(start)
					<< " path length left " << length_left);
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
		DEBUG("Going along edge of length " << g.length(edge));
		size_t kplus_one_mer_coverage = g.KPlusOneMerCoverage(edge);
		pair<vector<EdgeId> , int> path_and_coverage = BestPath(g,
				g.EdgeEnd(edge), end, length_left - g.length(edge));
		if (path_and_coverage.second >= 0 && path_and_coverage.second
				+ (int) kplus_one_mer_coverage > best_path_coverage) {
			best_path_coverage = path_and_coverage.second
					+ kplus_one_mer_coverage;
			best_path = path_and_coverage.first;
			best_path.push_back(edge);
		}
	}
	DEBUG(
			"Best path from vertex " << g.VertexNucls(start) << " is "
					<< PrintPath(g, best_path) << " with coverage "
					<< best_path_coverage);
	return make_pair(best_path, best_path_coverage);
}

}
#endif /* BULGE_REMOVER_HPP_ */

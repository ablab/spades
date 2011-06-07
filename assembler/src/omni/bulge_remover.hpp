/*
 * bulge_remover.hpp
 *
 *  Created on: Apr 13, 2011
 *      Author: sergey
 */

#ifndef BULGE_REMOVER_HPP_
#define BULGE_REMOVER_HPP_

#include <cmath>
#include "omni_utils.hpp"

namespace omnigraph {

/**
 * This class removes simple bulges from given graph with the following algorithm: it iterates through all edges of
 * the graph and for each edge checks if this edge is likely to be a simple bulge
 * if edge is judged to be one it is removed.
 */
template<class Graph>
class BulgeRemover {
public:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	/**
	 * Create BulgeRemover with specified parameters.
	 */
	BulgeRemover(size_t max_length, double max_coverage,
			double max_relative_coverage, double max_delta,
			double max_relative_delta) :
		max_length_(max_length), max_coverage_(max_coverage),
				max_relative_coverage_(max_relative_coverage),
				max_delta_(max_delta), max_relative_delta_(max_relative_delta) {
		assert(max_relative_coverage_ > 1.);
	}

	void RemoveBulges(Graph& g);

private:

	size_t max_length_;
	double max_coverage_;
	double max_relative_coverage_;
	double max_delta_;
	double max_relative_delta_;

	bool PossibleBulgeEdge(const Graph& g, EdgeId e);

	bool SimplePathCondition(const Graph& g, EdgeId edge, const vector<EdgeId>& path);

	/**
	 * Checks if alternative path is simple (doesn't contain conjugate edges, edge e or conjugate(e))
	 * and its average coverage is greater than max_relative_coverage_ * g.coverage(e)
	 */
	bool BulgeCondition(const Graph& g, EdgeId e, const vector<EdgeId>& path) {
		return PathAvgCoverage(g, path) > max_relative_coverage_ * g.coverage(
				e) && SimplePathCondition(g, e, path);
		//		return path_and_coverage.second > max_relative_coverage * g.kplus_one_mer_coverage(edge);
	}

	double PathAvgCoverage(const Graph& g, const vector<EdgeId> path);

	size_t PathLength(const Graph& g, const vector<EdgeId> path);

	/**
	 * Returns most covered path from start to the end such that its length doesn't exceed length_left.
	 * Returns pair of empty vector and -1 if no such path could be found.
	 * Edges are returned in reverse order!
	 */
	pair<vector<EdgeId> , int> BestPath(const Graph& g, VertexId start,
			VertexId end, int length_left);
};

template<class Graph>
bool BulgeRemover<Graph>::PossibleBulgeEdge(const Graph& g, EdgeId e) {
	return g.length(e) <= max_length_ && g.coverage(e) < max_coverage_;
}

template<class Graph>
size_t BulgeRemover<Graph>::PathLength(const Graph& g,
		const vector<EdgeId> path) {
	size_t length = 0;
	for (size_t i = 0; i < path.size(); ++i) {
		length += g.length(path[i]);
	}
	return length;
}

template<class Graph>
bool BulgeRemover<Graph>::SimplePathCondition(const Graph& g, EdgeId edge,
		const vector<EdgeId>& path) {
	for (size_t i = 0; i < path.size(); ++i)
		if (edge == path[i] || edge == g.conjugate(path[i]))
			return false;
	for (size_t i = 0; i < path.size(); ++i) {
		if (path[i] == g.conjugate(path[i])) {
			return false;
		}
		for (size_t j = i + 1; j < path.size(); ++j)
			if (path[i] == path[j] || path[i] == g.conjugate(path[j]))
				return false;
	}
	return true;
}

template<class Graph>
double BulgeRemover<Graph>::PathAvgCoverage(const Graph& g,
		const vector<EdgeId> path) {
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
	//todo add right comparator
	typedef SmartEdgeIterator<Graph> EdgeIter;

	TRACE("Bulge remove process started");

	for (EdgeIter iterator = g.SmartEdgeBegin(); !iterator.IsEnd(); ++iterator) {
		EdgeId edge = *iterator;
		TRACE(
				"Considering edge of length " << g.length(edge)
				<< " and avg coverage " << g.coverage(edge));
		TRACE("Is possible bulge " << PossibleBulgeEdge(g, edge));
		if (PossibleBulgeEdge(g, edge)) {
			size_t kplus_one_mer_coverage = std::floor(
					g.length(edge) * g.coverage(edge) + 0.5);
			TRACE(
					"Processing edge " << g.EdgeNucls(edge) << " and coverage "
					<< kplus_one_mer_coverage);

			VertexId start = g.EdgeStart(edge);
			TRACE("Start " << g.VertexNucls(start));

			VertexId end = g.EdgeEnd(edge);
			TRACE("End " << g.VertexNucls(end));

			size_t length_threshold = g.length(edge) + std::max(
					max_relative_delta_ * g.length(edge), max_delta_);
			TRACE(
					"Looking for best path between start and end shorter than "
					<< length_threshold);

			pair<vector<EdgeId> , int> path_and_coverage = BestPath(g, start,
					end, length_threshold);
			const vector<EdgeId>& path = path_and_coverage.first;
			double path_coverage = (double) path_and_coverage.second
					/ PathLength(g, path);
			TRACE(
					"Best path with coverage " << path_coverage << " is "
					<< PrintPath<Graph> (g, path));

			//if edge was returned, this condition will fail
			if (BulgeCondition(g, edge, path)) {
				TRACE("Satisfied condition");
				TRACE(
						"Deleting edge " << g.EdgeNucls(edge)
						<< " and compressing ends");
				g.DeleteEdge(edge);
				g.CompressVertex(start);
				g.CompressVertex(end);
			} else {
				TRACE("Didn't satisfy condition");
			}
		}
		TRACE("-----------------------------------");
	}
}

template<class Graph>
pair<vector<typename Graph::EdgeId> , int> BulgeRemover<Graph>::BestPath(
		const Graph& g, typename Graph::VertexId start,
		typename Graph::VertexId end, int length_left) {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	TRACE(
			"Recursive call for vertex " << g.VertexNucls(start)
			<< " path length left " << length_left);
	if (length_left < 0) {
		TRACE("Length threshold was exceeded, returning no path");
		return make_pair(vector<EdgeId> (0), -1);
	}
	if (start == end) {
		TRACE("Path found, backtracking");
		return make_pair(vector<EdgeId> (), 0);
	}
	vector<EdgeId> outgoing_edges = g.OutgoingEdges(start);
	int best_path_coverage = -1;
	vector<EdgeId> best_path(0);
	TRACE("Iterating through outgoing edges, finding best path");
	for (size_t i = 0; i < outgoing_edges.size(); ++i) {
		EdgeId edge = outgoing_edges[i];
		TRACE("Going along edge of length " << g.length(edge));
		size_t kplus_one_mer_coverage = std::floor(
				g.length(edge) * g.coverage(edge) + 0.5);//kplus_one_mer_coverage(edge);
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
	TRACE(
			"Best path from vertex " << g.VertexNucls(start) << " is "
			<< PrintPath(g, best_path) << " with coverage "
			<< best_path_coverage);
	return make_pair(best_path, best_path_coverage);
}

}
#endif /* BULGE_REMOVER_HPP_ */

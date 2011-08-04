/*
 * seeds.hpp
 *
 *  Created on: Aug 3, 2011
 *      Author: andrey
 */

#ifndef SEEDS_HPP_
#define SEEDS_HPP_

#include "lc_common.hpp"

namespace long_contigs {

using namespace debruijn_graph;

// ====== Seed functions ======

//Extends trivial path forward
//If a start of another trivial path is found, returns it
//Otherwise returns 0
EdgeId ExtendTrivialForward(Graph& g, BidirectionalPath& path, const std::map<EdgeId, BidirectionalPath>& starts) {
	if (path.empty()) {
		return 0;
	}

	VertexId currentVertex = g.EdgeEnd(path.back());
	while (g.CheckUniqueOutgoingEdge(currentVertex)) {
		EdgeId nextEdge = g.GetUniqueOutgoingEdge(currentVertex);

		if (starts.count(nextEdge) != 0) {
			return nextEdge;
		}

		path.push_back(nextEdge);
		currentVertex = g.EdgeEnd(nextEdge);
	}
	return 0;
}

//Previous one without checking for other seeds' starts
EdgeId ExtendTrivialForward(Graph& g, BidirectionalPath& path) {
	static std::map<EdgeId, BidirectionalPath> empty = std::map<EdgeId, BidirectionalPath>();
	return ExtendTrivialForward(g, path, empty);
}


//Trivially extend path backward
void ExtendTrivialBackward(Graph& g, BidirectionalPath& path) {
	if (path.empty()) {
		return;
	}

	VertexId currentVertex = g.EdgeStart(path.front());
	while (g.CheckUniqueIncomingEdge(currentVertex)) {
		EdgeId nextEdge = g.GetUniqueIncomingEdge(currentVertex);
		path.push_front(nextEdge);
		currentVertex = g.EdgeStart(nextEdge);
	}
}

//Glue second path to the first one
void JoinPaths(BidirectionalPath& path1, BidirectionalPath& path2) {
	if (path1 == path2) {
		INFO("Cannot join path with itself");
		return;
	}

	for (auto iter = path2.begin(); iter != path2.end(); ++iter) {
		path1.push_back(*iter);
	}
}

//Find all seeds as trivial paths
void FindSeeds(Graph& g, std::vector<BidirectionalPath>& seeds) {
	std::set<EdgeId> visited;
	std::map<EdgeId, BidirectionalPath> starts;
	int count = 0;

	INFO("Finding seeds started");
	for (auto iter = g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
		count++;
		EdgeId e = *iter;

		//Try to make seed starting from unvisited edge
		if (visited.count(e) == 0) {
			visited.insert(e);

			starts[e] = BidirectionalPath();
			BidirectionalPath& newPath = starts[e];
			newPath.push_back(e);

			//Extend trivially
			EdgeId nextStart = ExtendTrivialForward(g, newPath, starts);
			for (auto edgeInPath = newPath.begin(); edgeInPath != newPath.end(); ++edgeInPath) {
				visited.insert(*edgeInPath);
			}

			//If extended till another seed, than concatenate them
			if (nextStart != 0) {
				//INFO("Join paths");
				JoinPaths(newPath, starts[nextStart]);
				starts.erase(nextStart);
			}
		}
	}

	//Debug part
	for (auto pathIter = starts.begin(); pathIter != starts.end(); ++pathIter) {
		seeds.push_back(pathIter->second);
	}
	PrintPathCoverage(g, seeds);
	//End of debug part

	//Extending seed backward
	seeds.clear();
	seeds.reserve(starts.size());
	INFO("Extending seeds backward");
	for (auto pathIter = starts.begin(); pathIter != starts.end(); ++pathIter) {
		ExtendTrivialBackward(g, pathIter->second);
		seeds.push_back(pathIter->second);
	}
	INFO("Finding seeds finished");
}


} // namespace long_contigs


#endif /* SEEDS_HPP_ */

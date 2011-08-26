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
	static bool glueSeeds = lc_cfg::get().ss.glue_seeds;
	if (path.empty()) {
		return 0;
	}

	VertexId currentVertex = g.EdgeEnd(path.back());
	//TODO cycling -- tmp fix
	while (g.CheckUniqueOutgoingEdge(currentVertex) && path.size() < 10) {
		EdgeId nextEdge = g.GetUniqueOutgoingEdge(currentVertex);

		if (glueSeeds && starts.count(nextEdge) != 0) {
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
	//TODO fix cycling
	while (g.CheckUniqueIncomingEdge(currentVertex) && path.size() < 12) {
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
	std::map<EdgeId, BidirectionalPath> starts;
	int count = 0;

	INFO("Finding seeds started");
	for (auto iter = g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
		INFO(count);
		count++;
		EdgeId e = *iter;

		starts[e] = BidirectionalPath();
//		seeds.push_back(BidirectionalPath());
		BidirectionalPath& newPath = starts[e];
		newPath.push_back(e);

		//Extend trivially
		EdgeId nextStart = ExtendTrivialForward(g, newPath, starts);

		//If extended till another seed, than concatenate them
		if (nextStart != 0) {
			//INFO("Join paths");
			JoinPaths(newPath, starts[nextStart]);
			starts.erase(nextStart);
		}
	}

	//Debug part
//	for (auto pathIter = starts.begin(); pathIter != starts.end(); ++pathIter) {
//		seeds.push_back(pathIter->second);
//	}
//	PrintPathCoverage(g, seeds);
	//End of debug part

	//Extending seed backward
	seeds.clear();
	seeds.reserve(starts.size());
	INFO("Extending seeds backward");
	for (auto pathIter = starts.begin(); pathIter != starts.end(); ++pathIter) {
		ExtendTrivialBackward(g, pathIter->second);
		seeds.push_back(pathIter->second);
	}
//	for (auto pathIter = seeds.begin(); pathIter != seeds.end(); ++pathIter) {
//		ExtendTrivialBackward(g, *pathIter);
//	}

	INFO("Finding seeds finished");
}

} // namespace long_contigs


#endif /* SEEDS_HPP_ */

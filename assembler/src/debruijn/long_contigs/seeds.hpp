/*
 * seeds.hpp
 *
 *  Created on: Aug 3, 2011
 *      Author: andrey
 */

#ifndef SEEDS_HPP_
#define SEEDS_HPP_

#include "lc_common.hpp"
#include "path_utils.hpp"
#include "loop.hpp"


namespace long_contigs {

using namespace debruijn_graph;

// ====== Seed functions ======

//Extends trivial path forward
//If a start of another trivial path is found, returns it
//Otherwise returns 0
EdgeId ExtendTrivialForward(Graph& g, BidirectionalPath& path, LoopDetector& detector, const std::map<EdgeId, BidirectionalPath>& starts,
		PathLengths* lengths = 0) {
	static bool glueSeeds = lc_cfg::get().ss.glue_seeds;
	static bool maxCycles = lc_cfg::get().ss.max_cycles;

	if (path.empty()) {
		return 0;
	}

	VertexId currentVertex = g.EdgeEnd(path.back());
	while (g.CheckUniqueOutgoingEdge(currentVertex)) {
		EdgeId nextEdge = g.GetUniqueOutgoingEdge(currentVertex);

		if (glueSeeds && starts.count(nextEdge) != 0) {
			return nextEdge;
		}

		path.push_back(nextEdge);
		if (lengths != 0) {
			IncreaseLengths(g, *lengths, nextEdge, true);
		}
		currentVertex = g.EdgeEnd(nextEdge);

		detector.temp.clear();
		detector.temp.AddAlternative(nextEdge);
		detector.AddNewEdge(nextEdge, path.size() - 1);
		if (CheckCycle(path, nextEdge, detector, maxCycles)) {
			break;
		}
	}
	return 0;
}

//Previous one without checking for other seeds' starts
EdgeId ExtendTrivialForward(Graph& g, BidirectionalPath& path, LoopDetector& detector, PathLengths* lengths = 0) {
	static std::map<EdgeId, BidirectionalPath> empty = std::map<EdgeId, BidirectionalPath>();
	return ExtendTrivialForward(g, path, detector, empty, lengths);
}


//Trivially extend path backward
void ExtendTrivialBackward(Graph& g, BidirectionalPath& path, LoopDetector& detector, PathLengths* lengths = 0) {
	static bool maxCycles = lc_cfg::get().ss.max_cycles;

	if (path.empty()) {
		return;
	}

	VertexId currentVertex = g.EdgeStart(path.front());
	while (g.CheckUniqueIncomingEdge(currentVertex)) {
		EdgeId nextEdge = g.GetUniqueIncomingEdge(currentVertex);

		path.push_front(nextEdge);
		if (lengths != 0) {
			IncreaseLengths(g, *lengths, nextEdge, false);
		}
		currentVertex = g.EdgeStart(nextEdge);

		detector.temp.clear();
		detector.temp.AddAlternative(nextEdge);
		detector.AddNewEdge(nextEdge, path.size() - 1);
		if (CheckCycle(path, nextEdge, detector, maxCycles)) {
			break;
		}
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

void SimpleRecountDetectorForward(BidirectionalPath& path, LoopDetector& detector) {
	detector.clear();

	for (int i = 0; i < (int) path.size(); ++i) {
		detector.temp.clear();
		detector.temp.AddAlternative(path[i]);
		detector.AddNewEdge(path[i], i);
	}
}

void SimpleRecountDetectorBackward(BidirectionalPath& path, LoopDetector& detector) {
	detector.clear();

	for (int i = path.size() - 1; i >= 0; --i) {
		detector.temp.clear();
		detector.temp.AddAlternative(path[i]);
		detector.AddNewEdge(path[i], path.size() - 1 - i);
	}
}


//Find all seeds as trivial paths
void FindSeeds(Graph& g, std::vector<BidirectionalPath>& seeds) {
	std::map<EdgeId, BidirectionalPath> starts;
	LoopDetector detector;
	int count = 0;

	INFO("Finding seeds started");
	for (auto iter = g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
		count++;
		EdgeId e = *iter;

		detector.clear();
		detector.temp.clear();
		detector.temp.AddAlternative(e);
		detector.AddNewEdge(e, 0);

		starts[e] = BidirectionalPath();
		BidirectionalPath& newPath = starts[e];
		newPath.push_back(e);

		//Extend trivially
		EdgeId nextStart = ExtendTrivialForward(g, newPath, detector, starts);

		//If extended till another seed, than concatenate them
		if (nextStart != 0) {
			JoinPaths(newPath, starts[nextStart]);
			starts.erase(nextStart);
		}
	}

	//Extending seed backward
	seeds.clear();
	seeds.reserve(starts.size());

	INFO("Extending seeds backward");
	for (auto pathIter = starts.begin(); pathIter != starts.end(); ++pathIter) {
		SimpleRecountDetectorBackward(pathIter->second, detector);
		ExtendTrivialBackward(g, pathIter->second, detector);
		seeds.push_back(pathIter->second);
	}

	INFO("Finding seeds finished");
}

} // namespace long_contigs


#endif /* SEEDS_HPP_ */

/*
 * long_contigs.hpp
 *
 *  Created on: Jul 8, 2011
 *      Author: andrey
 */

#ifndef LONG_CONTIGS_HPP_
#define LONG_CONTIGS_HPP_

#include <vector>
#include <deque>
#include <set>
#include <map>

#include "../launch.hpp"
#include "../../omni/paired_info.hpp"

namespace long_contigs {

using namespace debruijn_graph;

//Heuristic constants here
static const int DISTANCE_DIV = 1;

typedef std::deque<EdgeId> BiderectionalPath;
typedef std::deque<double> PathLengths;

//Extends trivial path forward
//If a start of another trivial path is found, returns it
//Otherwise returns 0
EdgeId ExtendTrivialForward(Graph& g, BiderectionalPath& path, std::map<EdgeId, BiderectionalPath> starts) {
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


//Trivially extend path backward
void ExtendTrivialBackward(Graph& g, BiderectionalPath& path) {
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
void JoinPaths(BiderectionalPath& path1, BiderectionalPath path2) {
	for (auto iter = path2.begin(); iter != path2.end(); ++iter) {
		path2.push_back(*iter);
	}
}

//Find all seeds as trivial paths
void FindSeeds(Graph& g, std::vector<BiderectionalPath>& seeds) {
	std::set<EdgeId> visited;
	std::map<EdgeId, BiderectionalPath> starts;

	for (auto iter = g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
		EdgeId e = *iter;

		//Try to make seed starting from unvisited edge
		if (visited.count(e) != 0) {
			visited.insert(e);

			BiderectionalPath newPath;
			newPath.push_back(e);
			starts[e] = newPath;

			//Extend trivially
			EdgeId nextStart = ExtendTrivialForward(g, newPath, starts);
			for (auto edgeInPath = newPath.begin(); edgeInPath != newPath.end(); ++edgeInPath) {
				visited.insert(*edgeInPath);
			}

			//If extended till another seed, than concatenate them
			if (nextStart != 0) {
				JoinPaths(newPath, starts[nextStart]);
				starts.erase(nextStart);
			}
		}
	}

	//Extending seed backward
	seeds.clear();
	seeds.reserve(starts.size());
	for (auto pathIter = starts.begin(); pathIter != starts.end(); ++pathIter) {
		ExtendTrivialBackward(g, pathIter->second);
		seeds.push_back(pathIter->second);
	}
}

//Recounting lengths form all edges to path's end
void RecountLengthsForward(Graph& g, BiderectionalPath& path, PathLengths& lengths) {
	lengths.clear();
	double currentLength = 0;

	for(auto iter = path.rbegin(); iter != path.rend(); ++iter) {
		currentLength += g.length(*iter);
		lengths.push_front(currentLength);
	}
}

//Recounting lengths form all edges to path's start
void ReocuntLengthsBackward(Graph& g, BiderectionalPath& path, PathLengths& lengths) {
	lengths.clear();
	double currentLength = 0;

	for(auto iter = path.begin(); iter != path.end(); ++iter) {
		lengths.push_back(currentLength);
		currentLength += g.length(*iter);
	}
}

size_t PathLength(Graph& g, BiderectionalPath& path) {
	double currentLength = 0;

	for(auto iter = path.begin(); iter != path.end(); ++iter) {
		currentLength += g.length(*iter);
	}
	return currentLength;
}

//Calculate weight for particular path extension
double ExtentionWeight(Graph& g, BiderectionalPath& path, PathLengths& lengths, EdgeId e, PairedInfoIndex<Graph> pairedInfo, bool forward) {
	double weight = 0;
	int edgeLength = forward ? 0 : g.length(e);

	for(size_t i = 0; i < path.size(); ++i) {
		omnigraph::PairedInfoIndex<Graph>::PairInfos pairs = pairedInfo.GetEdgePairInfo(path[i], e);
		int distance = lengths[i] + edgeLength;

		for (auto iter = pairs.begin(); iter != pairs.end(); ++iter) {
			int pairedDistance = rounded_d(*iter);
			//Can be modified according to distance comparison
			if (pairedDistance >= distance - DISTANCE_DIV &&
					pairedDistance <= distance + DISTANCE_DIV) {
				weight += iter->weight;
			}
		}
	}

	return weight;
}


void ExtendPathForward(Graph& g, BiderectionalPath& path, PathLengths& lengths, PairedInfoIndex<Graph> pairedInfo) {
//	std::vector<EdgeId> edges = g.OutgoingEdges(path.back());
//
//	for (auto iter = edges.begin(); iter != edges.end(); ++iter) {
//
//	}
}





}  // namespace long_contigs

#endif /* LONG_CONTIGS_HPP_ */

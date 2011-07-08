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
#include "debruijn/launch.hpp"
#include "omni/paired_info.hpp"

namespace long_contigs {

typedef std::deque<EdgeId> BiderectionalPath;
typedef std::deque<double> PathLengths;

//Extends trivial path forward
//If a start of another trivial path is found, returns it
//Otherwise returns 0
EdgeId ExtendTrivialForward(Graph& g, BiderectionalPath& path, std::map<EdgeId, BiderectionalPath&> starts) {
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
void FindSeeds(Graph& g, std::vector<BiderectionalPath&>& seeds) {
	std::set<EdgeId> visited;
	std::map<EdgeId, BiderectionalPath&> starts;

	for (auto iter = g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
		EdgeId e = *iter;

		if (visited.count(e) != 0) {
			visited.insert(e);

			BiderectionalPath newPath;
			newPath.push_back(e);
			starts[e] = newPath;

			EdgeId nextStart = ExtendTrivialForward(g, newPath, starts);
			for (auto edgeInPath = newPath.begin(); edgeInPath != newPath.end(); ++edgeInPath) {
				visited.insert(*edgeInPath);
			}

			if (nextStart != 0) {
				JoinPaths(newPath, starts[nextStart]);
				starts.erase(nextStart);
			}
		}
	}

	seeds.clear();
	seeds.reserve(starts.size());
	for (auto pathIter = starts.begin(); pathIter != starts.end(); ++pathIter) {
		ExtendTrivialBackward(g, pathIter->second);
		seeds.push_back(pathIter->second);
	}
}


void ReocuntLengthsForward(Graph& g, BiderectionalPath& path, PathLengths& lengths) {
	lengths.clear();
	double currentLength = 0;

	for(auto iter = path.begin(); iter != path.end(); ++iter) {
		currentLength += g.length(*iter);
		lengths.push_back(currentLength);
	}
}


void RecountLengthsBackward(Graph& g, BiderectionalPath& path, PathLengths& lengths) {
	lengths.clear();
	double currentLength = 0;

	for(auto iter = path.rbegin(); iter != path.rend(); ++iter) {
		lengths.push_front(currentLength);
		currentLength += g.length(*iter);
	}
}


double ExtentionWeight(Graph& g, BiderectionalPath& path, PathLengths& lengths, EdgeId e, PairedInfoIndex<Graph> pairedInfo) {
	double weight = 0;
	for(auto iter = path.begin(); iter != path.end(); ++iter) {
		PairInfos pairs = pairedInfo.GetEdgePairInfo(*iter, e);
		g.

	}
	return weight;
}


bool ExtendPath(Graph& g, BiderectionalPath& path, PairedInfoIndex<Graph> pairedInfo) {

}





}  // namespace long_contigs





#endif /* LONG_CONTIGS_HPP_ */

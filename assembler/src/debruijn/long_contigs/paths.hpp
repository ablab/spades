/*
 * paths.hpp
 *
 *  Created on: Aug 3, 2011
 *      Author: andrey
 */

#ifndef PATHS_HPP_
#define PATHS_HPP_

#include "lc_common.hpp"

namespace long_contigs {

using namespace debruijn_graph;

// ====== Extension functions ======

//Calculate weight for particular path extension
double ExtentionWeight(Graph& g, BidirectionalPath& path, PathLengths& lengths, EdgeId e, PairedInfoIndex<Graph>& pairedInfo, bool forward) {
	double weight = 0;
	int edgeLength = forward ? 0 : g.length(e);

	for(size_t i = 0; i < path.size(); ++i) {
		EdgeId edge = path[i];
		omnigraph::PairedInfoIndex<Graph>::PairInfos pairs = forward ? pairedInfo.GetEdgePairInfo(edge, e) : pairedInfo.GetEdgePairInfo(e, edge);
		int distance = lengths[i] + edgeLength;

		for (auto iter = pairs.begin(); iter != pairs.end(); ++iter) {
			int pairedDistance = rounded_d(*iter);
			//Can be modified according to distance comparison
			if (pairedDistance >= distance - DISTANCE_DEV &&
					pairedDistance <= distance + DISTANCE_DEV) {
				weight += iter->weight;
			}
		}
	}

	return weight / (double) std::min(g.length(e), READ_SIZE);
}

//Check whether selected extension is good enough
bool ExtensionGoodEnough(double weight) {
	//Condition of passing threshold is to be done
	return weight > WEIGHT_TRESHOLD;
}

//Choose best matching extension
//Threshold to be discussed
EdgeId ChooseExtension(Graph& g, BidirectionalPath& path, const std::vector<EdgeId>& edges,
		PathLengths& lengths, PairedInfoIndex<Graph>& pairedInfo, double& maxWeight, bool forward) {
	INFO("Choosing extension " << (forward ? "forward" : "backward"));

	maxWeight = 0;
	EdgeId bestEdge = 0;

	for (auto iter = edges.begin(); iter != edges.end(); ++iter) {
		INFO("Calculating weight");
		double weight = ExtentionWeight(g, path, lengths, *iter, pairedInfo, forward);
		INFO("Weight " << weight);

		if (weight > maxWeight) {
			maxWeight = weight;
			bestEdge = *iter;
		}
	}
	INFO("Best " << maxWeight);

	return ExtensionGoodEnough(maxWeight) ? bestEdge : 0;
}

//Increase path lengths
void IncreaseLengths(Graph& g, PathLengths& lengths, EdgeId edge, bool forward) {
	size_t len = g.length(edge);
	for(auto iter = lengths.begin(); iter != lengths.end(); ++iter) {
		*iter += len;
	}

	if (forward) {
		lengths.push_back(len);
	} else {
		lengths.push_front(0);
	}
}

//Check for cycling
bool CheckForCyclesSimple(Graph& g, BidirectionalPath& path, size_t datasetLen) {
	return PathLength(g, path) > datasetLen;
}

//Add edge to cycle detector and check
bool CheckCycle(BidirectionalPath& path, EdgeId extension, CycleDetector& detector, double weight) {
	detector.insert(std::make_pair(extension, std::make_pair(path.size(), weight)));

	return detector.count(extension) > MAX_LOOPS;
}

//Edges to remove
size_t CountLoopEdges(EdgeId lastEdge, CycleDetector& detector, bool fullRemoval) {
	auto iter = detector.upper_bound(lastEdge);

	--iter;
	size_t loopSize = iter->second.first;
	--iter;
	loopSize -= iter->second.first;

	if (fullRemoval) {
		return MAX_LOOPS * loopSize;
	} else {
		return (MAX_LOOPS - 1) * loopSize + 1;
	}
}

//Cut loop forward
void RemoveLoopForward(BidirectionalPath& path, CycleDetector& detector, bool fullRemoval) {
	size_t edgesToRemove = CountLoopEdges(path.back(), detector, fullRemoval);
	INFO("Removing loop of " << edgesToRemove << " edges");

	for(size_t i = 0; i < edgesToRemove; ++i) {
		path.pop_back();
	}
}

void RemoveLoopBackward(BidirectionalPath& path, CycleDetector& detector, bool fullRemoval) {
	size_t edgesToRemove = CountLoopEdges(path.front(), detector, fullRemoval);
	INFO("Removing loop of " << edgesToRemove << " edges");

	for(size_t i = 0; i < edgesToRemove; ++i) {
		path.pop_front();
	}
}


//Extend path forward
bool ExtendPathForward(Graph& g, BidirectionalPath& path, PathLengths& lengths,
		CycleDetector& detector, PairedInfoIndex<Graph>& pairedInfo) {

	double w;
	std::vector<EdgeId> edges = g.OutgoingEdges(g.EdgeEnd(path.back()));
	EdgeId extension = ChooseExtension(g, path, edges, lengths, pairedInfo, w, true);

	if (extension == 0) {
		return false;
	}
	path.push_back(extension);
	IncreaseLengths(g, lengths, extension, true);

	if (CheckCycle(path, extension, detector, w)) {
		INFO("Loop found");
		PrintPath(g, path);
		PrintDetector(detector);

		int aaa;
		std::cin >> aaa;

		RemoveLoopForward(path, detector, FULL_LOOP_REMOVAL);
		PrintPath(g, path);

		std::cin >> aaa;
		return false;
	}

	//PrintPath(g, path, lengths);

	return true;
}

//And backward
bool ExtendPathBackward(Graph& g, BidirectionalPath& path, PathLengths& lengths,
		CycleDetector& detector, PairedInfoIndex<Graph>& pairedInfo) {

	double w;
	std::vector<EdgeId> edges = g.IncomingEdges(g.EdgeStart(path.front()));
	EdgeId extension = ChooseExtension(g, path, edges, lengths, pairedInfo, w, false);

	if (extension == 0) {
		return false;
	}

	path.push_front(extension);
	IncreaseLengths(g, lengths, extension, false);

	if (CheckCycle(path, extension, detector, w)) {
		INFO("Loop found");
		PrintPath(g, path);
		PrintDetector(detector);

		//MakeKeyPause();

		RemoveLoopBackward(path, detector, FULL_LOOP_REMOVAL);
		PrintPath(g, path);

		//MakeKeyPause();
		return false;
	}

	//PrintPath(g, path, lengths);

	return true;
}

//Grow selected seed in both directions
void GrowSeed(Graph& g, BidirectionalPath& seed, PairedInfoIndex<Graph>& pairedInfo) {
	PathLengths lengths;
	CycleDetector detector;

	RecountLengthsForward(g, seed, lengths);

	//PrintPath(g, seed, lengths);

	INFO("Extending forward");
	while (ExtendPathForward(g, seed, lengths, detector, pairedInfo)) {
		INFO("Added edge");
	}

	detector.clear();
	RecountLengthsBackward(g, seed, lengths);

	//PrintPath(g, seed, lengths);

	INFO("Extending backward");
	while (ExtendPathBackward(g, seed, lengths, detector, pairedInfo)) {
		INFO("Added edge");
	}

	INFO("Growing done");
}

//Metrics for choosing seeds
size_t SeedPriority(const BidirectionalPath& seed) {
	return seed.size();
}

//Find paths with given seeds
void FindPaths(Graph& g, std::vector<BidirectionalPath>& seeds, PairedInfoIndex<Graph>& pairedInfo, std::vector<BidirectionalPath>& paths) {
	std::multimap<size_t, BidirectionalPath> priorityQueue;

	INFO("Finding paths started");
	for(auto seed = seeds.begin(); seed != seeds.end(); ++seed) {
		priorityQueue.insert(std::make_pair(SeedPriority(*seed), *seed));
	}

	for(auto seed = priorityQueue.rbegin(); seed != priorityQueue.rend(); ++seed) {
		GrowSeed(g, seed->second, pairedInfo);
		paths.push_back(seed->second);

		if (!ALL_SEEDS && PathsCoverage(g, paths) > EDGE_COVERAGE_TRESHOLD && PathsLengthCoverage(g, paths) > LENGTH_COVERAGE_TRESHOLD) {
			break;
		}
	}

	INFO("Finding paths finished");
}


} // namespace long_contigs


#endif /* PATHS_HPP_ */

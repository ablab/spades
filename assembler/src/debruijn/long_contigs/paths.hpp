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

//Recounting lengths form all edges to path's end
void RecountLengthsForward(Graph& g, BidirectionalPath& path, PathLengths& lengths) {
	lengths.clear();
	double currentLength = 0;

	for(auto iter = path.rbegin(); iter != path.rend(); ++iter) {
		currentLength += g.length(*iter);
		lengths.push_front(currentLength);
	}
}

//Recounting lengths from path's start to all edges
void RecountLengthsBackward(Graph& g, BidirectionalPath& path, PathLengths& lengths) {
	lengths.clear();
	double currentLength = 0;

	for(auto iter = path.begin(); iter != path.end(); ++iter) {
		lengths.push_back(currentLength);
		currentLength += g.length(*iter);
	}
}

// ====== Extension functions ======


//Calculate weight for particular path extension
double ExtentionWeight(Graph& g, BidirectionalPath& path, PathLengths& lengths, EdgeId e, PairedInfoIndexLibrary& pairedInfoLibrary, bool forward) {
	double weight = 0;
	int edgeLength = forward ? 0 : g.length(e);
	static int DISTANCE_DEV = CONFIG.read<bool>("etalon_info_mode") ? LC_CONFIG.read<int>("etalon_distance_dev") : LC_CONFIG.read<int>("real_distance_dev");

	for(size_t i = 0; i < path.size(); ++i) {
		EdgeId edge = path[i];
		omnigraph::PairedInfoIndex<Graph>::PairInfos pairs =
				forward ? pairedInfoLibrary.pairedInfoIndex->GetEdgePairInfo(edge, e) : pairedInfoLibrary.pairedInfoIndex->GetEdgePairInfo(e, edge);
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

	return weight / (double) std::min(g.length(e), pairedInfoLibrary.readSize);
}

//Check whether selected extension is good enough
bool ExtensionGoodEnough(double weight) {
	static size_t WEIGHT_TRESHOLD = LC_CONFIG.read<size_t>("weight_threshold");

	//Condition of passing threshold is to be done
	return weight > WEIGHT_TRESHOLD;
}

//Choose best matching extension
//Threshold to be discussed
EdgeId ChooseExtension(Graph& g, BidirectionalPath& path, const std::vector<EdgeId>& edges,
		PathLengths& lengths, PairedInfoIndices& pairedInfo, double& maxWeight, bool forward) {
	//INFO("Choosing extension " << (forward ? "forward" : "backward"));

	maxWeight = 0;
	EdgeId bestEdge = 0;

	for (auto iter = edges.begin(); iter != edges.end(); ++iter) {
		//INFO("Calculating weight");
		double weight = 0;
		for (auto lib = pairedInfo.begin(); lib != pairedInfo.end(); ++lib) {
			weight += ExtentionWeight(g, path, lengths, *iter, *lib, forward);
			//INFO(weight);
		}
		//INFO("Weight " << weight);

		if (weight > maxWeight) {
			maxWeight = weight;
			bestEdge = *iter;
		}
	}
	//INFO("Best " << maxWeight);

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
	static size_t MAX_LOOPS = LC_CONFIG.read<size_t>("max_loops");
	detector.insert(std::make_pair(extension, std::make_pair(path.size(), weight)));

	return detector.count(extension) > MAX_LOOPS;
}

//Edges to remove
size_t CountLoopEdges(EdgeId lastEdge, CycleDetector& detector, bool fullRemoval) {
	static size_t MAX_LOOPS = LC_CONFIG.read<size_t>("max_loops");
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
	//INFO("Removing loop of " << edgesToRemove << " edges");

	for(size_t i = 0; i < edgesToRemove; ++i) {
		path.pop_back();
	}
}

void RemoveLoopBackward(BidirectionalPath& path, CycleDetector& detector, bool fullRemoval) {
	size_t edgesToRemove = CountLoopEdges(path.front(), detector, fullRemoval);
	//INFO("Removing loop of " << edgesToRemove << " edges");

	for(size_t i = 0; i < edgesToRemove; ++i) {
		path.pop_front();
	}
}


//Extend path forward
bool ExtendPathForward(Graph& g, BidirectionalPath& path, PathLengths& lengths,
		CycleDetector& detector, PairedInfoIndices& pairedInfo) {

	double w;
	static bool FULL_LOOP_REMOVAL = LC_CONFIG.read<bool>("full_loop_removal");
	std::vector<EdgeId> edges = g.OutgoingEdges(g.EdgeEnd(path.back()));
	EdgeId extension = ChooseExtension(g, path, edges, lengths, pairedInfo, w, true);

	if (extension == 0) {
		return false;
	}
	path.push_back(extension);
	IncreaseLengths(g, lengths, extension, true);

	if (CheckCycle(path, extension, detector, w)) {
		//INFO("Loop found");
		//PrintPath(g, path);
		//PrintDetector(detector);
		//MakeKeyPause();

		RemoveLoopForward(path, detector, FULL_LOOP_REMOVAL);

		//PrintPath(g, path);
		//MakeKeyPause();
		return false;
	}

	//PrintPath(g, path, lengths);

	return true;
}

//And backward
bool ExtendPathBackward(Graph& g, BidirectionalPath& path, PathLengths& lengths,
		CycleDetector& detector, PairedInfoIndices& pairedInfo) {

	double w;
	static bool FULL_LOOP_REMOVAL = LC_CONFIG.read<bool>("full_loop_removal");
	std::vector<EdgeId> edges = g.IncomingEdges(g.EdgeStart(path.front()));
	EdgeId extension = ChooseExtension(g, path, edges, lengths, pairedInfo, w, false);

	if (extension == 0) {
		return false;
	}

	path.push_front(extension);
	IncreaseLengths(g, lengths, extension, false);

	if (CheckCycle(path, extension, detector, w)) {
		//INFO("Loop found");
		//PrintPath(g, path);
		//PrintDetector(detector);

		//MakeKeyPause();

		RemoveLoopBackward(path, detector, FULL_LOOP_REMOVAL);
		//PrintPath(g, path);

		//MakeKeyPause();
		return false;
	}

	//PrintPath(g, path, lengths);

	return true;
}

//Grow selected seed in both directions
void GrowSeed(Graph& g, BidirectionalPath& seed, PairedInfoIndices& pairedInfo) {
	PathLengths lengths;
	CycleDetector detector;

	RecountLengthsForward(g, seed, lengths);

	//PrintPath(g, seed, lengths);

	//INFO("Extending forward");
	while (ExtendPathForward(g, seed, lengths, detector, pairedInfo)) {
		//INFO("Added edge");
	}

	detector.clear();
	RecountLengthsBackward(g, seed, lengths);

	//PrintPath(g, seed, lengths);

	//INFO("Extending backward");
	while (ExtendPathBackward(g, seed, lengths, detector, pairedInfo)) {
		//INFO("Added edge");
	}

	//INFO("Growing done");
}

//Metrics for choosing seeds
size_t SeedPriority(const BidirectionalPath& seed) {
	return seed.size();
}

//Find paths with given seeds
void FindPaths(Graph& g, std::vector<BidirectionalPath>& seeds, PairedInfoIndices& pairedInfo, std::vector<BidirectionalPath>& paths) {
	std::multimap<size_t, BidirectionalPath> priorityQueue;
	static bool ALL_SEEDS = LC_CONFIG.read<bool>("all_seeds");
	static double EDGE_COVERAGE_TRESHOLD = LC_CONFIG.read<double>("edge_coverage");
	static double LENGTH_COVERAGE_TRESHOLD = LC_CONFIG.read<double>("len_coverage");

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


//Remove duplicate paths
void RemoveDuplicate(const std::vector<BidirectionalPath>& paths, std::vector<BidirectionalPath>& output) {
	for (auto path = paths.begin(); path != paths.end(); ++path) {
		bool copy = true;
		for (auto iter = output.begin(); iter != output.end(); ++iter) {
			if (ComparePaths(*path, *iter)) {
					copy = false;
					break;
			}
		}

		if (copy) {
			output.push_back(*path);
		}
	}
}


} // namespace long_contigs


#endif /* PATHS_HPP_ */

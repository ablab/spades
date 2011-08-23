/*
 * paths.hpp
 *
 *  Created on: Aug 3, 2011
 *      Author: andrey
 */

#ifndef PATHS_HPP_
#define PATHS_HPP_

#include "lc_common.hpp"
#include "loop.hpp"

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

// ====== Weight functions ======
//Weight filter
double WeightFunction(double weight) {
	return weight != 0 ? 1 : 0;
}

//Calculate weight
double GetWeight(omnigraph::PairedInfoIndex<Graph>::PairInfos pairs, int distance, int distanceDev, bool useWeightFunction = false) {
	double weight = 0;

	for (auto iter = pairs.begin(); iter != pairs.end(); ++iter) {
		int pairedDistance = rounded_d(*iter);

		if (lc_cfg::get().es.use_delta_first) {
			if (iter->variance != 0) {
				distanceDev = iter->variance;
			}
		}

		//Can be modified according to distance comparison
		if (pairedDistance >= distance - distanceDev &&
				pairedDistance <= distance + distanceDev) {
			weight += iter->weight;
		}
	}

	return useWeightFunction ? WeightFunction(weight) : weight;
}

//Weight fixing coefficient, sum weight of ideal info
int FixingCoefficient(Graph& g, const BidirectionalPath& path, EdgeId edge, size_t edgesToExclude, PairedInfoIndexLibrary& pairedInfoLibrary, bool forward) {
	int pathLen = 0;
	size_t start = forward ? 0 : edgesToExclude;
	size_t end = forward ? path.size() - edgesToExclude : path.size();
	for (size_t i = start; i < end; ++i) {
		pathLen += g.length(path[i]);
	}
	int exclLen = PathLength(g, path) - pathLen;
	int edgeLen = g.length(edge);

	int is = pairedInfoLibrary.insertSize;
	int rs = pairedInfoLibrary.readSize;

	int right = std::min(is, exclLen + edgeLen + rs);
	int left = std::max(exclLen - is, - rs - pathLen) + is;

	int delta = right - left + 1 - K;

	return delta > 0 ? delta : -1;
}

//Fixing weight value
double WeightFixing(Graph& g, const BidirectionalPath& path, EdgeId edge, size_t edgesToExclude, PairedInfoIndexLibrary& pairedInfoLibrary, double weight, bool forward) {
	//return weight / (double) std::min(g.length(edge), pairedInfoLibrary.readSize);
	return weight/(double) FixingCoefficient(g, path, edge, edgesToExclude, pairedInfoLibrary, forward);
}

//Calculate weight for particular path extension from one library
double ExtentionWeight(Graph& g, BidirectionalPath& path, PathLengths& lengths, EdgeId e, PairedInfoIndexLibrary& pairedInfoLibrary,
		size_t edgesToExclude, bool forward, bool useWeightFunction = false) {
	double weight = 0;
	int edgeLength = forward ? 0 : g.length(e);
	size_t start = forward ? 0 : edgesToExclude;
	size_t end = forward ? path.size() - edgesToExclude : path.size();

	static int DISTANCE_DEV = cfg::get().etalon_info_mode ? lc_cfg::get().es.etalon_distance_dev : pairedInfoLibrary.var;

	for(size_t i = start; i < end; ++i) {
		EdgeId edge = path[i];
		omnigraph::PairedInfoIndex<Graph>::PairInfos pairs =
				forward ? pairedInfoLibrary.pairedInfoIndex->GetEdgePairInfo(edge, e) : pairedInfoLibrary.pairedInfoIndex->GetEdgePairInfo(e, edge);
		int distance = lengths[i] + edgeLength;

		weight += GetWeight(pairs, distance, DISTANCE_DEV, useWeightFunction);
	}

	return WeightFixing(g, path, e, edgesToExclude, pairedInfoLibrary, weight, forward);
}

//Weight from a set of libraries
double ExtentionWeight(Graph& g, BidirectionalPath& path, PathLengths& lengths, EdgeId e, PairedInfoIndices& pairedInfo,
		size_t edgesToExclude, bool forward, bool useWeightFunction = false) {

	double weight = 0;
	for (auto lib = pairedInfo.begin(); lib != pairedInfo.end(); ++lib) {
		weight += ExtentionWeight(g, path, lengths, e, *lib, edgesToExclude, forward, useWeightFunction);
	}
	return weight;
}

// ====== Extension functions ======

//Check whether selected extension is good enough
EdgeId ExtensionGoodEnough(EdgeId edge, double weight, double threshold) {
	//Condition of passing threshold is to be done
	return weight > threshold ? edge : 0;
}

//Check whether selected extension is good enough
EdgeId ExtensionGoodEnough(EdgeId edge, double weight, double threshold, Graph& g, BidirectionalPath& path, PathStopHandler& handler, bool forward) {
	//Condition of passing threshold is to be done
	if (weight > threshold) {
		return edge;
	} else {
		handler.AddStop(path, WEAK_EXTENSION, forward);
		return 0;
	}
}


//Select only best extensions
double FilterExtentions(Graph& g, BidirectionalPath& path, std::vector<EdgeId>& edges,
		PathLengths& lengths, PairedInfoIndices& pairedInfo, size_t edgesToExclude, bool forward, bool useWeightFunction = false) {

	std::multimap<double, EdgeId> weights;

	for (auto iter = edges.begin(); iter != edges.end(); ++iter) {
		double weight = ExtentionWeight(g, path, lengths, *iter, pairedInfo, edgesToExclude, forward, useWeightFunction);

		weights.insert(std::make_pair(weight, *iter));
	}

	DETAILED_INFO("Choosing weights " << (forward ? "forward" : "backward"))
	for (auto iter = weights.begin(); iter != weights.end(); ++iter) {
		DETAILED_INFO(iter->second << " (" << g.length(iter->second) << ") = " << iter->first);
	}

	//Filling maximum edges
	edges.clear();
	auto bestEdge = weights.lower_bound((--weights.end())->first / lc_cfg::get().es.priority_coeff);
	for (auto maxEdge = bestEdge; maxEdge != weights.end(); ++maxEdge) {
		edges.push_back(maxEdge->second);
	}

	return bestEdge->first;
}

//Choose best matching extension
//Threshold to be discussed
EdgeId ChooseExtension(Graph& g, BidirectionalPath& path, std::vector<EdgeId>& edges,
		PathLengths& lengths, PairedInfoIndices& pairedInfo, double& maxWeight, size_t edgesToExclude, bool forward,
		PathStopHandler& handler) {

	//INFO("Choosing extension " << (forward ? "forward" : "backward"));
	if (edges.size() == 0) {
		handler.AddStop(path, NO_EXTENSION, forward);
		return 0;
	}
	if (edges.size() == 1) {
		return edges.back();
	}

	EdgeId toReturn = 0;
	if (lc_cfg::get().rs.research_mode && lc_cfg::get().rs.force_to_cycle) {
		for (auto edge = edges.begin(); edge != edges.end(); ++edge) {
			if (g.length(*edge) == lc_cfg::get().rs.cycle_priority_edge) {
				toReturn = *edge;
			}
		}
	}

	static bool useWeightFunctionFirst = lc_cfg::get().es.use_weight_function_first;

	if (useWeightFunctionFirst) {
		FilterExtentions(g, path, edges, lengths, pairedInfo, edgesToExclude, forward, true);

		if (edges.size() == 1) {
			static double weightFunThreshold = lc_cfg::get().es.weight_fun_threshold;
			maxWeight = ExtentionWeight(g, path, lengths, edges.back(), pairedInfo, edgesToExclude, forward);

			return toReturn == 0 ? ExtensionGoodEnough(edges.back(), maxWeight, weightFunThreshold, g, path, handler, forward) : toReturn;
		}
	}

	maxWeight = FilterExtentions(g, path, edges, lengths, pairedInfo, edgesToExclude, forward);

	if (edges.size() == 1) {
		static double weightThreshold = lc_cfg::get().es.weight_threshold;
		return toReturn == 0 ? ExtensionGoodEnough(edges.back(), maxWeight, weightThreshold, g, path, handler, forward) : toReturn;
	}
	else if (edges.size() > 1) {
		DETAILED_INFO("Cannot choose extension, no obvious maximum");
		handler.AddStop(path, MANY_GOOD_EXTENSIONS, forward);
	}

	return toReturn;
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
	static size_t MAX_LOOPS = lc_cfg::get().lr.max_loops;
	detector.insert(std::make_pair(extension, std::make_pair(path.size(), weight)));

	return detector.count(extension) > MAX_LOOPS;
}

//Edges to remove
size_t CountLoopEdges(EdgeId lastEdge, CycleDetector& detector, bool fullRemoval) {
	static size_t MAX_LOOPS = lc_cfg::get().lr.max_loops;
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

//Count edges to be excluded
size_t EdgesToExcludeForward(Graph& g, BidirectionalPath& path) {
	VertexId currentVertex = g.EdgeEnd(path.back());
	size_t toExclude = 0;

	while (g.CheckUniqueIncomingEdge(currentVertex)) {
		currentVertex = g.EdgeStart(g.GetUniqueIncomingEdge(currentVertex));
		++toExclude;
	}

	//INFO("Edges to exclude backward " << toExclude);
	return toExclude;
}

//Count edges to be excludeD
size_t EdgesToExcludeBackward(Graph& g, BidirectionalPath& path) {
	VertexId currentVertex = g.EdgeStart(path.front());
	size_t toExclude = 0;

	while (g.CheckUniqueOutgoingEdge(currentVertex)) {
		currentVertex = g.EdgeEnd(g.GetUniqueOutgoingEdge(currentVertex));
		++toExclude;
	}

	return toExclude;
}


//Extend path forward
bool ExtendPathForward(Graph& g, BidirectionalPath& path, PathLengths& lengths,
		CycleDetector& detector, PairedInfoIndices& pairedInfo,
		PathStopHandler& handler) {

	double w = 0;
	static bool FULL_LOOP_REMOVAL = lc_cfg::get().lr.full_loop_removal;
	std::vector<EdgeId> edges = g.OutgoingEdges(g.EdgeEnd(path.back()));
	EdgeId extension = ChooseExtension(g, path, edges, lengths, pairedInfo, w, EdgesToExcludeForward(g, path), true, handler);

	if (extension == 0) {
		return false;
	}
	DETAILED_INFO("Chosen forward " << extension << " (" << g.length(extension) << ")");

	path.push_back(extension);
	IncreaseLengths(g, lengths, extension, true);
	if (lc_cfg::get().rs.detailed_output) {
		PrintPath(g, path, lengths);
	}

	if (CheckCycle(path, extension, detector, w)) {
		DETAILED_INFO("Cycle detected");
		handler.AddStop(path, LOOP, true);
		RemoveLoopForward(path, detector, FULL_LOOP_REMOVAL);
		if (lc_cfg::get().rs.detailed_output) {
			PrintPath(g, path, lengths);
		}
		return false;
	}

	return true;
}

//And backward
bool ExtendPathBackward(Graph& g, BidirectionalPath& path, PathLengths& lengths,
		CycleDetector& detector, PairedInfoIndices& pairedInfo,
		PathStopHandler& handler) {

	double w = 0;
	static bool FULL_LOOP_REMOVAL = lc_cfg::get().lr.full_loop_removal;
	std::vector<EdgeId> edges = g.IncomingEdges(g.EdgeStart(path.front()));
	EdgeId extension = ChooseExtension(g, path, edges, lengths, pairedInfo, w, EdgesToExcludeBackward(g, path), false, handler);

	if (extension == 0) {
		return false;
	}
	DETAILED_INFO("Chosen backward " << extension << " (" << g.length(extension) << ")");

	path.push_front(extension);
	IncreaseLengths(g, lengths, extension, false);
	if (lc_cfg::get().rs.detailed_output) {
		PrintPath(g, path, lengths);
	}

	if (CheckCycle(path, extension, detector, w)) {
		DETAILED_INFO("Cycle detected");
		handler.AddStop(path, LOOP, false);
		RemoveLoopBackward(path, detector, FULL_LOOP_REMOVAL);
		if (lc_cfg::get().rs.detailed_output) {
			PrintPath(g, path, lengths);
		}
		return false;
	}

	return true;
}

size_t GetMaxInsertSize(PairedInfoIndices& pairedInfo) {
	size_t maxIS = 0;
	for(auto lib = pairedInfo.begin(); lib != pairedInfo.end(); ++lib) {
		if (maxIS < lib->insertSize) {
			maxIS = lib->insertSize;
		}
	}
	return maxIS;
}

//Grow selected seed in both directions
void GrowSeed(Graph& g, BidirectionalPath& seed, PairedInfoIndices& pairedInfo, PathStopHandler& handler) {
	PathLengths lengths;
	CycleDetector detector;

	static size_t maxIS = GetMaxInsertSize(pairedInfo);
	int i = 0;
	bool stop = false;

	while (i < lc_cfg::get().es.max_iter && !stop) {
		RecountLengthsForward(g, seed, lengths);
		DETAILED_INFO("Before forward");
		if (lc_cfg::get().rs.detailed_output) {
			PrintPath(g, seed, lengths);
		}

		while (ExtendPathForward(g, seed, lengths, detector, pairedInfo, handler)) {
		}
		detector.clear();

		if (PathLength(g, seed) > maxIS) {
			stop = true;
		}

		RecountLengthsBackward(g, seed, lengths);
		DETAILED_INFO("Before backward");
		if (lc_cfg::get().rs.detailed_output) {
			PrintPath(g, seed, lengths);
		}
		while (ExtendPathBackward(g, seed, lengths, detector, pairedInfo, handler)) {
		}
		detector.clear();

		++i;
	}
}

//Metrics for choosing seeds
size_t SeedPriority(const BidirectionalPath& seed) {
	return seed.size();
}

//Find paths with given seeds
void FindPaths(Graph& g, std::vector<BidirectionalPath>& seeds, PairedInfoIndices& pairedInfo, std::vector<BidirectionalPath>& paths,
		PathStopHandler& handler) {
	std::multimap<size_t, BidirectionalPath> priorityQueue;
	static bool ALL_SEEDS = lc_cfg::get().sc.all_seeds;
	static double EDGE_COVERAGE_TRESHOLD = lc_cfg::get().sc.edge_coverage;
	static double LENGTH_COVERAGE_TRESHOLD = lc_cfg::get().sc.len_coverage;

	INFO("Finding paths started");
	for(auto seed = seeds.begin(); seed != seeds.end(); ++seed) {
		priorityQueue.insert(std::make_pair(SeedPriority(*seed), *seed));
	}

	for(auto seed = priorityQueue.rbegin(); seed != priorityQueue.rend(); ++seed) {
		GrowSeed(g, seed->second, pairedInfo, handler);
		paths.push_back(seed->second);

		if (!ALL_SEEDS && PathsCoverage(g, paths) > EDGE_COVERAGE_TRESHOLD && PathsLengthCoverage(g, paths) > LENGTH_COVERAGE_TRESHOLD) {
			break;
		}
	}

	INFO("Finding paths finished");
}


} // namespace long_contigs


#endif /* PATHS_HPP_ */

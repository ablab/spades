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
		PathLengths& lengths, PairedInfoIndices& pairedInfo, size_t edgesToExclude, bool forward,
		LoopDetector& detector,
		bool useWeightFunction = false) {

	std::multimap<double, EdgeId> weights;

	for (auto iter = edges.begin(); iter != edges.end(); ++iter) {
		double weight = ExtentionWeight(g, path, lengths, *iter, pairedInfo, edgesToExclude, forward, useWeightFunction);
		weights.insert(std::make_pair(weight, *iter));
		detector.temp.AddAlternative(*iter, weight);
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
		PathLengths& lengths, PairedInfoIndices& pairedInfo, double * maxWeight, size_t edgesToExclude, bool forward,
		LoopDetector& detector,
		PathStopHandler& handler) {

	detector.temp.clear();

	if (edges.size() == 0) {
		handler.AddStop(path, NO_EXTENSION, forward);
		return 0;
	}
	if (edges.size() == 1) {
		detector.temp.AddAlternative(edges.back(), 1);
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
		FilterExtentions(g, path, edges, lengths, pairedInfo, edgesToExclude, forward, detector, true);

		if (edges.size() == 1) {
			static double weightFunThreshold = lc_cfg::get().es.weight_fun_threshold;
			*maxWeight = ExtentionWeight(g, path, lengths, edges.back(), pairedInfo, edgesToExclude, forward);

			return toReturn == 0 ? ExtensionGoodEnough(edges.back(), *maxWeight, weightFunThreshold, g, path, handler, forward) : toReturn;
		}
	}

	*maxWeight = FilterExtentions(g, path, edges, lengths, pairedInfo, edgesToExclude, forward, detector);

	static double weightThreshold = lc_cfg::get().es.weight_threshold;
	if (edges.size() == 1) {
		return toReturn == 0 ? ExtensionGoodEnough(edges.back(), *maxWeight, weightThreshold, g, path, handler, forward) : toReturn;
	}
	else if (edges.size() > 1) {
		DETAILED_INFO("Cannot choose extension, no obvious maximum");
		handler.AddStop(path, MANY_GOOD_EXTENSIONS, forward);
	}

	return toReturn;
}

//Count edges to be excluded
size_t EdgesToExcludeForward(Graph& g, BidirectionalPath& path) {
	static bool maxCycles = lc_cfg::get().ss.max_cycles;
	static LoopDetector detector;
	detector.clear();

	VertexId currentVertex = g.EdgeEnd(path.back());
	size_t toExclude = 0;

	while (g.CheckUniqueIncomingEdge(currentVertex)) {
		EdgeId e = g.GetUniqueIncomingEdge(currentVertex);
		currentVertex = g.EdgeStart(e);
		++toExclude;

		detector.temp.clear();
		detector.temp.AddAlternative(e);
		detector.AddNewEdge(e, toExclude);
		if (CheckCycle(path, e, detector, maxCycles)) {
			INFO("Cycled trivial path");
			return 0;
		}
	}

	return toExclude;
}

//Count edges to be excludeD
size_t EdgesToExcludeBackward(Graph& g, BidirectionalPath& path) {
	static bool maxCycles = lc_cfg::get().ss.max_cycles;
	static LoopDetector detector;
	detector.clear();

	VertexId currentVertex = g.EdgeStart(path.front());
	size_t toExclude = 0;

	while (g.CheckUniqueOutgoingEdge(currentVertex)) {
		EdgeId e = g.GetUniqueOutgoingEdge(currentVertex);
		currentVertex = g.EdgeEnd(e);
		++toExclude;

		detector.temp.clear();
		detector.temp.AddAlternative(e);
		detector.AddNewEdge(e, toExclude);
		if (CheckCycle(path, e, detector, maxCycles)) {
			INFO("Cycled trivial path");
			return 0;
		}
	}

	return toExclude;
}

EdgeId FindExitFromLoop(BidirectionalPath& path, LoopDetector& detector) {
	INFO("Finding exit");
	if (CountLoopExits(path, path.back(), detector) > 1) {
		INFO("Many exists, will resolve in normal way");
		return 0;
	}
	return FindFirstFork(path, path.back(), detector);
}

void ImitateFork(Graph& g, BidirectionalPath& path, PathLengths& lengths,
		LoopDetector& detector, PairedInfoIndices& pairedInfo,
		EdgeId loopEdge, EdgeId loopExit, bool forward) {

	IncreaseLengths(g, lengths, loopEdge, forward);

	size_t edgesToExclude = forward ? EdgesToExcludeForward(g, path) : EdgesToExcludeBackward(g, path);

	detector.temp.clear();
	detector.temp.AddAlternative(loopExit, ExtentionWeight(g, path, lengths, loopExit, pairedInfo, edgesToExclude, forward));

	double w = ExtentionWeight(g, path, lengths, loopEdge, pairedInfo, edgesToExclude, forward);
	detector.temp.AddAlternative(loopEdge, w);
	detector.AddNewEdge(loopEdge, path.size(), w);
}

void ReducePathTo(BidirectionalPath& path, size_t newSize, bool forward) {
	if (path.size() < newSize) {
		return;
	}

	if (forward) {
		while (path.size() != newSize) {
			path.pop_back();
		}
	} else {
		while (path.size() != newSize) {
			path.pop_front();
		}
	}
}

//Find best loop path
void ResolveLoopForward(Graph& g, BidirectionalPath& path, PathLengths& lengths,
		LoopDetector& detector, PairedInfoIndices& pairedInfo,
		EdgeId loopEdge) {

	INFO("Resolving loop");
	size_t originalSize = path.size();
	size_t loopLength = 0;
	if (loopEdge == 0) {
		//Loop already found
		loopEdge = FindExitFromLoop(path, detector);

		if (loopEdge == 0) {
			INFO("Not found");
			return;
		}
		loopLength = CountLoopLength(g, path, detector, true);
	} else {
		loopLength = g.length(loopEdge) + g.length(path.back());
	}

	if (loopLength > GetMaxInsertSize(pairedInfo)) {
		return;
	}

	INFO("Finding forward fork");
	EdgeId loopExit = GetForwardFork(g, loopEdge);
	INFO("Calculating max loop count")
	size_t maxCycles = 2 * GetMaxInsertSize(pairedInfo) / loopLength + 2;

	size_t i = 0;
	INFO("Imitating loop")
	do {
		INFO("Extending trivially")
		ExtendTrivialForward(g, path, detector, &lengths);
		detector.print();

		path.push_back(loopEdge);
		INFO("Imitating fork")
		ImitateFork(g, path, lengths, detector, pairedInfo, loopEdge, loopExit, true);
		detector.print();

		++i;
	} while (i <= maxCycles && LoopBecameStable(loopEdge, detector));

	INFO("Counting proper loop count")
	size_t properSize = GetMaxExitIteration(loopEdge, loopExit, detector);
	size_t firstToExit = GetFirstExitIteration(loopEdge, loopExit, detector);
	if (firstToExit == std::numeric_limits<size_t>::max()) {
		firstToExit = GetFirstExitIteration(loopEdge, loopExit, detector, 1);
	}

	if (firstToExit == properSize) {
		INFO("Resolved fine");
	} else {
		INFO("Proper resolved: " << properSize << ", usual resolved: " << firstToExit << ", loop size " << CountLoopEdges(path.back(), detector));
	}
	if (firstToExit == 0) {
		ReducePathTo(path, properSize - 1, true);
		path.push_back(loopExit);
	} else {
		ReducePathTo(path, originalSize, true);
	}

	lengths.clear();
	RecountLengthsForward(g, path, lengths);
}


//Extend path forward
bool ExtendPathForward(Graph& g, BidirectionalPath& path, PathLengths& lengths,
		LoopDetector& detector, PairedInfoIndices& pairedInfo,
		PathStopHandler& handler) {

	double w = 0;
	static bool FULL_LOOP_REMOVAL = lc_cfg::get().lr.full_loop_removal;
	static size_t MAX_LOOPS = lc_cfg::get().lr.max_loops;
	static size_t LOOPS_TO_IVESTIGATE = lc_cfg::get().lr.loop_to_investigate;

	EdgeId loopEdge = 0;
	if (lc_cfg::get().lr.investigation) {
		loopEdge = IsEdgeInShortLoopForward(g, path.back());
		if (loopEdge != 0) {
			DETAILED_INFO("Seed already near loop");
			if (!ResolveLoopForward(g, path, lengths, detector, pairedInfo, loopEdge)) {
				return false;
			}
			loopEdge = 0;
		}
	}

	std::vector<EdgeId> edges = g.OutgoingEdges(g.EdgeEnd(path.back()));
	EdgeId extension = ChooseExtension(g, path, edges, lengths, pairedInfo, &w, EdgesToExcludeForward(g, path), true, detector, handler);
	if (extension == 0) {
		return false;
	}

	detector.AddNewEdge(extension, path.size(), w);
	IncreaseLengths(g, lengths, extension, true);
	path.push_back(extension);

	DETAILED_INFO("Chosen forward " << extension << " (" << g.length(extension) << ")");
	DetailedPrintPath(g, path, lengths);

	EdgeId loopEdge = IsEdgeInShortLoopForward(g, extension);
	if (lc_cfg::get().lr.investigation && (loopEdge != 0 || CheckCycle(path, extension, detector, LOOPS_TO_IVESTIGATE))) {
		ResolveLoopForward(g, path, lengths, detector, pairedInfo, loopEdge);
	}

	if (CheckCycle(path, extension, detector, MAX_LOOPS)) {
		RemoveLoopForward(path, detector, FULL_LOOP_REMOVAL, MAX_LOOPS);

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
		LoopDetector& detector, PairedInfoIndices& pairedInfo,
		PathStopHandler& handler) {

	double w = 0;
	static bool FULL_LOOP_REMOVAL = lc_cfg::get().lr.full_loop_removal;
	static size_t MAX_LOOPS = lc_cfg::get().lr.max_loops;
	static size_t LOOPS_TO_IVESTIGATE = lc_cfg::get().lr.loop_to_investigate;

	EdgeId loopEdge = 0;
	if (lc_cfg::get().lr.investigation) {
		loopEdge = IsEdgeInShortLoopBackward(g, path.front());
		if (loopEdge != 0) {
			DETAILED_INFO("Seed already near loop");
			if (!ResolveLoopBackward(g, path, lengths, detector, pairedInfo, loopEdge)){
				return false;
			}
			loopEdge = 0;
		}
	}

	EdgeId extension = ChooseExtension(g, path, edges, lengths, pairedInfo, &w, EdgesToExcludeBackward(g, path), false, detector, handler);
	if (extension == 0) {
		return false;
	}

	detector.AddNewEdge(extension, path.size(), w);
	IncreaseLengths(g, lengths, extension, false);
	path.push_front(extension);

	DETAILED_INFO("Chosen backward " << extension << " (" << g.length(extension) << ")");
	DetailedPrintPath(g, path, lengths);

	if (lc_cfg::get().lr.investigation) {
		loopEdge = IsEdgeInShortLoopBackward(g, extension);
		if ((loopEdge != 0 || CheckCycle(path, extension, detector, LOOPS_TO_IVESTIGATE)) &&
				!ResolveLoopBackward(g, path, lengths, detector, pairedInfo, loopEdge)) {
			return false;
		}
	}

	if (CheckCycle(path, extension, detector, MAX_LOOPS)) {
		RemoveLoopBackward(path, detector, FULL_LOOP_REMOVAL, MAX_LOOPS);

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

//Grow selected seed in both directions
void GrowSeed(Graph& g, BidirectionalPath& seed, PairedInfoIndices& pairedInfo, PathStopHandler& handler) {
	PathLengths lengths;
	LoopDetector detector;

	static size_t maxIS = GetMaxInsertSize(pairedInfo);
	int i = 0;
	bool stop = false;

	while (i < lc_cfg::get().es.max_iter && !stop) {
		RecountLengthsForward(g, seed, lengths);
		DETAILED_INFO("Before forward");
		DetailedPrintPath(g, seed, lengths);

		detector.clear();
		while (ExtendPathForward(g, seed, lengths, detector, pairedInfo, handler)) {
		}

		if (PathLength(g, seed) > maxIS) {
			stop = true;
		}

		RecountLengthsBackward(g, seed, lengths);
		DETAILED_INFO("Before backward");
		DetailedPrintPath(g, seed, lengths);

		detector.clear();
		while (ExtendPathBackward(g, seed, lengths, detector, pairedInfo, handler)) {
		}

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
	std::multimap<size_t, BidirectionalPath*, std::greater<size_t> > priorityQueue;
	static bool ALL_SEEDS = lc_cfg::get().sc.all_seeds;
	static double EDGE_COVERAGE_TRESHOLD = lc_cfg::get().sc.edge_coverage;
	static double LENGTH_COVERAGE_TRESHOLD = lc_cfg::get().sc.len_coverage;

	INFO("Finding paths started");
	for(auto seed = seeds.begin(); seed != seeds.end(); ++seed) {
		priorityQueue.insert(std::make_pair(SeedPriority(*seed), &(*seed)));
	}

	for(auto seed = priorityQueue.rbegin(); seed != priorityQueue.rend(); ++seed) {
		GrowSeed(g, *(seed->second), pairedInfo, handler);
		paths.push_back(*(seed->second));

		if (!ALL_SEEDS && PathsCoverage(g, paths) > EDGE_COVERAGE_TRESHOLD && PathsLengthCoverage(g, paths) > LENGTH_COVERAGE_TRESHOLD) {
			break;
		}
	}

	INFO("Finding paths finished");
}


} // namespace long_contigs


#endif /* PATHS_HPP_ */

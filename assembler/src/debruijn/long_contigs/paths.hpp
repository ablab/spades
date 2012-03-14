/*
 * paths.hpp
 *
 *  Created on: Aug 3, 2011
 *      Author: andrey
 */

#ifndef PATHS_HPP_
#define PATHS_HPP_

#include "extend.hpp"
#include "xmath.h"

namespace long_contigs {

using namespace debruijn_graph;

EdgeId FindExitFromLoop(BidirectionalPath& path, LoopDetector& detector,
		bool forward) {
	EdgeId lastEdge = forward ? path.back() : path.front();
	if (CountLoopExits(path, lastEdge, detector, forward) > 1) {
		DETAILED_DEBUG("Many exists, will resolve in normal way");
		return EdgeId(0);
	}
	return FindFirstFork(path, lastEdge, detector, forward);
}

void ImitateFork(const Graph& g, BidirectionalPath& path, PathLengths& lengths,
		LoopDetector& detector, PairedInfoIndices& pairedInfo, EdgeId loopEdge,
		EdgeId loopExit, bool forward, int excludeCycle = -1) {
	DEBUG("Imitating fork " << g.length(loopEdge) << " " << g.length(loopExit));
	size_t edgesToExclude =
			forward ?
					EdgesToExcludeForward(g, path) :
					EdgesToExcludeBackward(g, path);
	if (excludeCycle != -1) {
		edgesToExclude = excludeCycle;
	}

	detector.temp.clear();
	detector.temp.AddAlternative(
			loopExit,
			ExtentionWeight(g, path, lengths, loopExit, pairedInfo,
					edgesToExclude, forward));

	double w = ExtentionWeight(g, path, lengths, loopEdge, pairedInfo,
			edgesToExclude, forward);
	detector.temp.AddAlternative(loopEdge, w);
	detector.AddNewEdge(loopEdge, path.size(), w);
}

void RemoveFromDetector(LoopDetector& detector, EdgeId e, size_t iteration,
		EdgeId substitution = EdgeId(0)) {
	auto upper = --detector.data.upper_bound(e);
	auto lower = detector.data.upper_bound(e);

	for (auto iter = upper; iter != lower; --iter) {
		if (iter->second.iteration == iteration) {
			if (substitution == EdgeId(0)) {
				detector.data.erase(iter);
			} else {
				auto newData = std::make_pair(substitution, iter->second);
				detector.data.erase(iter);
				detector.data.insert(newData);
			}
			return;
		}
	}

	if (lower->second.iteration == iteration) {
		if (substitution == EdgeId(0)) {
			detector.data.erase(lower);
		} else {
			auto newData = std::make_pair(substitution, lower->second);
			detector.data.erase(lower);
			detector.data.insert(newData);
		}
		return;
	}

	DETAILED_DEBUG("Iteration not found in detector!");
}

void ReducePathTo(BidirectionalPath& path, LoopDetector& detector,
		size_t newSize, EdgeId loopExit, bool forward) {
	if (path.size() < newSize) {
		return;
	}

	if (forward) {
		while (path.size() > newSize + 1) {
			RemoveFromDetector(detector, path.back(), path.size() - 1);
			path.pop_back();
		}
		RemoveFromDetector(detector, path.back(), path.size() - 1, loopExit);
		path.pop_back();
	} else {
		while (path.size() > newSize + 1) {
			RemoveFromDetector(detector, path.front(), path.size() - 1);
			path.pop_front();
		}
		RemoveFromDetector(detector, path.front(), path.size() - 1, loopExit);
		path.pop_front();
	}
}

//Find length
bool CheckLoop(const Graph& g, BidirectionalPath& path, LoopDetector& detector,
		EdgeId& loopEdge, size_t& loopLength, bool forward, size_t& loopSize) {
	loopLength = 0;
	if (loopEdge == EdgeId(0)) {
		//Loop already found
		DETAILED_DEBUG("Loop was already found");
		loopLength = CountLoopLength(g, path, detector, forward);

		if (PathIsOnlyLoop(path, detector, forward)) {
			DETAILED_DEBUG("Not enough info");
			return false;
		}

		loopEdge = FindExitFromLoop(path, detector, forward);
		if (loopEdge == EdgeId(0)) {
			DETAILED_DEBUG("Not found");
			return false;
		}
		loopSize = CountLoopEdges(forward ? path.back() : path.front(),
				detector);
	} else {
		loopLength = g.length(loopEdge)
				+ g.length(forward ? path.back() : path.front());

		DETAILED_DEBUG("Short loop");
		if (PathIsOnlyLoop(path, loopEdge, forward)) {
			DETAILED_DEBUG("Not enough info");
			return false;
		}
		loopSize = 2;
	}
	return true;
}

bool MakeCorrectLoop(BidirectionalPath& path, LoopDetector& detector,
		EdgeId loopEdge, EdgeId loopExit, size_t originalSize, bool forward) {
	size_t properSize = GetMaxExitIteration(loopEdge, loopExit, detector,
			std::make_pair(originalSize - 1, path.size() - 1));
	size_t firstToExit = GetFirstExitIteration(loopEdge, loopExit, detector,
			std::make_pair(originalSize - 1, path.size() - 1));
//	if (firstToExit == std::numeric_limits<size_t>::max()) {
//		firstToExit = GetFirstExitIteration(loopEdge, loopExit, detector, 1);
//	}

	size_t loopLen = CountLoopEdges(forward ? path.back() : path.front(),
			detector);
	if (firstToExit == properSize) {
		DETAILED_DEBUG(
				"Resolved fine: " << properSize << ", usual resolved: " << firstToExit << ", loop size " << loopLen);
	} else {
		DEBUG(
				"Proper resolved: " << properSize << ", usual resolved: " << firstToExit << ", loop size " << loopLen);
	}

	if (properSize != 0) {
		ReducePathTo(path, detector, properSize, loopExit, forward);
	} else if (firstToExit != std::numeric_limits<size_t>::max()) {
		ReducePathTo(path, detector, firstToExit, loopExit, forward);
	} else {
		DEBUG("Cannot detect proper cycle exit!");
		ReducePathTo(path, detector, originalSize, EdgeId(0), forward);
		return false;
	}

	if (forward) {
		path.push_back(loopExit);
	} else {
		path.push_front(loopExit);
	}
	return true;
}

//Find best loop path
bool ResolveLoopForward(const Graph& g, BidirectionalPath& path, PathLengths& lengths,
		LoopDetector& detector, PairedInfoIndices& pairedInfo,
		EdgeId loopEdge) {

	DETAILED_DEBUG("Resolving loop forward");
	size_t originalSize = path.size();
	size_t loopLength = 0;
	size_t loopSize = 0;
	bool goodLoop = CheckLoop(g, path, detector, loopEdge, loopLength, true,
			loopSize);

	if (loopLength > GetMaxInsertSize(pairedInfo) - K) {
		DETAILED_DEBUG("Loop is too long");
		return !params.ps.lr.stop_on_long;
	}
	if (!goodLoop) {
		return true;
	}

	EdgeId loopExit = GetForwardFork(g, loopEdge);
	if (loopExit == EdgeId(0)) {
		return true;
	}
	size_t maxCycles = 2 * GetMaxInsertSize(pairedInfo) / loopLength + 2;

	size_t i = 0;
	do {
		ExtendTrivialForward(g, path, detector, &lengths);

		int excludeCycle =
				(params.ps.lr.exlude_cycle && loopSize == 2) ?
						loopSize * i + 1 : -1;
		ImitateFork(g, path, lengths, detector, pairedInfo, loopEdge, loopExit,
				true, excludeCycle);

		path.push_back(loopEdge);
		IncreaseLengths(g, lengths, loopEdge, true);

		++i;
	} while (i <= maxCycles && !LoopBecameStable(loopEdge, detector));

	detector.print(g);

	bool result = MakeCorrectLoop(path, detector, loopEdge, loopExit,
			originalSize, true);
	lengths.clear();
	RecountLengthsForward(g, path, lengths);
	DETAILED_DEBUG("Resolved");

	return result;
}

//Find best loop path
bool ResolveLoopBackward(const Graph& g, BidirectionalPath& path,
		PathLengths& lengths, LoopDetector& detector,
		PairedInfoIndices& pairedInfo, EdgeId loopEdge) {

	DETAILED_DEBUG("Resolving loop backward");
	size_t originalSize = path.size();
	size_t loopLength = 0;
	size_t loopSize = 0;
	bool goodLoop = CheckLoop(g, path, detector, loopEdge, loopLength, false,
			loopSize);

	if (loopLength > GetMaxInsertSize(pairedInfo) - K) {
		DETAILED_DEBUG("Loop is too long");
		return !params.ps.lr.stop_on_long;
	}
	if (!goodLoop) {
		return true;
	}

	EdgeId loopExit = GetBackwardFork(g, loopEdge);
	if (loopExit == EdgeId(0)) {
		return true;
	}
	size_t maxCycles = 2 * GetMaxInsertSize(pairedInfo) / loopLength + 2;

	size_t i = 0;
	DETAILED_DEBUG("Imitating loop backward")
	do {
		DEBUG("Extending trivially backward")
		ExtendTrivialBackward(g, path, detector, &lengths);
		int excludeCycle =
				(params.ps.lr.exlude_cycle && loopSize == 2) ?
						loopSize * i + 1 : -1;
		ImitateFork(g, path, lengths, detector, pairedInfo, loopEdge, loopExit,
				false, excludeCycle);

		path.push_front(loopEdge);
		IncreaseLengths(g, lengths, loopEdge, false);

		DetailedPrintPath(g, path, lengths);

		++i;
	} while (i <= maxCycles && !LoopBecameStable(loopEdge, detector));

	detector.print(g);

	bool result = MakeCorrectLoop(path, detector, loopEdge, loopExit,
			originalSize, false);
	lengths.clear();
	RecountLengthsBackward(g, path, lengths);
	DETAILED_DEBUG("Resolved");

	return result;
}

//Extend path forward
bool ExtendPathForward(const Graph& g, BidirectionalPath& path, PathLengths& lengths,
		LoopDetector& detector, PairedInfoIndices& pairedInfo,
		PathStopHandler& handler, JumpingHero<Graph>& hero) {

	DETAILED_DEBUG("Try to forward-extend path " << path.uid << ": " << g.str(path));

	if (path.empty()) {
		WARN("Extension of path " << path.uid << " ended unexpectedly");
		return false;
	}

	double w = 0;
	static bool FULL_LOOP_REMOVAL = params.ps.lr.full_loop_removal;
	static size_t MAX_LOOPS = params.ps.lr.max_loops;
	static size_t LOOPS_TO_IVESTIGATE = params.ps.lr.loop_to_investigate;

	EdgeId loopEdge(0);
	if (params.ps.lr.investigation) {
		loopEdge = IsEdgeInShortLoopForward(g, path.back());
		if (loopEdge != EdgeId(0)
				|| CheckCycle(path, path.back(), detector, LOOPS_TO_IVESTIGATE)) {
			DETAILED_DEBUG("Seed already near loop");
			if (!ResolveLoopForward(g, path, lengths, detector, pairedInfo,
					loopEdge)) {
				handler.AddStop(&path, LONG_LOOP, true);
				WARN("Extension of path " << path.uid << " ended unexpectedly");
				return false;
			}
			loopEdge = EdgeId(0);
		}
	}

	std::vector<EdgeId> edges = g.OutgoingEdges(g.EdgeEnd(path.back()));

//	cout << "Choosing extension among edges " << g.str(edges) << endl;
	DETAILED_DEBUG("Choosing extension after edge " << g.str(path.back()));

	EdgeId extension = ChooseExtension(g, path, edges, lengths, pairedInfo, &w,
			EdgesToExcludeForward(g, path), true, detector, handler, hero);
	if (extension == EdgeId(0)) {
		DETAILED_DEBUG("No good extension found after edge " << g.str(path.back()));
		return false;
	}

	DETAILED_DEBUG(
			"Chosen forward extension " << g.str(extension) << " after edge " << g.str(path.front()));

	detector.AddNewEdge(extension, path.size(), w);
	IncreaseLengths(g, lengths, extension, true);
	hero.ProcessEdge(extension);
	path.push_back(extension);

	DetailedPrintPath(g, path, lengths);

	if (params.ps.lr.investigation) {
		loopEdge = IsEdgeInShortLoopForward(g, extension);
		if (loopEdge != EdgeId(0)
				|| CheckCycle(path, extension, detector, LOOPS_TO_IVESTIGATE)) {
			if (!ResolveLoopForward(g, path, lengths, detector, pairedInfo,
					loopEdge)) {
				handler.AddStop(&path, LONG_LOOP, true);
				WARN("Extension of path " << path.uid << " ended unexpectedly");
				return false;
			}
		}
	}

	if (CheckCycle(path, extension, detector, MAX_LOOPS)) {
		detector.print(g);
		RemoveLoopForward(path, detector, FULL_LOOP_REMOVAL, MAX_LOOPS);

		DETAILED_DEBUG("Cycle detected");
		DetailedPrintPath(g, path, lengths);
		handler.AddStop(&path, LOOP, true);
		WARN("Extension of path " << path.uid << " ended unexpectedly");
		return false;
	}

	return true;
}

//And backward
bool ExtendPathBackward(const Graph& g, BidirectionalPath& path, PathLengths& lengths,
		LoopDetector& detector, PairedInfoIndices& pairedInfo,
		PathStopHandler& handler, JumpingHero<Graph>& hero) {

	DETAILED_DEBUG("Try to backward-extend path " << path.uid << ": " << g.str(path));

	if (path.empty()) {
		WARN("Extension of path " << path.uid << " ended unexpectedly");
		return false;
	}

	double w = 0;
	static bool FULL_LOOP_REMOVAL = params.ps.lr.full_loop_removal;
	static size_t MAX_LOOPS = params.ps.lr.max_loops;
	static size_t LOOPS_TO_IVESTIGATE = params.ps.lr.loop_to_investigate;

	EdgeId loopEdge(0);
	if (params.ps.lr.investigation) {
		loopEdge = IsEdgeInShortLoopBackward(g, path.front());
		if (loopEdge != EdgeId(0)
				|| CheckCycle(path, path.front(), detector, LOOPS_TO_IVESTIGATE)) {
			DETAILED_DEBUG("Seed already near loop");
			if (!ResolveLoopBackward(g, path, lengths, detector, pairedInfo,
					loopEdge)) {
				handler.AddStop(&path, LONG_LOOP, false);
				return false;
			}
			loopEdge = EdgeId(0);
		}
	}

	std::vector<EdgeId> edges = g.IncomingEdges(g.EdgeStart(path.front()));

	DETAILED_DEBUG("Choosing extension before edge " << g.str(path.front()));

	EdgeId extension = ChooseExtension(g, path, edges, lengths, pairedInfo, &w,
			EdgesToExcludeBackward(g, path), false, detector, handler, hero);
	if (extension == EdgeId(0)) {
		DETAILED_DEBUG("No good extension found before edge " << g.str(path.front()));
		return false;
	}

	DETAILED_DEBUG(
			"Chosen backward extension " << g.str(extension) << " before edge " << g.str(path.front()));

	detector.AddNewEdge(extension, path.size(), w);
	IncreaseLengths(g, lengths, extension, false);
	hero.ProcessEdge(extension);
	path.push_front(extension);

	DetailedPrintPath(g, path, lengths);

	if (params.ps.lr.investigation) {
		loopEdge = IsEdgeInShortLoopBackward(g, extension);
		if (loopEdge != EdgeId(0)
				|| CheckCycle(path, extension, detector, LOOPS_TO_IVESTIGATE)) {
			if (!ResolveLoopBackward(g, path, lengths, detector, pairedInfo,
					loopEdge)) {
				handler.AddStop(&path, LONG_LOOP, false);
				WARN("Extension of path " << path.uid << " ended unexpectedly");
				return false;
			}
		}
	}

	if (CheckCycle(path, extension, detector, MAX_LOOPS)) {
		detector.print(g);
		RemoveLoopBackward(path, detector, FULL_LOOP_REMOVAL, MAX_LOOPS);

		DETAILED_DEBUG("Cycle detected");
		DetailedPrintPath(g, path, lengths);
		handler.AddStop(&path, LOOP, false);
		WARN("Extension of path " << path.uid << " ended unexpectedly");
		return false;
	}

	return true;
}

template<class It>
boost::optional<EdgeId> FindFirstLongEdge(const Graph& g, double length_bound,
		It begin, It end) {
	for (It it = begin; it != end; ++it) {
		if (g.length(*it) > length_bound) {
			return boost::optional<EdgeId>(*it);
		}
	}
	return boost::none;
}

template<class It>
boost::optional<EdgeId> FindLastLongEdge(const Graph& g, double length_bound,
		It begin, It end) {
	boost::optional<EdgeId> answer;
	for (It it = begin; it != end; ++it) {
		if (g.length(*it) > length_bound) {
			answer.reset(*it);
		}
	}
	return answer;
}


//Grow selected seed in both directions
void GrowSeed(const Graph& g, BidirectionalPath& seed, PairedInfoIndices& pairedInfo,
		PathStopHandler& handler,
		const PairedInfoIndex<Graph>& jump_index) {
	PathLengths lengths;
	LoopDetector detector(g);

	static size_t maxIS = GetMaxInsertSize(pairedInfo);
	int i = 0;
	bool stop = false;

	while (i < params.ps.es.max_iter && !stop) {
		RecountLengthsForward(g, seed, lengths);

		DETAILED_DEBUG("Before forward");
		DetailedPrintPath(g, seed, lengths);

		RecountDetectorForward(g, seed, pairedInfo, detector);

		VERIFY(cfg::get().jump.weight_threshold > 0. && cfg::get().jump.weight_threshold < 10.);

		JumpingHero<Graph> forward_hero(g, seed, jump_index, 2000, /*invalidation length*/3500
				, true, cfg::get().jump.weight_threshold);

		while (ExtendPathForward(g, seed, lengths, detector, pairedInfo,
				handler, forward_hero)) {
		}

		if (PathLength(g, seed) > maxIS) {
			stop = true;
		}

		RecountLengthsBackward(g, seed, lengths);
		DETAILED_DEBUG("Before backward");
		DetailedPrintPath(g, seed, lengths);

		RecountDetectorBackward(g, seed, pairedInfo, detector);

		JumpingHero<Graph> backward_hero(g, seed, jump_index, 2000, 3500
				, false, cfg::get().jump.weight_threshold);

		while (ExtendPathBackward(g, seed, lengths, detector, pairedInfo,
				handler, backward_hero)) {
		}

		++i;
	}
}

//Metrics for choosing seeds
size_t SeedPriority(const BidirectionalPath& seed) {
	return seed.size();
}

//Find paths with given seeds
void FindPaths(const Graph& g, std::vector<BidirectionalPath>& seeds,
		PairedInfoIndices& pairedInfo, PathStopHandler& handler,
		const PairedInfoIndex<Graph>& jump_index) {

//	static bool ALL_SEEDS = params.ps.sc.all_seeds;
//	static double EDGE_COVERAGE_TRESHOLD = params.ps.sc.edge_coverage;
//	static double LENGTH_COVERAGE_TRESHOLD = params.ps.sc.len_coverage;

	INFO("Finding paths started");
	for (auto seed = seeds.begin(); seed != seeds.end(); ++seed) {
		GrowSeed(g, *seed, pairedInfo, handler, jump_index);
		DETAILED_DEBUG("Growing seed w/length " << PathLength(g, *seed));

//		if (!ALL_SEEDS && PathsCoverage(g, paths) > EDGE_COVERAGE_TRESHOLD && PathsLengthCoverage(g, paths) > LENGTH_COVERAGE_TRESHOLD) {
//			break;
//		}
	}

	INFO("Finding paths finished");
}

// === Totally symmertic mode ===

void CompareConjugateGrowth(const Graph& g, BidirectionalPath& path,
		PathLengths& lengths, LoopDetector& detector,
		BidirectionalPath& conjPath, PathLengths& conjLengths,
		LoopDetector& conjDetector) {

	DEBUG("After symmetric growth")
	if (!ComparePaths(path, conjPath)) {
		PrintPath(g, path, lengths);
		detector.print(g);
		PrintPath(g, conjPath, conjLengths);
		conjDetector.print(g);

	} else {
		DEBUG("Paths are equal!")
	}

}

void SymmetrizePaths() {

}

//todo ask andrew
//Grow selected seed in both directions
//void GrowSeedSymmetric(const Graph& g, BidirectionalPath& seed,
//		BidirectionalPath& conjSeed, PairedInfoIndices& pairedInfo,
//		PathStopHandler& handler) {
//
//	PathLengths lengths;
//	LoopDetector detector;
//	PathLengths conjLengths;
//	LoopDetector conjDetector;
//
//	static size_t maxIS = GetMaxInsertSize(pairedInfo);
//	int i = 0;
//	bool stop = false;
//	bool start = true;
//
//	while (i < params.ps.es.max_iter && !stop) {
//		if (!start || params.first_grow_forward) {
//			RecountLengthsForward(g, seed, lengths);
//			RecountLengthsBackward(g, conjSeed, conjLengths);
//
//			DETAILED_DEBUG("Before forward");
//			DetailedPrintPath(g, seed, lengths);
//			DetailedPrintPath(g, conjSeed, conjLengths);
//
//			RecountDetectorForward(g, seed, pairedInfo, detector);
//			RecountDetectorBackward(g, conjSeed, pairedInfo, conjDetector);
//
//			while (ExtendPathForward(g, seed, lengths, detector, pairedInfo,
//					handler)) {
//			}
//			while (ExtendPathBackward(g, conjSeed, conjLengths, conjDetector,
//					pairedInfo, handler)) {
//			}
//
//			if (PathLength(g, seed) > maxIS) {
//				stop = true;
//			}
//
//		}
//		start = false;
//
//		RecountLengthsBackward(g, seed, lengths);
//		RecountLengthsForward(g, conjSeed, conjLengths);
//
//		DETAILED_DEBUG("Before backward");
//		DetailedPrintPath(g, seed, lengths);
//		DetailedPrintPath(g, conjSeed, conjLengths);
//
//		RecountDetectorBackward(g, seed, pairedInfo, detector);
//		RecountDetectorForward(g, conjSeed, pairedInfo, conjDetector);
//
//		while (ExtendPathBackward(g, seed, lengths, detector, pairedInfo,
//				handler)) {
//		}
//		while (ExtendPathForward(g, conjSeed, conjLengths, conjDetector,
//				pairedInfo, handler)) {
//		}
//
//		++i;
//	}
//}

//todo ask andrew
//void FindPathsSymmetric(const Graph& g, std::vector<BidirectionalPath>& seeds,
//		PairedInfoIndices& pairedInfo, PathStopHandler& handler,
//		std::vector<int>& seedPairs) {
//
//	DEBUG("Finding paths started in totally symmetric way");
//	std::set<int> grown;
//	for (int i = 0; i < (int) seeds.size(); ++i) {
//		if (grown.count(i) == 0) {
//			GrowSeedSymmetric(g, seeds[i], seeds[seedPairs[i]], pairedInfo,
//					handler);
//			grown.insert(i);
//			grown.insert(seedPairs[i]);
//		}
//	}
//
//	DEBUG("Finding paths finished");
//}

} // namespace long_contigs

#endif /* PATHS_HPP_ */

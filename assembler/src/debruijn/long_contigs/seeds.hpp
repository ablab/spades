/*
 * seeds.hpp
 *
 *  Created on: Aug 3, 2011
 *      Author: andrey
 */

#ifndef SEEDS_HPP_
#define SEEDS_HPP_

#include "lc_common.hpp"
#include "lc_io.hpp"
#include "path_utils.hpp"
#include "loop.hpp"
#include "extend.hpp"


namespace long_contigs {

using namespace debruijn_graph;

// ====== Seed functions ======

//Extends trivial path forward
//If a start of another trivial path is found, returns it
//Otherwise returns 0
void ExtendTrivialForward(const Graph& g, BidirectionalPath& path, LoopDetector& detector,
		PathLengths* lengths = 0, PairedInfoIndices * pairedInfo = 0) {
	static bool maxCycles = params.ps.ss.max_cycles;

	if (path.empty()) {
		return;
	}

	VertexId currentVertex = g.EdgeEnd(path.back());
	while (g.CheckUniqueOutgoingEdge(currentVertex)) {
		EdgeId nextEdge = g.GetUniqueOutgoingEdge(currentVertex);

		if (pairedInfo != 0 &&  params.ps.ss.check_trusted) {

			int toExclude = 0;
			double weight =
					ExtentionWeight(g, path, *lengths, nextEdge, *pairedInfo, toExclude, true, false);

			DETAILED_DEBUG("Forward " << nextEdge << " (" << g.length(nextEdge) << "), weight " << weight);
			DetailedPrintPath(g, path, *lengths);

			if (ExtensionGoodEnough(nextEdge, weight, params.ps.ss.trusted_threshold) == EdgeId(0)) {
				break;
			}
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
			if (path.size() >= 2 && path.back() == path[path.size() - 2]) {
				path.pop_back();
				if (lengths != 0) {
					RecountLengthsForward(g, path, *lengths);
				}
				if (pairedInfo != 0) {
					RecountDetectorForward(g, path, *pairedInfo, detector);
				}
			}
			break;
		}
	}
	return;
}


//Trivially extend path backward
void ExtendTrivialBackward(const Graph& g, BidirectionalPath& path, LoopDetector& detector, PathLengths* lengths = 0, PairedInfoIndices * pairedInfo = 0) {
	static bool maxCycles = params.ps.ss.max_cycles;

	if (path.empty()) {
		return;
	}

	VertexId currentVertex = g.EdgeStart(path.front());
	while (g.CheckUniqueIncomingEdge(currentVertex)) {
		EdgeId nextEdge = g.GetUniqueIncomingEdge(currentVertex);

		if (pairedInfo != 0 &&  params.ps.ss.check_trusted) {

			int toExclude = 0;
			double weight =
					ExtentionWeight(g, path, *lengths, nextEdge, *pairedInfo, toExclude, false, false);

			DETAILED_DEBUG("Backward " << nextEdge << " (" << g.length(nextEdge) << "), weight " << weight);
			DetailedPrintPath(g, path, *lengths);

			if (ExtensionGoodEnough(nextEdge, weight, params.ps.ss.trusted_threshold) == EdgeId(0)) {
				break;
			}
		}

		path.push_front(nextEdge);
		if (lengths != 0) {
			IncreaseLengths(g, *lengths, nextEdge, false);
		}
		currentVertex = g.EdgeStart(nextEdge);

		detector.temp.clear();
		detector.temp.AddAlternative(nextEdge);
		detector.AddNewEdge(nextEdge, path.size() - 1);
		if (CheckCycle(path, nextEdge, detector, maxCycles)) {
			if (path.size() >= 2 && path.front() == path[1]) {
				path.pop_front();
				if (lengths != 0) {
					RecountLengthsBackward(g, path, *lengths);
				}
				if (pairedInfo != 0) {
					RecountDetectorBackward(g, path, *pairedInfo, detector);
				}
			}
			break;
		}
	}
}

//Glue second path to the first one
void JoinPaths(BidirectionalPath& path1, BidirectionalPath& path2) {
	if (path1 == path2) {
		DEBUG("Cannot join path with itself");
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
void FindSeeds(const Graph& g, std::vector<BidirectionalPath>& seeds, PairedInfoIndices * pairedInfo = 0) {
	LoopDetector detector(g);
	PathLengths lengths;
	set<EdgeId> edges;
	seeds.clear();

	INFO("Finding seeds started");
	for (auto iter = g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
		EdgeId e = *iter;

		if ((g.length(e) >= params.ps.ss.chimeric_len - params.ps.ss.chimeric_delta && g.length(e) <= params.ps.ss.chimeric_len + params.ps.ss.chimeric_delta) ||
				(g.length(e) <= params.ps.ss.short_single &&
						(ClassifyEdge(g, e) == TOTALY_ISOLATED || ClassifyEdge(g, e) == HAS_NEIGHBOUR)) ||
						g.coverage(e) < params.ps.ss.start_egde_coverage) {
			edges.insert(e);
			edges.insert(g.conjugate(e));
		}
		else if (edges.count(e) == 0) {
			edges.insert(e);
			detector.clear();
			detector.temp.clear();
			detector.temp.AddAlternative(e);
			detector.AddNewEdge(e, 0);

			BidirectionalPath newPath;

			newPath.push_back(e);
			newPath.uid = g.int_id(e);
			RecountLengthsForward(g, newPath, lengths);

			//Extend trivially
			ExtendTrivialForward(g, newPath, detector, &lengths, pairedInfo);

			e = g.conjugate(e);
			edges.insert(e);
			detector.clear();
			detector.temp.clear();
			detector.temp.AddAlternative(e);
			detector.AddNewEdge(e, 0);

			BidirectionalPath conjPath;
			conjPath.uid = g.int_id(e);
			conjPath.push_back(e);

			RecountLengthsForward(g, conjPath, lengths);

			//Extend trivially
			ExtendTrivialForward(g, conjPath, detector, &lengths, pairedInfo);

			AddPathPairToContainer(newPath, conjPath, seeds);
		}
	}

	for (auto iter = seeds.begin(); iter != seeds.end(); ++iter) {
		if (std::abs(iter->conj_id - iter->id) != 1) {
			DEBUG("Complement ids are wrong.");
		}
	}

	INFO("Extending seeds backward");
	for (auto pathIter = seeds.begin(); pathIter != seeds.end(); ++pathIter) {
		SimpleRecountDetectorBackward(*pathIter, detector);
		RecountLengthsBackward(g, *pathIter, lengths);

		ExtendTrivialBackward(g, *pathIter, detector, &lengths, pairedInfo);
	}

	INFO("Finding seeds finished");
}

} // namespace long_contigs


#endif /* SEEDS_HPP_ */

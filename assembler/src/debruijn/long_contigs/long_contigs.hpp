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
#include "../../common/logging.hpp"

#include "lc_config.hpp"

namespace long_contigs {

using namespace debruijn_graph;

//Heuristic constants
static size_t READ_SIZE;

static int DISTANCE_DEV;
static size_t WEIGHT_TRESHOLD;

static size_t MAX_LOOPS;
static bool FULL_LOOP_REMOVAL;

static bool ALL_SEEDS;
static double EDGE_COVERAGE_TRESHOLD;
static double LENGTH_COVERAGE_TRESHOLD;


//Deque used for extending path in both directions
typedef std::deque<EdgeId> BidirectionalPath;
typedef std::deque<double> PathLengths;
typedef std::multimap<EdgeId, std::pair<size_t, double> > CycleDetector;

// ====== Support functions ======

//Heuristic constants loader
static void LoadLCConstants() {
	checkFileExistenceFATAL(LC_CONFIG_FILENAME);
	checkFileExistenceFATAL(CONFIG_FILENAME);

	READ_SIZE = LC_CONFIG.read<size_t>("read_size");

	DISTANCE_DEV = CONFIG.read<bool>("etalon_info_mode") ? LC_CONFIG.read<int>("etalon_distance_dev") : LC_CONFIG.read<int>("real_distance_dev");
	WEIGHT_TRESHOLD = LC_CONFIG.read<size_t>("weight_threshold");

	MAX_LOOPS = LC_CONFIG.read<size_t>("max_loops");
	FULL_LOOP_REMOVAL = LC_CONFIG.read<bool>("full_loop_removal");

	ALL_SEEDS = LC_CONFIG.read<bool>("all_seeds");
	EDGE_COVERAGE_TRESHOLD = LC_CONFIG.read<double>("edge_coverage");
	LENGTH_COVERAGE_TRESHOLD = LC_CONFIG.read<double>("len_coverage");

}

//Pause to see output
void MakeKeyPause() {
	int v;
	std::cin >> v;
}

//Edge count in graph
size_t EdgeCount(Graph& g) {
	size_t edgeCount = 0;
	for (auto iter = g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
		++edgeCount;
	}
	return edgeCount;
}

//Path length
size_t PathLength(Graph& g, BidirectionalPath& path) {
	double currentLength = 0;

	for(auto iter = path.begin(); iter != path.end(); ++iter) {
		currentLength += g.length(*iter);
	}
	return currentLength;
}

//Statistic functions
//Prints number of edges and amount of paths with respective edge count
void PrintPathEdgeLengthStats(std::vector<BidirectionalPath>& paths) {
	std::map<size_t, size_t> lengthMap;

	for(auto iter = paths.begin(); iter != paths.end(); ++iter) {
		++lengthMap[iter->size()];
	}

	INFO("Edges in path : path count");
	for(auto iter = lengthMap.begin(); iter != lengthMap.end(); ++iter) {
		INFO(iter->first << " : " << iter->second);
	}
}

//Prints length and amount of paths with respective length
void PrintPathLengthStats(Graph& g, std::vector<BidirectionalPath>& paths) {
	std::map<size_t, size_t> lengthMap;

	for(auto iter = paths.begin(); iter != paths.end(); ++iter) {
		++lengthMap[PathLength(g, *iter)];
	}

	INFO("Path length : paths count");
	for(auto iter = lengthMap.begin(); iter != lengthMap.end(); ++iter) {
		INFO(iter->first << " : " << iter->second);
	}
}

//Prints coverage of all edges by given paths and total edge coverage
double PrintPathCoverage(Graph& g, std::vector<BidirectionalPath>& paths) {
	std::multiset<EdgeId> covered;

	for(auto path = paths.begin(); path != paths.end(); ++path) {
		for(auto iter = path->begin(); iter != path->end(); ++iter) {
			covered.insert(*iter);
		}
	}

	std::map<size_t, size_t> coveredTimes;
	size_t edgeCount = 0;
	for (auto iter = g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
		++coveredTimes[covered.count(*iter)];
		++edgeCount;
	}

	INFO("Covered times : edges")
	for(auto iter = coveredTimes.begin(); iter != coveredTimes.end(); ++iter) {
		INFO(iter->first << " : " << iter->second);
	}

	INFO("Total edge coverage: " << edgeCount - coveredTimes[0] << " out of " << edgeCount);

	return 1.0 - (double) (coveredTimes[0]) / double (edgeCount);
}

//Percentage of edges covered by paths
double PathsCoverage(Graph& g, std::vector<BidirectionalPath>& paths) {
	std::set<EdgeId> covered;

	for(auto path = paths.begin(); path != paths.end(); ++path) {
		for(auto iter = path->begin(); iter != path->end(); ++iter) {
			covered.insert(*iter);
		}
	}

	size_t edgeCount = 0;
	for (auto iter = g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
		++edgeCount;
	}

	return (double) covered.size() / (double) edgeCount;
}

//Percentage of nucleotides covered by paths
double PathsLengthCoverage(Graph& g, std::vector<BidirectionalPath>& paths) {
	std::set<EdgeId> covered;

	for(auto path = paths.begin(); path != paths.end(); ++path) {
		for(auto iter = path->begin(); iter != path->end(); ++iter) {
			covered.insert(*iter);
		}
	}

	size_t pathsLength = 0;
	for(auto iter = covered.begin(); iter != covered.end(); ++iter) {
		pathsLength += g.length(*iter);
	}

	size_t length = 0;
	for (auto iter = g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
		length += g.length(*iter);
	}

	return (double) pathsLength / (double) length;
}

//Print short info about paths paths (length and edge count)
void PrintPathsShort(Graph& g, std::vector<BidirectionalPath>& paths) {
	INFO("Seed length : edge count")
	for(auto iter = paths.begin(); iter != paths.end(); ++iter) {
		INFO(PathLength(g, *iter) << " : " << iter->size());
	}
}

//Print path with length from start / end to the every edge
void PrintPath(Graph& g, BidirectionalPath& path, PathLengths& lengths) {
	INFO("Path " << &path)
	INFO("#, edge, length, total length")
	for(size_t i = 0; i < path.size(); ++i) {
		INFO(i << ", " << path[i] << ", " << g.length(path[i]) << ", " << lengths[i]);
	}
}

//Print path
void PrintPath(Graph& g, BidirectionalPath& path) {
	INFO("Path " << &path)
	INFO("#, edge, length")
	for(size_t i = 0; i < path.size(); ++i) {
		INFO(i << ", " << path[i] << ", " << g.length(path[i]));
	}
}

//Print cycle detector data
void PrintDetector(CycleDetector& detector) {
	INFO("Detector data")
	for(auto iter = detector.begin(); iter != detector.end(); ++iter) {
		INFO("Edge " << iter->first << " comes when path length is " << iter->second.first << " with weight " << iter->second.second);
	}
}



// ====== Seed functions ======

//Extends trivial path forward
//If a start of another trivial path is found, returns it
//Otherwise returns 0
EdgeId ExtendTrivialForward(Graph& g, BidirectionalPath& path, const std::map<EdgeId, BidirectionalPath>& starts) {
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
	while (g.CheckUniqueIncomingEdge(currentVertex)) {
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
	std::set<EdgeId> visited;
	std::map<EdgeId, BidirectionalPath> starts;
	int count = 0;

	INFO("Finding seeds started")
	for (auto iter = g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
		count++;
		EdgeId e = *iter;

		//Try to make seed starting from unvisited edge
		if (visited.count(e) == 0) {
			visited.insert(e);

			starts[e] = BidirectionalPath();
			BidirectionalPath& newPath = starts[e];
			newPath.push_back(e);

			//Extend trivially
			EdgeId nextStart = ExtendTrivialForward(g, newPath, starts);
			for (auto edgeInPath = newPath.begin(); edgeInPath != newPath.end(); ++edgeInPath) {
				visited.insert(*edgeInPath);
			}

			//If extended till another seed, than concatenate them
			if (nextStart != 0) {
				INFO("Join paths");
				JoinPaths(newPath, starts[nextStart]);
				starts.erase(nextStart);
			}
		}
	}

	//Debug part
	for (auto pathIter = starts.begin(); pathIter != starts.end(); ++pathIter) {
		seeds.push_back(pathIter->second);
	}
	PrintPathCoverage(g, seeds);
	//End of debug part

	//Extending seed backward
	seeds.clear();
	seeds.reserve(starts.size());
	INFO("Extending seeds backward");
	for (auto pathIter = starts.begin(); pathIter != starts.end(); ++pathIter) {
		ExtendTrivialBackward(g, pathIter->second);
		seeds.push_back(pathIter->second);
	}
	INFO("Finding seeds finished")
}

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
		INFO("Calculating weight")
		double weight = ExtentionWeight(g, path, lengths, *iter, pairedInfo, forward);
		INFO("Weight " << weight)

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

	INFO("Growing done")
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


// ====== Quality functions ======
//Find bidirectional path in given genome path
size_t FindInGenomePath(BidirectionalPath& myPath, Path<Graph::EdgeId>& genomePath) {
	if (myPath.size() > genomePath.size()) {
		INFO("Warning, unexpected path length")
		return -1;
	}

	for (size_t i = 0; i < genomePath.size() - myPath.size() + 1; ++i) {
		bool found = true;

		size_t j = i;
		for(auto iter = myPath.begin(); iter != myPath.end(); ++iter) {
			if (*iter != genomePath[j]) {
				found = false;
				break;
			}
			++j;
		}

		if (found) {
			return i;
		}
	}
	return -1;
}

//Find inexact match to genome path
size_t FindInGenomeInexact(Graph& g, BidirectionalPath& myPath, Path<Graph::EdgeId>& genomePath, int& startPos, size_t& maxLengthMached) {
	if (myPath.size() > genomePath.size()) {
		INFO("Warning, unexpected path length")
		return -1;
	}

	size_t maxEdgesMatched = 0;
	for (size_t i = 0; i < genomePath.size() - myPath.size() + 1; ++i) {

		size_t j = i;
		size_t edgesMatched = 0;
		size_t lengthMatched = 0;

		for(auto iter = myPath.begin(); iter != myPath.end(); ++iter) {
			if (*iter == genomePath[j]) {
				++edgesMatched;
				lengthMatched += g.length(*iter);
			}
			++j;
		}

		if (edgesMatched > maxEdgesMatched) {
			maxEdgesMatched = edgesMatched;
			maxLengthMached = lengthMatched;
			startPos = i;
		}
	}

	return maxEdgesMatched;
}


//Count all paths in genome paths
template<size_t k>
size_t PathsInGenome(Graph& g, const EdgeIndex<k + 1, Graph>& index, const Sequence& genome, std::vector<BidirectionalPath>& paths) {
	Path<typename Graph::EdgeId> path1 = FindGenomePath<k> (genome, g, index);
	Path<typename Graph::EdgeId> path2 = FindGenomePath<k> (!genome, g, index);

	size_t pathCount = 0;
	for(auto iter = paths.begin(); iter != paths.end(); ++iter) {

		int s = FindInGenomePath(*iter, path1);
		if (s != -1) {
			++pathCount;
			INFO("Path of length " << PathLength(g, *iter)  << " with " << iter->size() << " edges is found in genome path starting from edge " << s)
		}

		else {
			s = FindInGenomePath(*iter, path2);
			if (s != -1) {
				++pathCount;
				INFO("Path of length " << PathLength(g, *iter) << " with " << iter->size() << " edges is found in !genome path starting from edge " << s)
			}

			else {
				int pos1 = 0, pos2 = 0;
				size_t len1 = 0, len2 = 0;
				size_t edges1 = FindInGenomeInexact(g, *iter, path1, pos1, len1);
				size_t edges2 = FindInGenomeInexact(g, *iter, path2, pos2, len2);

				if (edges1 > edges2) {
					INFO("Path partly found, percent of edges matched " << (double) edges1 / iter->size() <<
							", length percentage " << (double) len1 / PathLength(g, *iter));
				}
				else {
					INFO("Path partly found, percent of edges matched " << (double) edges2 / iter->size() <<
												", length percentage " << (double) len2 / PathLength(g, *iter));
				}

			}
		}
	}
	return pathCount;
}


}  // namespace long_contigs

#endif /* LONG_CONTIGS_HPP_ */

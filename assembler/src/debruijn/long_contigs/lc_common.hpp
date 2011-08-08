/*
 * lc_common.hpp
 *
 *  Created on: Aug 3, 2011
 *      Author: andrey
 */

#ifndef LC_COMMON_HPP_
#define LC_COMMON_HPP_

#include <vector>
#include <deque>
#include <set>
#include <map>

#include "../launch.hpp"
#include "../../common/logging.hpp"
#include "../new_debruijn.hpp"
#include "../config.hpp"

#include "lc_config.hpp"

namespace long_contigs {

using namespace debruijn_graph;

//Deque used for extending path in both directions
typedef std::deque<EdgeId> BidirectionalPath;
typedef std::deque<double> PathLengths;
typedef std::multimap<EdgeId, std::pair<size_t, double> > CycleDetector;

struct PairedInfoIndexLibrary {

	PairedInfoIndexLibrary(size_t readS, size_t insS, PairedInfoIndex<Graph>* index): readSize(readS), insertSize(insS) , pairedInfoIndex(index) {
	}

	size_t readSize;
	size_t insertSize;
	PairedInfoIndex<Graph>* pairedInfoIndex;
};

typedef std::vector<PairedInfoIndexLibrary> PairedInfoIndices;

// ====== Support functions ======
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

	INFO("Total paths " << paths.size());
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

	INFO("Total paths " << paths.size());
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
	INFO("Total paths " << paths.size());
	INFO("Path length : edge count")
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
template<class PathType>
void PrintPath(Graph& g, PathType& path) {
	INFO("Path " << &path)
	INFO("#, edge, length")
	for(size_t i = 0; i < path.size(); ++i) {
		INFO(i << ", " << path[i] << ", " << g.length(path[i]));
	}
}

//Print path
template<class PathType>
void PrintPathWithVertices(Graph& g, PathType& path) {
	INFO("Path " << &path)
	INFO("#, edge, length")

	for(size_t i = 0; i < path.size(); ++i) {
		INFO(g.EdgeStart(path[i]));
		INFO(i << ", " << path[i] << ", " << g.length(path[i]));
	}
}

//[,)
void PrintPathFromTo(Graph& g, Path<Graph::EdgeId>& path, size_t startPos = 0, size_t endPos = 0) {
	if (startPos >= path.size() || endPos > path.size()) return;

	if (endPos == 0) {
		endPos = path.size();
	}

	INFO("Path of length " << path.size() << " from " << startPos << " to " << endPos);
	for (size_t i = startPos; i < endPos; ++i) {
		INFO(i << ", " << path[i] << ", " << g.length(path[i]));
	}
}

//[,)
void PrintPathFromTo(Graph& g, BidirectionalPath& path, size_t startPos = 0, size_t endPos = 0) {
	if (startPos >= path.size() || endPos > path.size()) return;

	if (endPos == 0) {
		endPos = path.size();
	}

	INFO("Path of length " << path.size() << " from " << startPos << " to " << endPos);
	for (size_t i = startPos; i < endPos; ++i) {
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

bool ComparePaths(const BidirectionalPath& path1, const BidirectionalPath& path2) {
	if (path1.size() != path2.size()) {
		return false;
	}

	for (size_t i = 0; i < path1.size(); ++i) {
		if (path1[i] != path2[i]) {
			return false;
		}
	}
	return true;
}

bool ContainsPath(const BidirectionalPath& path, const BidirectionalPath& sample) {
	if (path.size() < sample.size()) {
		return false;
	}

	for (size_t i = 0; i < path.size() - sample.size() + 1 ; ++i) {
		bool found = true;

		for (size_t j = 0; j < sample.size(); ++j) {
			if (sample[j] != path[i + j]) {
				found = false;
				break;
			}
		}

		if (found) {
			return true;
		}
	}
	return false;
}

//Find coverage of worst covered edge
double PathMinReadCoverage(Graph& g, BidirectionalPath& path) {
	if (path.empty()) {
		return 0;
	}

	double minCov = g.coverage(path[0]);

	for (auto edge = path.begin(); edge != path.end(); ++edge) {
		double cov = g.coverage(*edge);
		if (minCov > cov) {
			minCov = cov;
		}
	}

	return minCov;
}

//Remove paths with low covered edges
void FilterLowCovered(Graph& g, std::vector<BidirectionalPath>& paths, double threshold) {
	for (auto path = paths.begin(); path != paths.end(); ) {
		if (PathMinReadCoverage(g, *path) < threshold) {
			paths.erase(path);
		} else {
			++path;
		}
	}
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

//Remove subpaths paths
void RemoveSubpaths(const std::vector<BidirectionalPath>& paths, std::vector<BidirectionalPath>& output) {
	for (auto path = paths.begin(); path != paths.end(); ++path) {
		bool copy = true;
		for (auto iter = output.begin(); iter != output.end(); ++iter) {
			if (ContainsPath(*iter, *path)) {
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


#endif /* LC_COMMON_HPP_ */

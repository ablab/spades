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

#include "lc_config_struct.hpp"
#include "../launch.hpp"
#include "logging.hpp"
#include "../new_debruijn.hpp"

namespace long_contigs {

using namespace debruijn_graph;
using debruijn::K;

//Deque used for extending path in both directions
typedef std::deque<EdgeId> BidirectionalPath;

//Forward declaration
size_t PathLength(Graph& g, const BidirectionalPath& path);

//Path cumulative lengths
typedef std::deque<double> PathLengths;

//Simple cycle detector
typedef std::multimap<EdgeId, std::pair<size_t, double> > CycleDetector;

//Paired infi library
struct PairedInfoIndexLibrary {

	PairedInfoIndexLibrary(size_t readS, size_t insS, size_t v, PairedInfoIndex<Graph>* index): readSize(readS), insertSize(insS), var(v), pairedInfoIndex(index) {
	}

	size_t readSize;
	size_t insertSize;
	size_t var;
	PairedInfoIndex<Graph>* pairedInfoIndex;
};

typedef std::vector<PairedInfoIndexLibrary> PairedInfoIndices;


//Statistics
enum StopReason { LOOP, NO_EXTENSION, NO_GOOD_EXTENSION, MANY_GOOD_EXTENSIONS, WEAK_EXTENSION };

struct PathStatData {
	StopReason reason;
	size_t pathLength;
	std::string message;

	PathStatData(StopReason rsn, size_t pathLen, const std::string& msg): reason(rsn), pathLength(pathLen), message(msg) {
	}
};

class PathStopHandler {
private:
	Graph& g_;
	std::multimap<const BidirectionalPath*, PathStatData> forward_;
	std::multimap<const BidirectionalPath*, PathStatData> backward_;

public:
	PathStopHandler(Graph& g):  g_(g), forward_(), backward_() {}

	void AddStop(const BidirectionalPath& path, StopReason reason, bool forward) {
		std::string msg;
		switch (reason) {
		case LOOP: {
			msg = "cycle detected";
			break;
		}
		case NO_EXTENSION: {
			msg = "no extension found";
			break;
		}
		case MANY_GOOD_EXTENSIONS: {
			msg = "two or more possible extensions are very similar";
			break;
		}
		case NO_GOOD_EXTENSION: {
			msg = "all extensions are bad";
			break;
		}
		case WEAK_EXTENSION: {
			msg = "bast extension's weight does not pass the threshold";
			break;
		}
		default: {
			msg = "";
		}
		}

		AddStop(path, reason, forward, msg);
	}

	void AddStop(const BidirectionalPath& path, StopReason reason, bool forward, const std::string& msg) {
		if (forward) {
			forward_.insert(std::make_pair(&path, PathStatData(reason, PathLength(g_, path), msg)));
		} else {
			backward_.insert(std::make_pair(&path, PathStatData(reason, PathLength(g_, path), msg)));
		}
	}

	void print(const BidirectionalPath& path) {
		INFO("Stats for path with " << path.size() << " edges and length " <<  PathLength(g_, path));
		INFO("Stoppages forward (" << forward_.count(&path) << "):");
		for (auto iter = forward_.lower_bound(&path); iter != forward_.upper_bound(&path); ++iter) {
			INFO("Stop reason at length " << iter->second.pathLength << ", reason: " << iter->second.message);
		}
		INFO("Stoppages backward (" << backward_.count(&path) << "):");
		for (auto iter = backward_.lower_bound(&path); iter != backward_.upper_bound(&path); ++iter) {
			INFO("Stop reason at length " << iter->second.pathLength << ", reason: " << iter->second.message);
		}
	}

	void print() {
		for (auto iter = forward_.begin(); iter != forward_.end(); ++iter) {
			print(*(iter->first));
		}
	}
};

//Detailed output for research mode
#define DETAILED_INFO(message) { if (lc_cfg::get().rs.detailed_output) { INFO(message) } }

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
size_t PathLength(Graph& g, const BidirectionalPath& path) {
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

//Print path with length from start / end to the every edge
void DetailedPrintPath(Graph& g, BidirectionalPath& path, PathLengths& lengths) {
	if (lc_cfg::get().rs.detailed_output) {
		PrintPath(g, path, lengths);
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
void DetailedPrintPath(Graph& g, PathType& path) {
	if (lc_cfg::get().rs.detailed_output) {
		PrintPath(g, path);
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

//Ouput only edges with specified length
void PrintEdgeNuclsByLength(Graph& g, size_t edgeLen) {
	for (auto edge = g.SmartEdgeBegin(); !edge.IsEnd(); ++edge) {
		if (g.length(*edge) == edgeLen) {
			INFO("Length " << edgeLen << ", Data: " << g.EdgeNucls(*edge).Subseq(0, g.length(*edge) + 1).str());
		}
	}
}

} // namespace long_contigs


#endif /* LC_COMMON_HPP_ */

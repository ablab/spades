//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * lc_common.hpp
 *
 *  Created on: Aug 3, 2011
 *      Author: andrey
 */

#ifndef LC_COMMON_HPP_
#define LC_COMMON_HPP_


#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <deque>

#include "../standard.hpp"
#include "../config_struct.hpp"
#include "lc_config_struct.hpp"


#include "logger/logger.hpp"

#include "io/reader.hpp"
#include "io/easy_reader.hpp"
#include "io/rc_reader_wrapper.hpp"
#include "io/cutting_reader_wrapper.hpp"
#include "io/multifile_reader.hpp"
#include "io/careful_filtering_reader_wrapper.hpp"

#include "../standard.hpp"
#include "../new_debruijn.hpp"
#include "../graphio.hpp"
#include "../graph_construction.hpp"
#include "../graph_simplification.hpp"

#include "omni/distance_estimation.hpp"
//#include "omni/advanced_distance_estimation.hpp"

namespace long_contigs {

DECL_LOGGER("lc");

using namespace debruijn_graph;

//Deque used for extending path in both directions
class BidirectionalPath : public std::deque<EdgeId> {
public:
	 std::deque<int> gaps;

	 int uid;

	 int id;
	 int conj_id;

};

//Path length
template<class T>
size_t PathLength(const Graph& g, T& path) {
	double currentLength = 0;

	for(int i = 0; i < (int) path.size(); ++i) {
		currentLength += g.length(path[i]);
	}
	return currentLength;
}

const Path<EdgeId> as_simple_path(const Graph& g, const BidirectionalPath& path) {
	 return Path<EdgeId>(vector<EdgeId>(path.begin(), path.end()), 0, g.length(path.back()));
}

const MappingPath<EdgeId> as_trivial_mapping_path(const Graph& g, const Path<EdgeId>& path) {
	vector<MappingRange> ranges;
	size_t cum_length = 0;
	for (auto it = path.begin(); it != path.end(); ++it) {
		ranges.push_back(MappingRange(Range(cum_length, cum_length + g.length(*it))
				, Range(0, g.length(*it))));
		cum_length += g.length(*it);
	}
	return MappingPath<EdgeId>(vector<EdgeId>(path.begin(), path.end()), ranges);
}

//Path cumulative lengths
typedef std::deque<double> PathLengths;

//Simple cycle detector
typedef std::multimap<EdgeId, std::pair<size_t, double> > CycleDetector;

//Paired infi library
struct PairedInfoIndexLibrary {

	Graph const * g_;
	size_t readSize;
	size_t insertSize;
	size_t is_delta;
	size_t deDelta;
	size_t var;
	PairedInfoIndex<Graph> * pairedInfoIndex;

	bool has_advanced;
	PairedInfoIndexLibrary* advanced;

	PairedInfoIndex<Graph>* raw;

	PairedInfoIndexLibrary(const Graph& g, size_t readS, size_t insS, size_t delta, size_t ded, size_t v, PairedInfoIndex<Graph>* index):
		g_(&g), readSize(readS), insertSize(insS), is_delta(delta), deDelta(ded), var(v), pairedInfoIndex(index), has_advanced(false), advanced(0), raw(0) {
	}

	double NormalizeWeight(const PairInfo<EdgeId>& pair_info) {
		double w = 0.;
		if (math::eq(pair_info.d, 0.) && pair_info.first == pair_info.second) {
			w = 0. + g_->length(pair_info.first) - insertSize
					+ 2 * readSize + 1 - K;
		} else {
			EdgeId e1 =	(math::ge(pair_info.d, 0.)) ?
							pair_info.first : pair_info.second;
			EdgeId e2 =	(math::ge(pair_info.d, 0.)) ?
							pair_info.second : pair_info.first;

			int gap_len = std::abs(rounded_d(pair_info)) - g_->length(e1);
			int right = std::min(insertSize, gap_len + g_->length(e2) + readSize);
			int left = std::max(gap_len, int(insertSize) - int(readSize) - int(g_->length(e1)));
			w = 0. + right - left + 1 - K;
		}

		double result_weight = pair_info.weight;
		if (math::gr(w, 0.)) {
			result_weight /= w;
		} else {
			result_weight = 1;
		}

		return result_weight;
	}
};

typedef std::vector<PairedInfoIndexLibrary> PairedInfoIndices;


class SimplePathComparator {
private:
	const Graph& g_;

public:
	SimplePathComparator(const Graph& g): g_(g) {}

	bool operator() (const BidirectionalPath& path1, const BidirectionalPath& path2) const {
		return PathLength(g_, path1) > PathLength(g_, path2);
	}

	bool operator() (const BidirectionalPath* path1, const BidirectionalPath* path2) const {
		return PathLength(g_, *path1) > PathLength(g_, *path2);
	}
};

template <class T1, class T2>
class SimplePairComparator {
public:

	bool operator() (const std::pair<T1, T2>& p1, const std::pair<T1, T2>& p2) const {
		return p1.first > p2.first;
	}

	bool operator() (const std::pair<T1, T2>* p1, const std::pair<T1, T2>* p2) const {
		return p1->first > p2->first;
	}
};


//Statistics
enum StopReason { LOOP, LONG_LOOP, NO_EXTENSION, NO_GOOD_EXTENSION, MANY_GOOD_EXTENSIONS, WEAK_EXTENSION };

struct PathStatData {
	StopReason reason;
	size_t pathLength;
	std::string message;

	PathStatData(StopReason rsn, size_t pathLen, const std::string& msg): reason(rsn), pathLength(pathLen), message(msg) {
	}
};

class PathStopHandler {
private:
	const Graph& g_;
	std::multimap<BidirectionalPath*, PathStatData> forward_;
	std::multimap<BidirectionalPath*, PathStatData> backward_;

public:
	PathStopHandler(const Graph& g):  g_(g), forward_(), backward_() {}

	void AddStop(BidirectionalPath* path, StopReason reason, bool forward) {
		std::string msg;
		switch (reason) {
		case LOOP: {
			msg = "cycle detected";
			break;
		}
		case LONG_LOOP: {
			msg = "loop is too long to resolve correctly";
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

	void AddStop(BidirectionalPath* path, StopReason reason, bool forward, const std::string& msg) {
		if (forward) {
			forward_.insert(std::make_pair(path, PathStatData(reason, PathLength(g_, *path), msg)));
		} else {
			backward_.insert(std::make_pair(path, PathStatData(reason, PathLength(g_, *path), msg)));
		}
	}

	void print(BidirectionalPath* path) {
		DEBUG("Stats for path " << path << " with " << path->size() << " edges and length " <<  PathLength(g_, *path));
		DEBUG("Stoppages forward (" << forward_.count(path) << "):");
		for (auto iter = forward_.lower_bound(path); iter != forward_.upper_bound(path); ++iter) {
			DEBUG("Stop reason at length " << iter->second.pathLength << ", reason: " << iter->second.message);
		}
		DEBUG("Stoppages backward (" << backward_.count(path) << "):");
		for (auto iter = backward_.lower_bound(path); iter != backward_.upper_bound(path); ++iter) {
			DEBUG("Stop reason at length " << iter->second.pathLength << ", reason: " << iter->second.message);
		}
	}

	void print() {
		std::set<BidirectionalPath*> printed;

		for (auto iter = forward_.begin(); iter != forward_.end(); ++iter) {
			if (printed.count(iter->first) == 0) {
				printed.insert(iter->first);
				print(iter->first);
			}
		}
	}
};

//Detailed output for research mode
#define DETAILED_DEBUG(message) { /*if (lc_cfg::get().params.rs.detailed_output) {*/ DEBUG(message) /*}*/ }

// ====== Support functions ======
//Pause to see output
void MakeKeyPause() {
	int v;
	std::cin >> v;
}

//Edge count in graph
size_t EdgeCount(const Graph& g) {
	size_t edgeCount = 0;
	for (auto iter = g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
		++edgeCount;
	}
	return edgeCount;
}


//Prints number of edges and amount of paths with respective edge count
void PrintPathEdgeLengthStats(std::vector<BidirectionalPath>& paths) {
	std::map<size_t, size_t> lengthMap;

	for(auto iter = paths.begin(); iter != paths.end(); ++iter) {
		++lengthMap[iter->size()];
	}

	DEBUG("Total paths " << paths.size());
	DEBUG("Edges in path : path count");
	for(auto iter = lengthMap.begin(); iter != lengthMap.end(); ++iter) {
		DEBUG(iter->first << " : " << iter->second);
	}
}

//Prints length and amount of paths with respective length
void PrintPathLengthStats(const Graph& g, std::vector<BidirectionalPath>& paths) {
	std::map<size_t, size_t> lengthMap;

	for(auto iter = paths.begin(); iter != paths.end(); ++iter) {
		++lengthMap[PathLength(g, *iter)];
	}

	DEBUG("Total paths " << paths.size());
	DEBUG("Path length : paths count");
	for(auto iter = lengthMap.begin(); iter != lengthMap.end(); ++iter) {
		DEBUG(iter->first << " : " << iter->second);
	}
}



//Prints coverage of all edges by given paths and total edge coverage
double PrintPathCoverage(const Graph& g, std::vector<BidirectionalPath>& paths) {
	std::multiset<EdgeId> covered;

	for(auto path = paths.begin(); path != paths.end(); ++path) {
		for(auto iter = path->begin(); iter != path->end(); ++iter) {
			covered.insert(*iter);
		}
	}

	DEBUG("EDGE STAT")
	for (auto iter = covered.begin(); iter != covered.end(); ++iter) {
	    DEBUG(g.int_id(*iter) << " (" << g.length(*iter) << ") = " << covered.count(*iter));

	}

	std::map<size_t, size_t> coveredTimes;
	size_t edgeCount = 0;
	for (auto iter = g.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
		++coveredTimes[covered.count(*iter)];
		++edgeCount;
	}

	DEBUG("Covered times : edges")
	for(auto iter = coveredTimes.begin(); iter != coveredTimes.end(); ++iter) {
		DEBUG(iter->first << " : " << iter->second);
	}

	DEBUG("Total edge coverage: " << edgeCount - coveredTimes[0] << " out of " << edgeCount);

	return 1.0 - (double) (coveredTimes[0]) / double (edgeCount);
}

//Percentage of edges covered by paths
double PathsCoverage(const Graph& g, std::vector<BidirectionalPath>& paths) {
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
double PathsLengthCoverage(const Graph& g, std::vector<BidirectionalPath>& paths) {
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
void PrintPathsShort(const Graph& g, const std::vector<BidirectionalPath>& paths) {
	DEBUG("Total paths " << paths.size());
	DEBUG("Path length : edge count")
	for(auto iter = paths.begin(); iter != paths.end(); ++iter) {
		DEBUG(PathLength(g, *iter) << " : " << iter->size());
	}
}

//Print path with length from start / end to the every edge
void PrintPath(const Graph& g, BidirectionalPath& path, PathLengths& lengths) {
	DEBUG("Path " << path.uid)
	DEBUG("#, edge, length, total length")
	for(size_t i = 0; i < path.size(); ++i) {
		DEBUG(i << ", " << g.int_id(path[i]) << ", " << g.length(path[i]) << ", " << lengths[i]);
	}
}

//Print path with length from start / end to the every edge
void DetailedPrintPath(const Graph& g, BidirectionalPath& path, PathLengths& lengths) {
//	if (lc_cfg::get().params.rs.detailed_output) {
		PrintPath(g, path, lengths);
//	}
}

//Print path
template<class PathType>
void PrintAnyPath(const Graph& g, const PathType& path) {
	DEBUG("Path " << &path << " with length " << PathLength(g, path));
	DEBUG("#, edge, length")
	for(int i = 0; i < (int) path.size(); ++i) {
		DEBUG(i << ", " << g.int_id(path[i]) << ", " << g.length(path[i]));
	}
}

void PrintPath(const Graph& g, const BidirectionalPath& path) {
    INFO("Path " << path.uid << " with length " << PathLength(g, path));
    INFO("#, edge, length, coverage")
    for(int i = 0; i < (int) path.size(); ++i) {
        INFO(i << ", " << g.int_id(path[i]) << ", len: " << g.length(path[i]) << ", cov: " << g.coverage(path[i]));
    }
}

//Print path
template<class PathType>
void DetailedPrintAnyPath(const Graph& g, PathType& path) {
	if (lc_cfg::get().params.rs.detailed_output) {
		PrintPath(g, path);
	}
}

void DetailedPrintPath(const Graph& g, BidirectionalPath& path) {
    if (lc_cfg::get().params.rs.detailed_output) {
        PrintPath(g, path);
    }
}

//Print path
template<class PathType>
void PrintPathWithVertices(const Graph& g, PathType& path) {
	DEBUG("Path " << path.uid)
	DEBUG("#, edge, length")

	for(size_t i = 0; i < path.size(); ++i) {
		DEBUG(g.int_id(g.EdgeStart(path[i])));
		DEBUG(i << ", " << g.int_id(path[i]) << ", " << g.length(path[i]));
	}
}

void CountPathLengths(const Graph& g, std::vector<BidirectionalPath>& paths, std::vector<size_t>& lengths) {
	lengths.clear();
	for (auto path = paths.begin(); path != paths.end(); ++path) {
			lengths.push_back(PathLength(g, *path));
	}
}


//[,)
void PrintPathFromTo(const Graph& g, Path<Graph::EdgeId>& path, size_t startPos = 0, size_t endPos = 0) {
	if (startPos >= path.size() || endPos > path.size()) return;

	if (endPos == 0) {
		endPos = path.size();
	}

	DEBUG("Path of length " << path.size() << " from " << startPos << " to " << endPos);
	for (size_t i = startPos; i < endPos; ++i) {
		DEBUG(i << ", " << path[i] << ", " << g.length(path[i]));
	}
}

//[,)
void PrintPathFromTo(const Graph& g, BidirectionalPath& path, size_t startPos = 0, size_t endPos = 0) {
	if (startPos >= path.size() || endPos > path.size()) return;

	if (endPos == 0) {
		endPos = path.size();
	}

	DEBUG("Path of length " << path.size() << " from " << startPos << " to " << endPos);
	for (size_t i = startPos; i < endPos; ++i) {
		DEBUG(i << ", " << path[i] << ", " << g.length(path[i]));
	}
}

//Print cycle detector data
void PrintDetector(CycleDetector& detector) {
	DEBUG("Detector data");
	for(auto iter = detector.begin(); iter != detector.end(); ++iter) {
		DEBUG("Edge " << iter->first << " comes when path length is " << iter->second.first << " with weight " << iter->second.second);
	}
}

//Ouput only edges with specified length
void PrintEdgeNuclsByLength(const Graph& g, size_t edgeLen) {
	for (auto edge = g.SmartEdgeBegin(); !edge.IsEnd(); ++edge) {
		if (g.length(*edge) == edgeLen) {
			DEBUG("Length " << edgeLen << ", Data: " << g.EdgeNucls(*edge).Subseq(0, g.length(*edge) + 1).str());
		}
	}
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

size_t GetMinGapSize(PairedInfoIndices& pairedInfo) {
	size_t minIS = 100000;
	for(auto lib = pairedInfo.begin(); lib != pairedInfo.end(); ++lib) {
		if (minIS > lib->insertSize - 2 * lib->readSize) {
			minIS = lib->insertSize - 2 * lib->readSize;
		}
	}
	return minIS;
}


template <class T>
void AddPathPairToContainer(BidirectionalPath p1, BidirectionalPath p2, T& paths) {
	int newId = paths.size();
	p1.id = newId;
	p2.conj_id = newId;
	p2.id = newId + 1;
	p1.conj_id = newId + 1;
	paths.push_back(p1);
	paths.push_back(p2);
}

//Increase path lengths
void IncreaseLengths(const Graph& g, PathLengths& lengths, EdgeId edge, bool forward) {
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

} // namespace long_contigs


#endif /* LC_COMMON_HPP_ */

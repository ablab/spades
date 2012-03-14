/*
 * loop.hpp
 *
 *  Created on: Aug 19, 2011
 *      Author: andrey
 */

#ifndef LOOP_HPP_
#define LOOP_HPP_

#include "lc_common.hpp"

namespace long_contigs {

//Cycle detector data
struct LoopDetectorData {
	size_t iteration;
	double selfWeight;

	map<EdgeId, double, Graph::Comparator> weights;

	LoopDetectorData(size_t iter, double weight, const Graph::Comparator &comparator): iteration(iter), selfWeight(weight), weights(comparator)  {
	}

	LoopDetectorData(const LoopDetectorData& d) : iteration(d.iteration), selfWeight(d.selfWeight), weights(d.weights){
	}

	LoopDetectorData(const Graph::Comparator &comparator) : weights(comparator) {
	}

	void SetSelectedEdge(size_t iter, double w) {
		iteration = iter;
		selfWeight = w;
	}

	void AddAlternative(EdgeId e, double w = 1) {
		weights.insert(std::make_pair(e,w));
	}

	void clear() {
		weights.clear();
		iteration = 0;
		selfWeight = 0;
	}

	bool operator==(const LoopDetectorData& d) {
		if (selfWeight != d.selfWeight) {
			return false;
		}

		if(this->weights != d.weights) {
			return false;
		}

		return true;
	}
};

struct LoopDetector {
	LoopDetectorData temp;
	std::multimap<EdgeId, LoopDetectorData, Graph::Comparator> data;

	LoopDetector(const Graph &graph) : temp(graph.ReliableComparatorInstance()), data(graph.ReliableComparatorInstance()) {
	}

	void AddNewEdge(EdgeId e, size_t iter, double weight = 1) {
		temp.SetSelectedEdge(iter, weight);
		data.insert(std::make_pair(e, temp));
	}

	void clear() {
		data.clear();
		temp.clear();
	}

	void print(const Graph& g) {
		DEBUG("== Detector data ==");
		for (auto iter = data.begin(); iter != data.end(); ++iter) {
			DEBUG("Edge " << g.length(iter->first) << ", weight " << iter->second.selfWeight << ", iteration " << iter->second.iteration);
			for(auto alt = iter->second.weights.begin(); alt != iter->second.weights.end(); ++alt) {
				DEBUG("Edge " << g.length(alt->first) << ", weight " << alt->second);
			}
		}
	}

};


//Add edge to cycle detector and check
bool CheckCycle(BidirectionalPath& path, EdgeId extension, LoopDetector& detector, size_t loopCount) {
	return detector.data.count(extension) > loopCount;
}

size_t CountLoopEdges(EdgeId lastEdge, LoopDetector& detector) {
	auto iter = detector.data.upper_bound(lastEdge);
	--iter;
	size_t loopSize = iter->second.iteration;
	--iter;
	loopSize -= iter->second.iteration;

	return loopSize;
}


bool PathIsOnlyLoop(BidirectionalPath& path, EdgeId loopEdge, bool forward) {
	EdgeId secondEdge = forward ? path.back() : path.front();
	for (auto edge = path.begin(); edge != path.end(); ++edge) {
		if (*edge != secondEdge && *edge != loopEdge) {
			return false;
		}
	}
	return true;
}

bool PathIsOnlyLoop(BidirectionalPath& path, LoopDetector& detector, bool forward) {
	EdgeId lastEdge = forward ? path.back() : path.front();
	size_t loopSize = CountLoopEdges(lastEdge, detector);
	int start = forward ? path.size() - 1 : loopSize - 1;
	int end = forward ? path.size() - loopSize : 0;

	restricted::set<EdgeId> loopEdges;

	for (int i = start; i >= end; --i) {
		loopEdges.insert(path[i]);
	}

	for (int i = 0; i < (int) path.size(); ++i) {
		if (loopEdges.count(path[i]) == 0) {
			return false;
		}
	}
	return true;
}


size_t CountLoopLength(const Graph& g, BidirectionalPath& path, LoopDetector& detector, bool forward) {
	EdgeId lastEdge = forward ? path.back() : path.front();
	size_t loopSize = CountLoopEdges(lastEdge, detector);

	size_t length = 0;

	size_t start = forward ? path.size() - loopSize : 0;
	size_t end = forward ? path.size() : loopSize;

	for (size_t i = start; i != end; ++i) {
		length += g.length(path[i]);
	}
	return length;
}

//Edges to remove
size_t CountEdgesToRemove(BidirectionalPath& path, EdgeId lastEdge, LoopDetector& detector, bool fullRemoval, size_t loopCount, bool forward) {
	size_t loopSize = CountLoopEdges(lastEdge, detector);
	bool onlyCycle = PathIsOnlyLoop(path, detector, forward);

	if (onlyCycle || path.size() <= loopCount * loopSize + 1) {
		DEBUG("Only loop, loop size " << loopSize << ", loop count " << loopCount);
		return path.size() - loopSize;
	}

	if (fullRemoval) {
		DEBUG("Loop size " << loopSize << ", loop count " << loopCount);
		return loopCount * loopSize + 1;
	} else {
		DEBUG("loop size " << loopSize << ", loop count " << loopCount);
		return (loopCount - 1) * loopSize + 1;
	}
}

//Cut loop forward
void RemoveLoopForward(BidirectionalPath& path, LoopDetector& detector, bool fullRemoval, size_t loopCount) {
	size_t edgesToRemove = CountEdgesToRemove(path, path.back(), detector, fullRemoval, loopCount, true);

	for(size_t i = 0; i < edgesToRemove; ++i) {
		path.pop_back();
	}
}

void RemoveLoopBackward(BidirectionalPath& path, LoopDetector& detector, bool fullRemoval, size_t loopCount) {
	size_t edgesToRemove = CountEdgesToRemove(path, path.front(), detector, fullRemoval, loopCount, false);

	for(size_t i = 0; i < edgesToRemove; ++i) {
		path.pop_front();
	}
}

bool LoopBecameStable(EdgeId e, LoopDetector& detector) {
	if (detector.data.count(e) < 2) {
		DETAILED_DEBUG("Loop still unstable");
		return false;
	}
	auto iter = detector.data.upper_bound(e);
	auto last = --iter;
	auto prev = --iter;

	bool res = prev->second == last->second;

	if (res) {
		DETAILED_DEBUG("Loop became stable");
	} else {
		DETAILED_DEBUG("Loop still unstable");
	}
	return res;
}

size_t CountLoopExits(BidirectionalPath& path, EdgeId e, LoopDetector& detector, bool forward) {
	size_t loopSize = CountLoopEdges(e, detector);
	size_t exits = 0;
	int start = forward ? path.size() - 1 : loopSize - 1;
	int end = forward ? path.size() - loopSize : 0;

	for (int i = start; i >= end; --i) {
		LoopDetectorData& data = detector.data.find(path[i])->second;

		exits += data.weights.size() - 1;
	}
	return exits;
}

EdgeId FindFirstFork(BidirectionalPath& path, EdgeId e, LoopDetector& detector, bool forward) {
	size_t loopSize = CountLoopEdges(e, detector);
	int start = forward ? path.size() - 1 : loopSize - 1;
	int end = forward ? path.size() - loopSize : 0;

	for (int i = start; i >= end; --i) {
		LoopDetectorData& data = detector.data.find(path[i])->second;

		if (data.weights.size() == 2) {
			return path[i];
		}
	}
	return EdgeId(0);
}

EdgeId GetForwardFork(const Graph& g, EdgeId e) {
	VertexId v = g.EdgeStart(e);
	if (g.OutgoingEdgeCount(v) != 2) {
		return EdgeId(0);
	}
	std::vector<EdgeId> edges = g.OutgoingEdges(v);
	if (edges[1] == e) {
		return edges[0];
	} else {
		return edges[1];
	}
}

EdgeId GetBackwardFork(const Graph& g, EdgeId e) {
	VertexId v = g.EdgeEnd(e);
	if (g.IncomingEdgeCount(v) != 2) {
		return EdgeId(0);
	}
	std::vector<EdgeId> edges = g.IncomingEdges(v);
	if (edges[1] == e) {
		return edges[0];
	} else {
		return edges[1];
	}
}

bool EdgesMakeShortLoop(const Graph& g, EdgeId e1, EdgeId e2) {
	return g.EdgeStart(e1) == g.EdgeEnd(e2) && g.EdgeStart(e2) == g.EdgeEnd(e1);
}

EdgeId IsEdgeInShortLoopForward(const Graph& g, EdgeId e) {
	VertexId v = g.EdgeEnd(e);
	auto edges = g.OutgoingEdges(v);
	EdgeId result(0);

	for (auto edge = edges.begin(); edge != edges.end(); ++edge) {
		if (g.EdgeEnd(*edge) == g.EdgeStart(e)) {
			result = *edge;
		}
	}

	if (g.OutgoingEdgeCount(v) == 1 && result != EdgeId(0)) {
		DEBUG("Seems no fork backward: edge " << g.length(e) << ", loops with " << g.length(result) << ". " << g.OutgoingEdgeCount(v));
	}

	return result;
}

EdgeId IsEdgeInShortLoopBackward(const Graph& g, EdgeId e) {
	VertexId v = g.EdgeStart(e);
	auto edges = g.IncomingEdges(v);
	EdgeId result(0);

	for (auto edge = edges.begin(); edge != edges.end(); ++edge) {
		if (g.EdgeStart(*edge) == g.EdgeEnd(e)) {
			result = *edge;
		}
	}

	if (g.IncomingEdgeCount(v) == 1 && result != EdgeId(0)) {
		DEBUG("Seems no fork backward: edge " << g.length(e) << ", loops with " << g.length(result) << ". " << g.IncomingEdgeCount(v));
	}

	return result;
}

size_t GetMaxExitIteration(EdgeId loopEdge, EdgeId loopExit, LoopDetector& detector, std::pair<size_t, size_t> iterRange) {
	auto range = detector.data.equal_range(loopEdge);

	size_t maxIter = 0;
	double maxWeight = 0;
	for (auto iter = range.first; iter != range.second; ++iter) {
		double w = iter->second.weights[loopExit];
		if (w > maxWeight &&
				iter->second.iteration >= iterRange.first && iter->second.iteration <= iterRange.second) {

			maxIter = iter->second.iteration;
			maxWeight = w;
		}
	}
	return maxIter;
}

size_t GetFirstExitIteration(EdgeId loopEdge, EdgeId loopExit, LoopDetector& detector,
		std::pair<size_t, size_t> iterRange, double coeff = params.ps.es.priority_coeff) {
	auto range = detector.data.equal_range(loopEdge);

	size_t maxIter = std::numeric_limits<size_t>::max();
	for (auto iter = range.first; iter != range.second; ++iter) {
		if (iter->second.weights[loopExit] * coeff > iter->second.weights[loopEdge] && maxIter > iter->second.iteration	&&
				iter->second.iteration >= iterRange.first && iter->second.iteration <= iterRange.second) {

			maxIter = iter->second.iteration;
		}
	}
	return maxIter;
}


} //namespace long_contigs


#endif /* LOOP_HPP_ */

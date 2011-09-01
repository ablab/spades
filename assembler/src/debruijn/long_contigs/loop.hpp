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

	std::map<EdgeId, double> weights;

	LoopDetectorData(size_t iter, double weight): iteration(iter), selfWeight(weight), weights()  {
	}

	LoopDetectorData(const LoopDetectorData& d) {
		iteration = d.iteration;
		selfWeight = d.selfWeight;

		d.weights.insert(weights.begin(), weights.end());
	}

	LoopDetectorData(): weights()  {
	}

	void SetSelectedEdge(size_t iter, double w) {
		iteration = iter;
		selfWeight = w;
	}

	void AddAlternative(EdgeId e, double w) {
		weights.insert(std::make_pair(e,w));
	}

	void clear() {
		weights.clear();
		iteration = 0;
		selfWeight = 0;
	}

	bool operator==(const LoopDetectorData& d) {
		if (selfWeight != d.selfWeight || weights.size() != d.weights.size()) {
			return false;
		}

		auto iter2 = d.weights.begin();
		for (auto iter1 = weights.begin(); iter2 != d.weights.end() && iter1 != weights.end(); ++iter1, ++iter2) {
			if (iter1->first != iter2->first || iter1->second != iter2->second) {
				return false;
			}
		}

		return true;
	}
};

struct LoopDetector {
	LoopDetectorData temp;
	std::multimap<EdgeId, LoopDetectorData> data;

	void AddNewEdge(EdgeId e, size_t iter, double weight = 0) {
		temp.SetSelectedEdge(iter, weight);
		LoopDetectorData
		data.insert(std::make_pair(e, ));
	}

	void clear() {
		data.clear();
		temp.clear();
	}

	void print() {
		for (auto iter = data.begin(); iter != data.end(); ++iter) {
			INFO("Edge " << iter->first << ", weight " << iter->second.selfWeight << ", iteration " << iter->second.iteration);
			for(auto alt = iter->second.weights.begin(); alt != iter->second.weights.end(); ++alt) {
				INFO("Edge " << alt->first << ", weight " << alt->second);
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

size_t CountLoopLength(Graph& g, BidirectionalPath& path, LoopDetector& detector, bool forward) {
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
size_t CountEdgesToRemove(EdgeId lastEdge, LoopDetector& detector, bool fullRemoval, size_t loopCount) {
	size_t loopSize = CountLoopEdges(lastEdge, detector);

	if (fullRemoval) {
		return loopCount * loopSize;
	} else {
		return (loopCount - 1) * loopSize + 1;
	}
}

//Cut loop forward
void RemoveLoopForward(BidirectionalPath& path, LoopDetector& detector, bool fullRemoval, size_t loopCount) {
	size_t edgesToRemove = CountEdgesToRemove(path.back(), detector, fullRemoval, loopCount);

	for(size_t i = 0; i < edgesToRemove; ++i) {
		path.pop_back();
	}
}

void RemoveLoopBackward(BidirectionalPath& path, LoopDetector& detector, bool fullRemoval, size_t loopCount) {
	size_t edgesToRemove = CountEdgesToRemove(path.front(), detector, fullRemoval, loopCount);

	for(size_t i = 0; i < edgesToRemove; ++i) {
		path.pop_front();
	}
}

bool LoopBecameStable(EdgeId e, LoopDetector& detector) {
	auto iter = detector.data.upper_bound(e);
	auto last = --iter;
	auto prev = --iter;

	return prev->second == last->second;
}

size_t CountLoopExits(BidirectionalPath& path, EdgeId e, LoopDetector& detector) {
	size_t loopSize = CountLoopEdges(e, detector);
	size_t exits = 0;

	for (int i = (int) (path.size() - 1); i >= (int) (path.size() - loopSize); --i) {
		LoopDetectorData& data = detector.data.find(path[i])->second;

		exits += data.weights.size() - 1;
	}
	return exits;
}

EdgeId FindFirstFork(BidirectionalPath& path, EdgeId e, LoopDetector& detector) {
	size_t loopSize = CountLoopEdges(e, detector);

	for (int i = (int) (path.size() - 1); i >= (int) (path.size() - loopSize); --i) {
		LoopDetectorData& data = detector.data.find(path[i])->second;

		if (data.weights.size() == 2) {
			return path[i];
		}
	}
	return 0;
}

EdgeId GetForwardFork(Graph& g, EdgeId e) {
	VertexId v = g.EdgeStart(e);
	if (g.OutgoingEdgeCount(v) != 2) {
		return 0;
	}
	std::vector<EdgeId> edges = g.OutgoingEdges(v);
	if (edges[1] == e) {
		return edges[0];
	} else {
		return edges[1];
	}
}

EdgeId GetBackwardFork(Graph& g, EdgeId e) {
	VertexId v = g.EdgeEnd(e);
	if (g.IncomingEdgeCount(v) != 2) {
		return 0;
	}
	std::vector<EdgeId> edges = g.IncomingEdges(v);
	if (edges[1] == e) {
		return edges[0];
	} else {
		return edges[1];
	}
}

bool EdgesMakeShortLoop(Graph& g, EdgeId e1, EdgeId e2) {
	return g.EdgeStart(e1) == g.EdgeEnd(e2) && g.EdgeStart(e2) == g.EdgeEnd(e1);
}

EdgeId IsEdgeInShortLoopForward(Graph& g, EdgeId e) {
	VertexId v = g.EdgeEnd(e);
	auto edges = g.OutgoingEdges(v);
	EdgeId result = 0;

	for (auto edge = edges.begin(); edge != edges.end(); ++edge) {
		if (g.EdgeEnd(*edge) == g.EdgeStart(e)) {
			result = *edge;
		}
	}

	if (result != 0 && g.OutgoingEdgeCount(v) == 1) {
		result = e;
	}

	return result;
}

EdgeId IsEdgeInShortLoopBackward(Graph& g, EdgeId e) {
	VertexId v = g.EdgeStart(e);
	auto edges = g.IncomingEdges(v);
	EdgeId result = 0;

	for (auto edge = edges.begin(); edge != edges.end(); ++edge) {
		if (g.EdgeStart(*edge) == g.EdgeEnd(e)) {
			result = *edge;
		}
	}

	if (result != 0 && g.IncomingEdgeCount(v) == 1) {
		result = e;
	}

	return result;
}

size_t GetMaxExitIteration(EdgeId loopEdge, EdgeId loopExit, LoopDetector& detector) {
	auto range = detector.data.equal_range(loopEdge);

	size_t maxIter = 0;
	double maxWeight = 0;
	for (auto iter = range.first; iter != range.second; ++iter) {
		if (iter->second.weights[loopExit] > maxWeight) {
			maxIter = iter->second.iteration;
		}
	}
	return maxIter;
}

size_t GetFirstExitIteration(EdgeId loopEdge, EdgeId loopExit, LoopDetector& detector, double coeff = lc_cfg::get().es.priority_coeff) {
	auto range = detector.data.equal_range(loopEdge);

	size_t maxIter = std::numeric_limits<size_t>::max();
	for (auto iter = range.first; iter != range.second; ++iter) {
		if (iter->second.weights[loopExit] * coeff > iter->second.weights[loopExit] && maxIter > iter->second.iteration) {
			maxIter = iter->second.iteration;
		}
	}
	return maxIter;
}


} //namespace long_contigs


#endif /* LOOP_HPP_ */

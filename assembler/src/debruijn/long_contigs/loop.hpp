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

	std::vector<EdgeId> alternatives;
	std::vector<double> weights;

	LoopDetectorData(size_t iter, double weight): iteration(iter), selfWeight(weight), alternatives(), weights()  {
	}

	LoopDetectorData(): alternatives(), weights()  {
	}

	void SetSelectedEdge(size_t iter, double w) {
		iteration = iter;
		selfWeight = w;
	}

	void AddAlternative(EdgeId e, double w) {
		alternatives.push_back(e);
		weights.push_back(w);
	}

	void clear() {
		alternatives.clear();
		weights.clear();
		iteration = 0;
		selfWeight = 0;
	}

	bool operator==(const LoopDetectorData& d) {
		if (selfWeight != d.selfWeight || weights.size() != d.weights.size()) {
			return false;
		}

		for (size_t i = 0; i != weights.size(); ++i) {
			if (alternatives[i] != d.alternatives[i] ||  weights[i] != d.weights[i]) {
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
		temp.SetSelectedEdge(iter, w);
		data.insert(std::make_pair(e, temp));
	}

	void clear() {
		data.clear();
		temp.clear();
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
	EdgeId lastEdge = fowrard ? path.back() : path.front();
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
	size_t edgesToRemove = CountEdgesToRemove(path.back(), detector, fullRemoval);

	for(size_t i = 0; i < edgesToRemove; ++i) {
		path.pop_back();
	}
}

void RemoveLoopBackward(BidirectionalPath& path, LoopDetector& detector, bool fullRemoval, size_t loopCount) {
	size_t edgesToRemove = CountEdgesToRemove(path.front(), detector, fullRemoval);

	for(size_t i = 0; i < edgesToRemove; ++i) {
		path.pop_front();
	}
}

bool LoopBecameStable(EdgeId e, LoopDetector& detector, size_t maxImitateCount) {
	auto iter = detector.data.upper_bound(lastEdge);
	auto last = --iter;
	auto prev = --iter;

	return *prev == *last;
}

size_t GetMaxExitIteration(EdgeId e, LoopDetector& detector);

size_t GetFirstExitIteration();


} //namespace long_contigs


#endif /* LOOP_HPP_ */

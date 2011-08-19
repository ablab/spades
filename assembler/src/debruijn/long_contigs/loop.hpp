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
};

struct LoopDetector {
	LoopDetectorData temp;
	std::multimap<EdgeId, LoopDetectorData> detector;

	void AddNewEdge(EdgeId e, size_t iteration, double weight) {
		detector.insert(std::make_pair(e, LoopDetectorData(iteration, weight)));
	}

	void AddNewEdge(EdgeId e) {
		detector.insert(std::make_pair(e, temp));
	}

};


////Add edge to cycle detector and check
//bool CheckCycle(BidirectionalPath& path, EdgeId extension, CycleDetector& detector, double weight) {
//	static size_t MAX_LOOPS = lc_cfg::get().lr.max_loops;
//	detector.insert(std::make_pair(extension, std::make_pair(path.size(), weight)));
//
//	return detector.count(extension) > MAX_LOOPS;
//}
//
////Edges to remove
//size_t CountLoopEdges(EdgeId lastEdge, CycleDetector& detector, bool fullRemoval) {
//	static size_t MAX_LOOPS = lc_cfg::get().lr.max_loops;
//	auto iter = detector.upper_bound(lastEdge);
//
//	--iter;
//	size_t loopSize = iter->second.first;
//	--iter;
//	loopSize -= iter->second.first;
//
//	if (fullRemoval) {
//		return MAX_LOOPS * loopSize;
//	} else {
//		return (MAX_LOOPS - 1) * loopSize + 1;
//	}
//}
//
////Cut loop forward
//void RemoveLoopForward(BidirectionalPath& path, CycleDetector& detector, bool fullRemoval) {
//	size_t edgesToRemove = CountLoopEdges(path.back(), detector, fullRemoval);
//	//INFO("Removing loop of " << edgesToRemove << " edges");
//
//	for(size_t i = 0; i < edgesToRemove; ++i) {
//		path.pop_back();
//	}
//}
//
//void RemoveLoopBackward(BidirectionalPath& path, CycleDetector& detector, bool fullRemoval) {
//	size_t edgesToRemove = CountLoopEdges(path.front(), detector, fullRemoval);
//	//INFO("Removing loop of " << edgesToRemove << " edges");
//
//	for(size_t i = 0; i < edgesToRemove; ++i) {
//		path.pop_front();
//	}
//}

bool LoopBecameStable();

size_t GetExitIteration();

size_t ImitateLoop();

} //namespace long_contigs


#endif /* LOOP_HPP_ */

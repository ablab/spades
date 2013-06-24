/*
 * loop_traverser.hpp
 *
 *  Created on: Jan 28, 2013
 *      Author: ira
 */

#ifndef LOOP_TRAVERSER_H_
#define LOOP_TRAVERSER_H_

#include "path_extender.hpp"
#include "pe_resolver.hpp"
#include "path_visualizer.hpp"

namespace path_extend {

class LoopTraverser {

	Graph& g_;
	PathContainer& paths_;
	GraphCoverageMap& covMap_;
	PathExtender* extender_;
private:
	EdgeId FindStart(const set<VertexId>& component_set) const{
		EdgeId result;
		for (auto it = component_set.begin(); it != component_set.end(); ++it) {
			for (auto eit = g_.in_begin(*it); eit != g_.in_end(*it); ++eit) {
				if (component_set.count(g_.EdgeStart(*eit)) == 0) {
					if (result != EdgeId()) {
						return EdgeId();
					}
					result = *eit;
				}
			}
		}
		return result;
	}

	EdgeId FindFinish(const set<VertexId>& component_set) {
		EdgeId result;
		for (auto it = component_set.begin(); it != component_set.end(); ++it) {
			for (auto I = g_.out_begin(*it), E = g_.out_end(*it);
					I != E; ++I) {
				if (component_set.count(g_.EdgeEnd(*I)) == 0) {
					if (result != EdgeId()) {
						return EdgeId();
					}
					result = *I;
				}
			}
		}
		return result;
	}

	void TraverseLoop(EdgeId start, EdgeId end, set<VertexId>& component, PathContainer paths) {
		std::set<BidirectionalPath*> coveredStartPaths =
				covMap_.GetCoveringPaths(start);
		std::set<BidirectionalPath*> coveredEndPaths =
				covMap_.GetCoveringPaths(end);
		for (auto startPath = coveredStartPaths.begin();
				startPath != coveredStartPaths.end(); ++startPath) {
			if ((*startPath)->FindAll(end).size() > 0) {
				return;
			}
		}
		if (coveredStartPaths.size() < 1 or coveredEndPaths.size() < 1) {
			DEBUG(
					"TraverseLoop STRANGE SITUATION: start " << coveredStartPaths.size() << " end " << coveredEndPaths.size());
			return;
		}
		BidirectionalPath* startPath = *coveredStartPaths.begin();
		BidirectionalPath* endPath = *coveredEndPaths.begin();
		if ((*startPath) == endPath->Conjugate()){
			return;
		}
		extender_->GrowPath(*startPath);
		extender_->GrowPath(*endPath->getConjPath());

		size_t commonSize = startPath->CommonEndSize(*endPath);
		size_t nLen = 0;
		if (commonSize == 0) {
			DijkstraSearcher pathSeacher(g_);
			VertexId lastVertex = g_.EdgeEnd(
					startPath->At(startPath->Size() - 1));
			VertexId firstVertex = g_.EdgeStart(endPath->At(0));
			vector<EdgeId> pathToAdd = pathSeacher.FindShortestPathsFrom(
					lastVertex)[firstVertex];
			for (size_t i = 0; i < pathToAdd.size(); ++i){
				nLen += g_.length(pathToAdd[i]);
			}
			nLen += g_.k();
			/*for (size_t i = 0; i < pathToAdd.size(); ++i) {
				startPath->PushBack(pathToAdd[i]);
			}*/
		}
		if (commonSize < endPath->Size()){
			startPath->PushBack(endPath->At(commonSize), nLen);
		}
		for (size_t i = commonSize + 1; i < endPath->Size(); ++i) {
			startPath->PushBack(endPath->At(i));
		}

		startPath->clearOverlapedEnd();
		endPath->Clear();
	}

public:
	LoopTraverser(Graph& g, PathContainer& paths, GraphCoverageMap& coverageMap, PathExtender* extender) :
			g_(g), paths_(paths), covMap_(coverageMap), extender_(extender) {
	}

	void TraverseAllLoops() {
		DEBUG("TraverseAllLoops");
		LongEdgesExclusiveSplitter<Graph> splitter(g_, 1000);
		while (!splitter.Finished()) {
			vector<VertexId> component = splitter.NextComponent();
			if (component.size() > 10)
				continue;
			set<VertexId> component_set(component.begin(), component.end());
			EdgeId start = FindStart(component_set);
			EdgeId finish = FindFinish(component_set);
			if (start == EdgeId() || finish == EdgeId()) {
				continue;
			}
			TraverseLoop(start, finish, component_set, paths_);
		}

	}
};

}

#endif /* LOOP_TRAVERSER_H_ */

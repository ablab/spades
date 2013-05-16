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
	Graph& graph_;
	PathContainer& paths_;
	GraphCoverageMap& coverageMap_;
	PathExtender* extender_;
private:
	EdgeId FindStart(const set<VertexId>& component_set) const{
		EdgeId result;
		for (auto it = component_set.begin(); it != component_set.end(); ++it) {
			for (auto eit = graph_.in_begin(*it); eit != graph_.in_end(*it); ++eit) {
				if (component_set.count(graph_.EdgeStart(*eit)) == 0) {
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
			for (auto I = graph_.out_begin(*it), E = graph_.out_end(*it);
					I != E; ++I) {
				if (component_set.count(graph_.EdgeEnd(*I)) == 0) {
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
				coverageMap_.GetCoveringPaths(start);
		std::set<BidirectionalPath*> coveredEndPaths =
				coverageMap_.GetCoveringPaths(end);
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
		if (startPath->hasOverlapedEnd()) {
			BidirectionalPath* overlap = *startPath->getOverlapedEnd().begin();
			DEBUG("overlap");
			if ((*overlap) != startPath->Conjugate()) {
				for (size_t i = 0; i < overlap->Size(); ++i) {
					startPath->PushBack(overlap->At(i));
				}
				if ((overlap->FindAll(end)).size() != 0) {
					overlap->Clear();
					return;
				}
				overlap->Clear();
			}
		}
		extender_->GrowPath(*startPath);
		extender_->GrowPath(*endPath->getConjPath());

		size_t commonSize = startPath->CommonEndSize(*endPath);
		if (commonSize == 0) {
			DijkstraSearcher pathSeacher(graph_);
			VertexId lastVertex = graph_.EdgeEnd(
					startPath->At(startPath->Size() - 1));
			VertexId firstVertex = graph_.EdgeStart(endPath->At(0));
			vector<EdgeId> pathToAdd = pathSeacher.FindShortestPathsFrom(
					lastVertex)[firstVertex];
			for (size_t i = 0; i < pathToAdd.size(); ++i) {
				startPath->PushBack(pathToAdd[i]);
			}
		}

		for (size_t i = commonSize; i < endPath->Size(); ++i) {
			startPath->PushBack(endPath->At(i));
		}

		startPath->clearOverlapedEnd();
		startPath->addOverlapedEnd(endPath->getOverlapedEnd());
		endPath->Clear();
	}

public:
	LoopTraverser(Graph& g, PathContainer& paths, GraphCoverageMap& coverageMap, PathExtender* extender) :
			graph_(g), paths_(paths), coverageMap_(coverageMap), extender_(extender) {

	}

	void TraverseAllLoops() {
		DEBUG("TraverseAllLoops");
		LongEdgesExclusiveSplitter<Graph> splitter(graph_, 500);
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

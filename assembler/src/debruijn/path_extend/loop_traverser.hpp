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
	GraphCoverageMap& covMap_;
	ContigsMaker* extender_;
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

	void TraverseLoop(EdgeId start, EdgeId end) {
	    DEBUG("start " << g_.int_id(start) << " end " << g_.int_id(end));
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
			DEBUG("TraverseLoop STRANGE SITUATION: start " << coveredStartPaths.size() << " end " << coveredEndPaths.size());
			return;
		}
		BidirectionalPath* startPath = *coveredStartPaths.begin();
		BidirectionalPath* endPath = *coveredEndPaths.begin();
		if ((*startPath) == endPath->Conjugate()){
			return;
		}
		extender_->GrowPath(*startPath);
		extender_->GrowPath(*endPath->GetConjPath());

		size_t commonSize = startPath->CommonEndSize(*endPath);
		size_t nLen = 0;
		if (commonSize == 0) {
		    //TODO: use another seacher(omnigraph::PathProcessor), delete this DijkstraSearcher
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
			startPath->PushBack(endPath->At(commonSize), (int) nLen);
		}
		for (size_t i = commonSize + 1; i < endPath->Size(); ++i) {
			startPath->PushBack(endPath->At(i));
		}
		DEBUG("travers");
		startPath->Print();
		endPath->Print();
		DEBUG("conj");
		endPath->GetConjPath()->Print();
		DEBUG("hear1");
		endPath->Clear();
		DEBUG("hear2");
	}

public:
	LoopTraverser(Graph& g, GraphCoverageMap& coverageMap, ContigsMaker* extender) :
			g_(g), covMap_(coverageMap), extender_(extender) {
	}

	void TraverseAllLoops() {
		DEBUG("TraverseAllLoops");
		shared_ptr<GraphSplitter<Graph>> splitter = LongEdgesExclusiveSplitter<Graph>(g_, 1000);
		while (splitter->HasNext()) {
			GraphComponent<Graph> component = splitter->Next();
			if (component.v_size() > 10)
				continue;
			set<VertexId> component_set(component.v_begin(), component.v_end());
			EdgeId start = FindStart(component_set);
			EdgeId finish = FindFinish(component_set);
			if (start == EdgeId() || finish == EdgeId()) {
				continue;
			}
			TraverseLoop(start, finish);
		}

	}
protected:
    DECL_LOGGER("LoopTraverser");
};

}

#endif /* LOOP_TRAVERSER_H_ */

//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

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

	const Graph& g_;
	GraphCoverageMap& covMap_;
	shared_ptr<ContigsMaker> extender_;
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

	void TryToGrow(BidirectionalPath* path, EdgeId component_entrance) {
		BidirectionalPath clone = *path;
		extender_->GrowPath(*path);
		if (!path->Contains(component_entrance)) {
			DEBUG("Grown paths do not contain initial edges, rolling back");
			path->Clear();
			path->PushBack(clone);
		}
	}

    bool IsEndInsideComponent(const BidirectionalPath &path,
                              const set <VertexId> &component_set) {
        if (component_set.count(g_.EdgeStart(path.Front())) == 0) {
            return false;
        }
        for (size_t i = 0; i < path.Size(); ++i) {
            if (component_set.count(g_.EdgeEnd(path.At(i))) == 0)
                return false;
        }
        return true;
    }


    bool IsEndInsideComponent(const BidirectionalPath &path, EdgeId component_entrance,
                              const set <VertexId> &component_set,
                              bool conjugate = false) {
        int i = path.FindLast(component_entrance);
        VERIFY_MSG(i != -1, "Component edge is not found in the path")

        if ((size_t) i == path.Size() - 1) {
            if (conjugate)
                return component_set.count(g_.conjugate(g_.EdgeEnd(path.Back()))) > 0;
            else
                return component_set.count(g_.EdgeEnd(path.Back())) > 0;
        }

        if (conjugate)
            return IsEndInsideComponent(path.SubPath((size_t) i + 1).Conjugate(), component_set);
        else
            return IsEndInsideComponent(path.SubPath((size_t) i + 1), component_set);
    }

	void TraverseLoop(EdgeId start, EdgeId end, const set<VertexId>& component_set) {
	    DEBUG("start " << g_.int_id(start) << " end " << g_.int_id(end));
	    BidirectionalPathSet coveredStartPaths =
				covMap_.GetCoveringPaths(start);
	    BidirectionalPathSet coveredEndPaths =
				covMap_.GetCoveringPaths(end);

		for (auto it_path = coveredStartPaths.begin();
				it_path != coveredStartPaths.end(); ++it_path) {
			if ((*it_path)->FindAll(end).size() > 0) {
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

		TryToGrow(startPath, start);
		TryToGrow(endPath->GetConjPath(), g_.conjugate(end));

		//Checking that paths ends are within component
        if (!IsEndInsideComponent(*startPath, start, component_set) ||
                !IsEndInsideComponent(*endPath->GetConjPath(), g_.conjugate(end), component_set, true)) {
            DEBUG("Some path goes outside of the component")
            return;
        }

		size_t commonSize = startPath->CommonEndSize(*endPath);
		size_t nLen = 0;
        DEBUG("Str " << startPath->Size() << ", end" << endPath->Size());
        if (commonSize == 0 && !startPath->Empty() > 0 && !endPath->Empty()) {
            DEBUG("Estimating gap size");
            VertexId lastVertex = g_.EdgeEnd(startPath->Back());
            VertexId firstVertex = g_.EdgeStart(endPath->Front());

            if (firstVertex == lastVertex) {
                nLen = 0;
            } else {
                DijkstraHelper<Graph>::BoundedDijkstra dijkstra(DijkstraHelper<Graph>::CreateBoundedDijkstra(g_, 1000, 3000));
                dijkstra.run(lastVertex);
                vector<EdgeId> shortest_path = dijkstra.GetShortestPathTo(g_.EdgeStart(endPath->Front()));

                if (shortest_path.size() == 0) {
                    DEBUG("Failed to find closing path");
                    return;
                } else if (!IsEndInsideComponent(BidirectionalPath(g_, shortest_path), component_set)) {
                    DEBUG("Closing path is outside the component");
                    return;
                } else {
                    for (size_t i = 0; i < shortest_path.size(); ++i) {
                        nLen += g_.length(shortest_path[i]);
                    }
                    nLen += g_.k();
                }
            }
        }
		if (commonSize < endPath->Size()){
			startPath->PushBack(endPath->At(commonSize), (int) nLen);
		}
		for (size_t i = commonSize + 1; i < endPath->Size(); ++i) {
            startPath->PushBack(endPath->At(i), endPath->GapAt(i));
		}
		DEBUG("travers");
		startPath->Print();
		endPath->Print();
		DEBUG("conj");
		endPath->GetConjPath()->Print();
		endPath->Clear();
	}

public:
	LoopTraverser(const Graph& g, GraphCoverageMap& coverageMap, shared_ptr<ContigsMaker> extender) :
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
			TraverseLoop(start, finish, component_set);
		}

	}
protected:
    DECL_LOGGER("LoopTraverser");
};

}

#endif /* LOOP_TRAVERSER_H_ */

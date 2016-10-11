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
    size_t long_edge_limit_;
    size_t component_size_limit_;
    size_t shortest_path_limit_;
    static const size_t DIJKSTRA_LIMIT = 3000;
private:
    bool AnyTipsInComponent(const GraphComponent<Graph>& component) const{
        for(auto e : component.edges()) {
            if (g_.IncomingEdgeCount(g_.EdgeStart(e)) == 0 || g_.OutgoingEdgeCount(g_.EdgeEnd(e)) == 0)
                return true;
        }
        return false;
    }

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

    bool TraverseLoop(EdgeId start, EdgeId end, const set<VertexId>& component_set) {
        DEBUG("start " << g_.int_id(start) << " end " << g_.int_id(end));
        BidirectionalPathSet coveredStartPaths =
                covMap_.GetCoveringPaths(start);
        BidirectionalPathSet coveredEndPaths =
                covMap_.GetCoveringPaths(end);

        for (auto it_path = coveredStartPaths.begin();
                it_path != coveredStartPaths.end(); ++it_path) {
            if ((*it_path)->FindAll(end).size() > 0) {
                return false;
            }
        }
        if (coveredStartPaths.size() < 1 or coveredEndPaths.size() < 1) {
            DEBUG("TraverseLoop STRANGE SITUATION: start " << coveredStartPaths.size() << " end " << coveredEndPaths.size());
            return false;
        }

        if (coveredStartPaths.size() > 1 or coveredEndPaths.size() > 1) {
            DEBUG("Ambiguous situation in path joining, quitting");
            return false;
        }

        BidirectionalPath* startPath = *coveredStartPaths.begin();
        BidirectionalPath* endPath = *coveredEndPaths.begin();
        if ((*startPath) == endPath->Conjugate()){
            return false;
        }

        //Checking that paths ends are within component
        if (!IsEndInsideComponent(*startPath, start, component_set) ||
                !IsEndInsideComponent(*endPath->GetConjPath(), g_.conjugate(end), component_set, true)) {
            DEBUG("Some path goes outside of the component")
            return false;
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
                DijkstraHelper<Graph>::BoundedDijkstra dijkstra(DijkstraHelper<Graph>::CreateBoundedDijkstra(g_, shortest_path_limit_,
                                                                                                             DIJKSTRA_LIMIT));
                dijkstra.Run(lastVertex);
                vector<EdgeId> shortest_path = dijkstra.GetShortestPathTo(g_.EdgeStart(endPath->Front()));

                if (shortest_path.empty()) {
                    DEBUG("Failed to find closing path");
                    return false;
                } else if (!IsEndInsideComponent(BidirectionalPath(g_, shortest_path), component_set)) {
                    DEBUG("Closing path is outside the component");
                    return false;
                } else {
                    nLen = CumulativeLength(g_, shortest_path);
                }
            }
        }
        if (commonSize < endPath->Size()){
            startPath->PushBack(endPath->At(commonSize), (int) nLen);
        }
        for (size_t i = commonSize + 1; i < endPath->Size(); ++i) {
            startPath->PushBack(endPath->At(i), endPath->GapAt(i), endPath->TrashPreviousAt(i), endPath->TrashCurrentAt(i));
        }
        DEBUG("travers");
        startPath->Print();
        endPath->Print();
        DEBUG("conj");
        endPath->GetConjPath()->Print();
        endPath->Clear();
        return true;
    }

    bool ContainsLongEdges(const GraphComponent<Graph>& component) const {
        for(auto e : component.edges()) {
            if(g_.length(e) > long_edge_limit_) {
                return true;
            }
        }
        return false;
    }

public:
    LoopTraverser(const Graph& g, GraphCoverageMap& coverageMap, size_t long_edge_limit, size_t component_size_limit, size_t shortest_path_limit) :
            g_(g), covMap_(coverageMap), long_edge_limit_(long_edge_limit), component_size_limit_(component_size_limit), shortest_path_limit_(shortest_path_limit) {
    }

    size_t TraverseAllLoops() {
        DEBUG("TraverseAllLoops");
        size_t traversed = 0;
        shared_ptr<GraphSplitter<Graph>> splitter = LongEdgesExclusiveSplitter<Graph>(g_, long_edge_limit_);
        while (splitter->HasNext()) {
            GraphComponent<Graph> component = splitter->Next();
            if (component.v_size() > component_size_limit_)
                continue;
            if (ContainsLongEdges(component))
                continue;
            if (AnyTipsInComponent(component))
                continue;

            set<VertexId> component_set(component.v_begin(), component.v_end());
            EdgeId start = FindStart(component_set);
            EdgeId finish = FindFinish(component_set);
            if (start == EdgeId() || finish == EdgeId()) {
                continue;
            }
            if (TraverseLoop(start, finish, component_set))
                ++traversed;
        }
        return traversed;
    }

protected:
    DECL_LOGGER("LoopTraverser");
};

}

#endif /* LOOP_TRAVERSER_H_ */

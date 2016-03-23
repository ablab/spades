//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "../utils/bulge_utils.hpp"
#include "../dipspades_config.hpp"

using namespace debruijn_graph;

namespace dipspades {

class BaseBulgeCorrectionCondition{
protected:
    Graph &graph_;
public:
    BaseBulgeCorrectionCondition(Graph &graph) : graph_(graph) { }
    virtual bool IsBulgeCorrect(shared_ptr<BaseBulge> bulge) = 0;
    virtual bool IsPathBulgeSide(vector<EdgeId> path) = 0;
    virtual ~BaseBulgeCorrectionCondition(){ }
};

class RelatedVerticesCondition : public BaseBulgeCorrectionCondition {
    bool TwoVerticesRelated(VertexId v1, VertexId v2){
        return graph_.RelatedVertices(v1, v2);
    }

    bool PathContainsNoRelatedToVertex(vector<EdgeId> path, VertexId vertex,
            bool check_start_vertex = false){
        VERIFY(path.size() != 0);

        if(check_start_vertex)
            if(TwoVerticesRelated(graph_.EdgeStart(path[0]), vertex))
                return false;

        for(auto e = path.begin(); e != path.end(); e++)
            if(TwoVerticesRelated(vertex, graph_.EdgeEnd(*e)))
                return false;
        return true;
    }

    bool PathContainsNoRelatedVertices(vector<EdgeId> path){
        if(!PathContainsNoRelatedToVertex(path, graph_.EdgeStart(path[0])))
            return false;
        for(auto e1 = path.begin(); e1 != path.end(); e1++)
            for(auto e2 = e1 + 1; e2 != path.end(); e2++)
                if(TwoVerticesRelated(graph_.EdgeEnd(*e1), graph_.EdgeEnd(*e2)))
                    return false;
        return true;
    }

    bool PathsContainNoRelatedVertices(shared_ptr<BaseBulge> bulge){
        auto path1 = bulge->path1();
        auto path2 = bulge->path2();
        for(auto e1 = path1.begin(); e1 != path1.end(); e1++)
            for(auto e2 = path2.begin(); e2 != path2.end(); e2++)
                if((e1 != path1.end() - 1) && (e2 != path2.end() - 1))
                    if(TwoVerticesRelated(graph_.EdgeEnd(*e1), graph_.EdgeEnd(*e2)))
                        return false;
        return true;
    }

public:
    RelatedVerticesCondition(Graph &graph) : BaseBulgeCorrectionCondition(graph) { }
    bool IsBulgeCorrect(shared_ptr<BaseBulge> bulge){
        if(!PathContainsNoRelatedVertices(bulge->path1()) ||
                !PathContainsNoRelatedVertices(bulge->path2()))
            return false;
        return PathsContainNoRelatedVertices(bulge);
    }

    bool IsPathBulgeSide(vector<EdgeId> path){
        return PathContainsNoRelatedVertices(path);
    }
};

class AdjacencyToAutoRCEdges : public BaseBulgeCorrectionCondition {

public:
    AdjacencyToAutoRCEdges(Graph &graph) : BaseBulgeCorrectionCondition(graph) { }

    bool IsBulgeCorrect(shared_ptr<BaseBulge> bulge){
        return IsPathBulgeSide(bulge->path1()) && IsPathBulgeSide(bulge->path2());
    }

    bool IsPathBulgeSide(vector<EdgeId> path){
        return !PathAdjacentRelatedEdges(graph_, path);
    }
};

class DiploidyCondition : public BaseBulgeCorrectionCondition {
    double rel_length_;
    double rel_align_;
public:
    DiploidyCondition(Graph &graph,
            double rel_length,
            double rel_align) :
        BaseBulgeCorrectionCondition(graph),
        rel_length_(rel_length),
        rel_align_(rel_align) { }

    bool IsBulgeCorrect(shared_ptr<BaseBulge> bulge){
        return bulge->IsBulgeDiploid(rel_length_, rel_align_);
    }

    bool IsPathBulgeSide(vector<EdgeId>){
        return true;
    }
};

class CorrectSplitCondition : public BaseBulgeCorrectionCondition {
public:
    CorrectSplitCondition(Graph &graph) : BaseBulgeCorrectionCondition(graph) { }

    bool IsBulgeCorrect(shared_ptr<BaseBulge> bulge){
        return bulge->path1().size() == bulge->path2().size();
    }

    bool IsPathBulgeSide(vector<EdgeId>){
        return true;
    }
};

}

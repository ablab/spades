//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "../utils/bulge_utils.hpp"
#include "../utils/path_routines.hpp"
#include "../dipspades_config.hpp"

using namespace debruijn_graph;

namespace dipspades {


class GluingDirectionDefiner {
protected:
    Graph &graph_;
public:
    GluingDirectionDefiner(Graph &graph) : graph_(graph) { }
    virtual glue_direction Define(shared_ptr<BaseBulge>) {
        return undefined;
    }
    virtual ~GluingDirectionDefiner() { }
};

class RelatedBaseGlueDirectionDefiner : public GluingDirectionDefiner{
public:
    RelatedBaseGlueDirectionDefiner(Graph &graph) : GluingDirectionDefiner(graph) { }

    glue_direction Define(shared_ptr<BaseBulge> bulge){
        bool rel_edges_path1 = PathAdjacentRelatedEdges(this->graph_, bulge->path1());
        bool rel_edges_path2 = PathAdjacentRelatedEdges(this->graph_, bulge->path2());
        if(rel_edges_path1 && rel_edges_path2)
            return undefined;

        // if only path2 contains related edges
        // we need gluing path2 to path1
        if(rel_edges_path2)
            return reverse_gluing;
        return direct_gluing;
    }
};

class CoverageBaseGlueDirectionDefiner : public GluingDirectionDefiner{
public:
    CoverageBaseGlueDirectionDefiner(Graph &graph) : GluingDirectionDefiner(graph) { }

    glue_direction Define(shared_ptr<BaseBulge>){
        // todo implement me
        return direct_gluing;
    }
};

class CompositeGlueDirectionDefiner : public GluingDirectionDefiner {
    vector<shared_ptr<GluingDirectionDefiner> > &definers_;
public:
    CompositeGlueDirectionDefiner(Graph &graph,
            vector<shared_ptr<GluingDirectionDefiner> > &definers) :
                GluingDirectionDefiner(graph),
                definers_(definers) { }

    glue_direction Define(shared_ptr<BaseBulge> bulge){
        set<glue_direction> directions;
        for(auto it = definers_.begin(); it != definers_.end(); it++)
            directions.insert((*it)->Define(bulge));
        if(directions.size() == 1)
            return *(directions.begin());
        return undefined;
    }
};

}

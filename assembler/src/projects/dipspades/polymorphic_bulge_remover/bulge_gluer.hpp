//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "glue_direction_definer.hpp"
#include "gluing_vertices_definer.hpp"
#include "bulge_splitter.hpp"
#include "bulge_correction_condition.hpp"
#include "../utils/edge_gluer.hpp"

using namespace debruijn_graph;

namespace dipspades {

template<class GlueDirectionDefiner, class GluingVericesDefiner, class BulgeSplitter>
class ComplexBulgeGluer {

    Graph &graph_;
    GlueDirectionDefiner glue_dir_definer_;
    GluingVericesDefiner glue_definer_;
    BulgeSplitter splitter_;

    bool IsSplittedBulgeCorrect(shared_ptr<BaseBulge> splitted_bulge){
        return !splitted_bulge->IsEmpty() && CorrectSplitCondition(graph_).IsBulgeCorrect(splitted_bulge) &&
                RelatedVerticesCondition(graph_).IsBulgeCorrect(splitted_bulge);
    }

    void GlueSplittedBulge(shared_ptr<BaseBulge> splitted_bulge){
        size_t bulge_edge_len = splitted_bulge->path1().size();
        EdgeGluer edge_gluer(graph_);
        TRACE("Edge gluer starts");
        for(size_t i = 0; i < bulge_edge_len - 1; i++){
            auto edge1 = splitted_bulge->path1()[i];
            auto edge2 = splitted_bulge->path2()[i];
            auto next_edge1 = splitted_bulge->path1()[i + 1];
            TRACE("edge1 - " << graph_.str(edge1) << ", edge2 - " << graph_.str(edge2) <<
                    ", next_edge1 - " << graph_.str(next_edge1));
            vector<EdgeId> tmp = {edge1, edge2, next_edge1};
            edge_gluer.MoveEdgesFromVertexToVertex(
                    graph_.EdgeEnd(edge1),
                    graph_.EdgeEnd(edge2),
                    tmp);
            graph_.GlueEdges(edge1, edge2);
            TRACE("Edges were moved");
        }
        graph_.GlueEdges(splitted_bulge->path1()[bulge_edge_len - 1],
                splitted_bulge->path2()[bulge_edge_len - 1]);
        TRACE("Gluing was completed");
    }

public:
    ComplexBulgeGluer(Graph &graph, GlueDirectionDefiner glue_dir_definer,
            GluingVericesDefiner glue_definer, BulgeSplitter splitter) :
        graph_(graph),
        glue_dir_definer_(glue_dir_definer),
        glue_definer_(glue_definer),
        splitter_(splitter) { }

    bool GlueBulge(shared_ptr<BaseBulge> bulge){
        auto glue_dir = glue_dir_definer_.Define(bulge);
        TRACE("Gluing direction - " << glue_dir);
        if(glue_dir == undefined)
            return false;

        shared_ptr<BaseBulge> directed_bulge(new DirectedBulge(graph_, bulge, glue_dir));
        TRACE("Glue vertices definer starts");
        auto glue_def_res = glue_definer_.Run(directed_bulge);
        TRACE("Bulge splitter starts");
        auto splitted_bulge = splitter_.SplitBulge(directed_bulge, glue_def_res);

        if(IsSplittedBulgeCorrect(splitted_bulge)){
            TRACE("Splitted bulge correct");
            GlueSplittedBulge(splitted_bulge);
            return true;
        }
        return false;
    }

private:
    DECL_LOGGER("ComplexBulgeGluer");
};

}

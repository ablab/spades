//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "bulge_correction_condition.hpp"
#include "bulge_gluer.hpp"
#include "../utils/histogram.hpp"

using namespace debruijn_graph;

namespace dipspades {

class SimpleBulgeRemover{
    Graph &graph_;
    BaseHistogram<size_t> &bulge_len_hist_;
    RelatedVerticesCondition rel_bulge_checker_;
    DiploidyCondition dip_bulge_checker_;
public:
    SimpleBulgeRemover(Graph &graph,
            BaseHistogram<size_t> &bulge_len_hist,
            const dipspades_config::polymorphic_br &pbr_config) :
        graph_(graph),
        bulge_len_hist_(bulge_len_hist),
        rel_bulge_checker_(graph),
        dip_bulge_checker_(graph, pbr_config.rel_bulge_length, pbr_config.rel_bulge_align) {}

    size_t Run(){
        size_t glued_edges_count = 0;
        for(auto e = graph_.SmartEdgeBegin(); !e.IsEnd(); ++e){
            vector<EdgeId> edges = graph_.GetEdgesBetween(graph_.EdgeStart(*e),
                    graph_.EdgeEnd(*e));
                if(edges.size() >= 2){
                    auto bulge = shared_ptr<BaseBulge>(new Bulge(graph_, graph_.k(), edges[0], edges[1]));
                    if(rel_bulge_checker_.IsBulgeCorrect(bulge) &&
                            dip_bulge_checker_.IsBulgeCorrect(bulge)){
                        bulge_len_hist_.Add(max<size_t>(graph_.length(edges[0]), graph_.length(edges[1])));
                        graph_.GlueEdges(edges[0], edges[1]);
                        glued_edges_count++;
                    }
                }
        }
        return glued_edges_count;
    }
};

}

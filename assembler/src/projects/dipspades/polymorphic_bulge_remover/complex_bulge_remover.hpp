//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pipeline/config_struct.hpp"
#include "pipeline/graphio.hpp"
#include "stages/construction.hpp"

#include "utils/path_routines.hpp"
#include "utils/element_printers.hpp"
#include "utils/histogram.hpp"

#include "bulge_correction_condition.hpp"
#include "bulge_gluer.hpp"
#include "diploid_bulge_finder.hpp"

#include "io/reads/splitting_wrapper.hpp"

#include <stdlib.h>
#include <memory.h>
#include <math.h>

namespace dipspades {

bool EdgeExists(Graph &graph_, size_t edge_id){
    for(auto e = graph_.SmartEdgeBegin(); !e.IsEnd(); ++e)
        if(graph_.int_id(*e) == edge_id)
            return true;
    return false;
}

template<class BulgePathsSearcher, class BulgeGluer>
class BulgeRemoverAlgorithm{
    typedef vector<vector<EdgeId> > paths;
protected:
    Graph &graph_;
    BulgeGluer bulge_gluer_;
    BaseHistogram<size_t> &hist_;
    const dipspades_config::polymorphic_br &pbr_config_;

    DiploidBulgeFinder bulge_finder_;
    DiploidyCondition dip_bulge_checker_;
    RelatedVerticesCondition rel_bulge_checker_;

    bool BulgeExistFrom(VertexId start){
        return graph_.OutgoingEdgeCount(start) > 1;
    }

    bool BulgeExistTo(VertexId end){
        return graph_.IncomingEdgeCount(end) > 1;
    }

    void FillHistogram(shared_ptr<BaseBulge> bulge){
        hist_.Add(max<size_t>(GetPathLength(graph_, bulge->path1()),
                GetPathLength(graph_, bulge->path2())));
    }

    bool FindGlueBulge(paths &bulge_paths){
        TRACE("Bulge finder from " << bulge_paths.size() << " paths starts");
        auto bulge = bulge_finder_.Find(bulge_paths);
        if(bulge->IsEmpty()){
            TRACE("Paths do not form a bulge");
            return false;
        }
        TRACE("Paths form a bulge");
        TRACE("Bulge gluing starts");
        if(!rel_bulge_checker_.IsBulgeCorrect(bulge)/* ||
                !dip_bulge_checker_.IsBulgeCorrect(bulge)*/){
            TRACE("Bulge do not successed diploid condition");
            return false;
        }

        TRACE("Correct bulge:");
        TRACE("Path1:" << SimplePathWithVerticesToString(graph_, bulge->path1()));
        TRACE("Path2:" << SimplePathWithVerticesToString(graph_, bulge->path2()));

        FillHistogram(bulge);
        TRACE("Diploid condition was passed");
        if(!bulge_gluer_.GlueBulge(bulge))
            return false;

        TRACE("Bulge gluing ends");
        return true;
    }

public:
    BulgeRemoverAlgorithm(Graph &graph,
            BulgeGluer bulge_gluer,
            BaseHistogram<size_t> &hist,
            const dipspades_config::polymorphic_br &pbr_config) :
                graph_(graph),
                bulge_gluer_(bulge_gluer),
                hist_(hist),
                pbr_config_(pbr_config),
                bulge_finder_(graph, pbr_config.rel_bulge_length, pbr_config.rel_bulge_align),
                dip_bulge_checker_(graph, pbr_config.rel_bulge_length, pbr_config.rel_bulge_align),
                rel_bulge_checker_(graph) { }

    size_t Run(){
        size_t num_merged_paths = 0;
        BulgePathsSearcher paths_searcher(graph_,
                max<size_t>(hist_.max(), pbr_config_.max_bulge_nucls_len),
                pbr_config_.max_neigh_number);
        INFO("Maximal length of glued bulge: " << hist_.max());
        TRACE("BulgeRemoverAlgorithm starts");
        for(auto v = graph_.SmartVertexBegin(); !v.IsEnd(); ++v){
            TRACE("Processing vertex " << graph_.str(*v));
            if(BulgeExistFrom(*v)){
                auto reached_vertices = paths_searcher.VerticesReachedFrom(*v);
                TRACE("Number of neigs - " << reached_vertices.size());
                for(auto neigh = SmartSetIterator<Graph, VertexId>(graph_,
                        reached_vertices.begin(), reached_vertices.end());
                         !neigh.IsEnd(); ++neigh){
                    if(*neigh != *v && BulgeExistTo(*neigh)){
                        TRACE("Bulge can be found");
                        TRACE("Processing neigh " << graph_.str(*neigh));
                        auto bulge_paths = paths_searcher.GetAllPathsTo(*v, *neigh);

                        TRACE("Bulge paths:");
                        for(auto p = bulge_paths.begin(); p != bulge_paths.end(); p++)
                            TRACE(SimplePathWithVerticesToString(graph_, *p));

                        if(FindGlueBulge(bulge_paths)){
                            num_merged_paths++;
                            TRACE("Bulge was glued");
                            break;
                        }
                    }
                }
            }
        }
        TRACE(num_merged_paths << " bulges were glued");
        return num_merged_paths;
    }

private:
    DECL_LOGGER("PolymorphicBulgeRemover");
};

}

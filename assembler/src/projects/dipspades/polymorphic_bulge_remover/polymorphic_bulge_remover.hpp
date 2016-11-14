//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "simple_bulge_remover.hpp"
#include "complex_bulge_remover.hpp"
#include "iterative_tails_gluing.hpp"
#include "assembly_graph/stats/picture_dump.hpp"

#include "visualization/visualization.hpp"
#include "assembly_graph/handlers/edges_position_handler.hpp"
#include "assembly_graph/components/graph_component.hpp"
#include "modules/simplification/compressor.hpp"

using namespace debruijn_graph;

namespace dipspades {

class PolymorphicBulgeRemoverHelper {
public:
    typedef ComplexBulgeGluer<RelatedBaseGlueDirectionDefiner, GluingVericesDefiner, BulgeSplitter> BaseBulgeGluer;
    static BaseBulgeGluer CreateBaseBulgeGluer(Graph &graph, double rel_len_threshold){
        return BaseBulgeGluer(graph, RelatedBaseGlueDirectionDefiner(graph),
                GluingVericesDefiner(graph, rel_len_threshold), BulgeSplitter(graph));
    }
};

class PolymorphicBulgeRemover {
    conj_graph_pack         &graph_pack_;
    BaseHistogram<size_t>     &bulge_len_hist_;

    typedef BulgeRemoverAlgorithm<DijkstraBulgePathsSearcher,
            PolymorphicBulgeRemoverHelper::BaseBulgeGluer> LightBulgeRemover;
    typedef BulgeRemoverAlgorithm<PathProcessorBulgeSearcher,
            PolymorphicBulgeRemoverHelper::BaseBulgeGluer> HardBulgeRemover;

    void RunSimpleBRCycle(){
        INFO("Simple polymorphic bulge remover runs");
        SimpleBulgeRemover spath_br(graph_pack_.g, bulge_len_hist_, dsp_cfg::get().pbr);
        size_t num_glued_bulges = 1;
        for(size_t num_iter = 1; num_glued_bulges > 0; num_iter++){
            num_glued_bulges = spath_br.Run();
            CompressAllVertices(graph_pack_.g, false);
            INFO(ToString(num_iter) + " iteration: " + ToString(num_glued_bulges) + " simple bulges were glued");
        }
        INFO("Simple polymorphic bulge remover ends");
    }

    template<class BulgeRemover>
    void BulgeRemoverCycle(string bulge_remover_name, size_t num_iters){
        INFO(bulge_remover_name + " starts");
        INFO("Maximal number of iterations: " << num_iters);
        BulgeRemover br(graph_pack_.g,
                PolymorphicBulgeRemoverHelper::CreateBaseBulgeGluer(graph_pack_.g,
                        dsp_cfg::get().pbr.paired_vert_rel_threshold),
                bulge_len_hist_,
                dsp_cfg::get().pbr);
        size_t num_glued_bulges = 1;
        for(size_t i = 0; (i < num_iters) && (num_glued_bulges != 0); i++){
            num_glued_bulges = br.Run();
            CompressAllVertices(graph_pack_.g, false);
            INFO(ToString(i + 1) + " iteration: " + ToString(num_glued_bulges) + " complex bulges were glued");
        }
        INFO(bulge_remover_name + " ends");
    }

    void WriteComponents(string component_dir) {
        if(!dsp_cfg::get().rp.developer_mode)
            return;

        graph_pack_.EnsureDebugInfo();
        make_dir(dsp_cfg::get().io.output_dir + "components/");
        visualization::graph_labeler::DefaultLabeler<Graph> labeler(graph_pack_.g, graph_pack_.edge_pos);
        make_dir(dsp_cfg::get().io.output_dir + "components/" + component_dir + "/");
        visualization::visualization_utils::WriteComponents(graph_pack_.g,
                dsp_cfg::get().io.output_dir + "components/" + component_dir + "/",
                omnigraph::ReliableSplitter<Graph>(graph_pack_.g),
                visualization::graph_colorer::DefaultColorer(graph_pack_.g, Path<EdgeId>(), Path<EdgeId>()),
                labeler);
    }

public:
    PolymorphicBulgeRemover(conj_graph_pack &graph_pack,
            BaseHistogram<size_t> &bulge_len_hist) :
        graph_pack_(graph_pack),
        bulge_len_hist_(bulge_len_hist) { }

    void Run(){
        if(!dsp_cfg::get().pbr.enabled)
            return ;
        WriteComponents("before_pbr");
        graph_pack_.kmer_mapper.SetUnsafeMode(true);
        INFO("Polymorphic bulge remover starts");
        RunSimpleBRCycle();
        BulgeRemoverCycle<LightBulgeRemover>("LightBulgeRemover", dsp_cfg::get().pbr.num_iters_lbr);
        INFO("Index refilling");
        graph_pack_.index.Refill();
        INFO("Polymorphic ends remover ends");
        WriteComponents("after_pbr");
        graph_pack_.kmer_mapper.SetUnsafeMode(false);
    }
};

}

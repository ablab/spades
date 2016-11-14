//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "utils/logger/logger.hpp"
#include "assembly_graph/stats/picture_dump.hpp"
#include "visualization/graph_labeler.hpp"
#include "paired_info/distance_estimation.hpp"
#include "paired_info/smoothing_distance_estimation.hpp"
#include "modules/path_extend/path_extend_launch.hpp"
#include "assembly_graph/graph_support/contig_output.hpp"
#include "visualization/position_filler.hpp"
#include "modules/alignment/long_read_storage.hpp"
#include "repeat_resolving.hpp"

namespace debruijn_graph {

static void PEResolving(conj_graph_pack& gp) {
    std::string scaffolds_name = cfg::get().mode == config::pipeline_type::rna ? "transcripts" : "scaffolds";
    bool output_broke_scaffolds = cfg::get().mode != config::pipeline_type::rna;

    path_extend::PathExtendParamsContainer params(cfg::get().pe_params,
                                                  cfg::get().output_dir,
                                                  "final_contigs",
                                                  scaffolds_name,
                                                  cfg::get().mode,
                                                  cfg::get().uneven_depth,
                                                  cfg::get().avoid_rc_connections,
                                                  cfg::get().use_scaffolder,
                                                  output_broke_scaffolds);



    path_extend::ResolveRepeatsPe(cfg::get().ds, params, gp);
}

static bool HasValidLibs() {
    for (const auto& lib : cfg::get().ds.reads) {
        if (!lib.is_repeat_resolvable())
            continue;

        if (!lib.is_paired() ||
            !math::eq(lib.data().mean_insert_size, 0.0)) {
            return true;
        }
    }
    
    return false;
}

void RepeatResolution::run(conj_graph_pack &gp, const char*) {
    if (cfg::get().developer_mode)
        stats::PrepareForDrawing(gp);

    visualization::graph_labeler::DefaultLabeler<Graph> labeler(gp.g, gp.edge_pos);
    stats::detail_info_printer printer(gp, labeler, cfg::get().output_dir);
    printer(config::info_printer_pos::before_repeat_resolution);

    //todo awful hack to get around PE using cfg::get everywhere...
    auto tmp_params_storage = cfg::get().pe_params;
    if (preliminary_) {
        INFO("Setting up preliminary path extend settings")
        cfg::get_writable().pe_params = *cfg::get().prelim_pe_params;
    }
    OutputContigs(gp.g, cfg::get().output_dir + "before_rr", false);
    OutputContigsToFASTG(gp.g, cfg::get().output_dir + "assembly_graph", gp.components);

    bool no_valid_libs = !HasValidLibs();

    bool use_single_reads = cfg::get().use_single_reads;
    if (cfg::get().rr_enable && no_valid_libs && !use_single_reads)
        WARN("Insert size was not estimated for any of the paired libraries, repeat resolution module will not run.");

    if ((no_valid_libs || cfg::get().rm == config::resolving_mode::none) && !use_single_reads) {
        OutputContigs(gp.g, cfg::get().output_dir + "final_contigs", false);
        return;
    }
    if (cfg::get().rm == config::resolving_mode::path_extend) {
        INFO("Using Path-Extend repeat resolving");
        PEResolving(gp);
    } else {
        INFO("Unsupported repeat resolver");
        OutputContigs(gp.g, cfg::get().output_dir + "final_contigs", false);
    }
    if (preliminary_) {
        INFO("Restoring initial path extend settings")
        cfg::get_writable().pe_params = tmp_params_storage;
    }
}

void ContigOutput::run(conj_graph_pack &gp, const char*) {
    OutputContigs(gp.g, cfg::get().output_dir + "simplified_contigs", cfg::get().use_unipaths);
    OutputContigs(gp.g, cfg::get().output_dir + "before_rr", false);
    OutputContigsToFASTG(gp.g, cfg::get().output_dir + "assembly_graph", gp.components);
    OutputContigs(gp.g, cfg::get().output_dir + "final_contigs", false);

}

} // debruijn_graph

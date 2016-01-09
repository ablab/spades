//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "standard.hpp"

#include "logger/logger.hpp"
#include "stats/debruijn_stats.hpp"
#include "omni/visualization/graph_labeler.hpp"
#include "de/distance_estimation.hpp"
#include "de/smoothing_distance_estimation.hpp"
#include "omni/omni_utils.hpp"
#include "path_extend/path_extend_launch.hpp"
#include "contig_output.hpp"
#include "positions.hpp"
#include "long_read_storage.hpp"
#include "repeat_resolving.hpp"

namespace debruijn_graph {

void PEResolving(conj_graph_pack& gp) {
    vector<size_t> indexes;
    std::string name = "scaffolds";
    bool traverse_loops = true;
    if (!(cfg::get().use_scaffolder && cfg::get().pe_params.param_set.scaffolder_options.on)) {
        name = "final_contigs";
        traverse_loops = false;
    }
    path_extend::ResolveRepeatsPe(gp, cfg::get().output_dir, name, traverse_loops, boost::optional<std::string>("final_contigs"));
}

inline bool HasValidLibs() {
    for (const auto& lib : cfg::get().ds.reads) {
        if (lib.is_repeat_resolvable()) {
            if (!lib.is_paired() || !math::eq(lib.data().mean_insert_size, 0.0)) {
                return true;
            } 
        }
    }
    return false;
}

void RepeatResolution::run(conj_graph_pack &gp, const char*) {
    if (cfg::get().developer_mode) {
        stats::PrepareForDrawing(gp);
    }

    omnigraph::DefaultLabeler<Graph> labeler(gp.g, gp.edge_pos);
    stats::detail_info_printer printer(gp, labeler, cfg::get().output_dir);
    printer(ipp_before_repeat_resolution);

    //todo awful hack to get around PE using cfg::get everywhere...
    auto tmp_params_storage = cfg::get().pe_params;
    if (preliminary_) {
        INFO("Setting up preliminary path extend settings")
        cfg::get_writable().pe_params = cfg::get().prelim_pe_params;
    }

    OutputContigs(gp.g, cfg::get().output_dir + "before_rr");
    OutputContigsToFASTG(gp.g, cfg::get().output_dir + "assembly_graph");

    bool no_valid_libs = !HasValidLibs();

    bool use_single_reads = cfg::get().use_single_reads;
    if (cfg::get().rr_enable && no_valid_libs && !use_single_reads)
        WARN("Insert size was not estimated for any of the paired libraries, repeat resolution module will not run.");

    if ((no_valid_libs || cfg::get().rm == debruijn_graph::resolving_mode::rm_none) && !use_single_reads) {
        OutputContigs(gp.g, cfg::get().output_dir + "final_contigs");
        return;
    }
    if (cfg::get().rm == debruijn_graph::resolving_mode::rm_path_extend) {
        INFO("Using Path-Extend repeat resolving");
        PEResolving(gp);
    } else {
        INFO("Unsupported repeat resolver");
        OutputContigs(gp.g, cfg::get().output_dir + "final_contigs");
    }
    if (preliminary_) {
        INFO("Restoring initial path extend settings")
        cfg::get_writable().pe_params = tmp_params_storage;
    }
}

void ContigOutput::run(conj_graph_pack &gp, const char*) {
    OutputContigs(gp.g, cfg::get().output_dir + "simplified_contigs", cfg::get().use_unipaths,
                  cfg::get().simp.tec.plausibility_length);
    OutputContigs(gp.g, cfg::get().output_dir + "before_rr");
    OutputContigsToFASTG(gp.g, cfg::get().output_dir + "assembly_graph");
    OutputContigs(gp.g, cfg::get().output_dir + "final_contigs");
}

} // debruijn_graph




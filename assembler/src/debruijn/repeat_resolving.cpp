//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

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

void RepeatResolution::run(conj_graph_pack &gp, const char*) {
    OutputContigs(gp.g, cfg::get().output_dir + "before_rr", true);
    if (cfg::get().developer_mode) {
        FillPos(gp, gp.genome, "ref0");
        FillPos(gp, !gp.genome, "ref1");
    }

    bool no_valid_libs = true;
    for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
        if (cfg::get().ds.reads[i].data().mean_insert_size != 0.0 ||
            cfg::get().ds.reads[i].is_pacbio_alignable()) {
            no_valid_libs = false;
            break;
        }
    }

    bool use_single_reads = cfg::get().use_single_reads;
    if (cfg::get().rr_enable && no_valid_libs && !use_single_reads)
        WARN("Insert size was not estimated for any of the paired libraries, repeat resolution module will not run.");

    if ((no_valid_libs || cfg::get().rm == debruijn_graph::resolving_mode::rm_none) && !use_single_reads) {
        OutputContigs(gp.g, cfg::get().output_dir + "final_contigs", true);
        return;
    }
    if (cfg::get().rm == debruijn_graph::resolving_mode::rm_path_extend) {
        INFO("Using Path-Extend repeat resolving");
        PEResolving(gp);
    } else {
        INFO("Unsupported repeat resolver");
        OutputContigs(gp.g, cfg::get().output_dir + "final_contigs", true);
    }
}

void ContigOutput::run(conj_graph_pack &gp, const char*) {
    OutputContigs(gp.g, cfg::get().output_dir + "simplified_contigs", false, cfg::get().use_unipaths,
                  cfg::get().simp.tec.plausibility_length);
    OutputContigs(gp.g, cfg::get().output_dir + "before_rr", true);
    OutputContigs(gp.g, cfg::get().output_dir + "final_contigs", true);
}


} // debruijn_graph

//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "standard.hpp"

#include "logger/logger.hpp"
#include "debruijn_stats.hpp"
#include "omni_labelers.hpp"
#include "de/distance_estimation.hpp"
#include "de/smoothing_distance_estimation.hpp"
#include "omni/omni_utils.hpp"
#include "path_extend/path_extend_launch.hpp"
#include "contig_output.hpp"
#include "long_read_storage.hpp"
#include "repeat_resolving.hpp"

namespace debruijn_graph {

//TODO: get rid of this conversion
void ConvertLongReads(LongReadContainerT& single_long_reads, vector<PathStorageInfo<Graph> > &long_reads_libs) {
    for (size_t i = 0; i < single_long_reads.size(); ++i) {
        DEBUG("converting " << i)
        PathStorage<Graph>& storage = single_long_reads[i];
        vector<PathInfo<Graph> > paths = storage.GetAllPaths();
        PathStorageInfo<Graph> single_storage(paths,
                cfg::get().pe_params.long_reads.single_reads.filtering,
                cfg::get().pe_params.long_reads.single_reads.weight_priority,
                cfg::get().pe_params.long_reads.single_reads.unique_edge_priority);
        long_reads_libs.push_back(single_storage);
        DEBUG("done " << i)
    }
}

void PEResolving(conj_graph_pack& gp) {
    vector<size_t> indexes;
    vector<PathStorageInfo<Graph> > long_reads_libs;
    ConvertLongReads(gp.single_long_reads, long_reads_libs);
    std::string name = "scaffolds.fasta";
    bool traverse_loops = true;
    if (!(cfg::get().use_scaffolder && cfg::get().pe_params.param_set.scaffolder_options.on)) {
        name = "final_contigs.fasta";
        traverse_loops = false;
    }
    path_extend::ResolveRepeatsPe(
            gp, long_reads_libs, cfg::get().output_dir, name, traverse_loops,
            boost::optional<std::string>("final_contigs.fasta"));
}

void RepeatResolution::run(conj_graph_pack &gp, const char*) {
    OutputContigs(gp.g, cfg::get().additional_contigs, cfg::get().use_unipaths,
                  cfg::get().simp.tec.plausibility_length);
    OutputContigs(gp.g, cfg::get().output_dir + "before_rr.fasta");

    bool no_valid_libs = true;
    for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i)
        if (cfg::get().ds.reads[i].data().mean_insert_size != 0.0) {
            no_valid_libs = false;
            break;
        }

    if (cfg::get().rr_enable && no_valid_libs && !cfg::get().long_single_mode)
        WARN("Insert size was not estimated for any of the paired libraries, repeat resolution module will not run.");

    if ((no_valid_libs ||
            cfg::get().rm == debruijn_graph::resolving_mode::rm_none) &&
            !cfg::get().long_single_mode) {
        OutputContigs(gp.g, cfg::get().output_dir + "final_contigs.fasta");
        return;
    }

    // Repeat resolving begins
    if (cfg::get().rm == debruijn_graph::resolving_mode::rm_path_extend) {
        INFO("Path-Extend repeat resolving");
        PEResolving(gp);
    } else {
        INFO("Unsupported repeat resolver");
        OutputContigs(gp.g, cfg::get().output_dir + "final_contigs.fasta");
    }
}

void ContigOutput::run(conj_graph_pack &gp, const char*) {
    OutputContigs(gp.g, cfg::get().additional_contigs, cfg::get().use_unipaths,
                  cfg::get().simp.tec.plausibility_length);
    OutputContigs(gp.g, cfg::get().output_dir + "before_rr.fasta");
    OutputContigs(gp.g, cfg::get().output_dir + "final_contigs.fasta");
}


} // debruijn_graph

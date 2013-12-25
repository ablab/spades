//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
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

//TODO: get rid of this conversion
void ConvertLongReads(LongReadContainerT& single_long_reads, vector<PathStorageInfo<Graph> > &long_reads_libs) {
    for (size_t i = 0; i < single_long_reads.size(); ++i) {

        DEBUG("converting " << i)
        PathStorage<Graph>& storage = single_long_reads[i];
        vector<PathInfo<Graph> > paths = storage.GetAllPaths();
        auto single_read_param_set = cfg::get().pe_params.long_reads.single_reads;
        auto type = cfg::get().ds.reads[i].type();

        if (cfg::get().ds.reads[i].type() == io::LibraryType::PacBioReads  ||
                type == io::LibraryType::SangerReads) {
            single_read_param_set = cfg::get().pe_params.long_reads.pacbio_reads;
        } else if (type == io::LibraryType::TrustedContigs ||
                type == io::LibraryType::UntrustedContigs) {
            single_read_param_set = cfg::get().pe_params.long_reads.contigs;
        }

        auto tmp = single_read_param_set.unique_edge_priority;
        if (cfg::get().ds.single_cell &&  (cfg::get().ds.reads[i].type() == io::LibraryType::PacBioReads  ||
                type == io::LibraryType::SangerReads)) tmp = 10000.0;

        PathStorageInfo<Graph> single_storage(paths,
                single_read_param_set.filtering,
                single_read_param_set.weight_priority,
                tmp);
        long_reads_libs.push_back(single_storage);
        DEBUG("done " << i)
    }
}

void PEResolving(conj_graph_pack& gp) {
    vector<size_t> indexes;
    vector<PathStorageInfo<Graph> > long_reads_libs;
    ConvertLongReads(gp.single_long_reads, long_reads_libs);
    std::string name = "scaffolds";
    bool traverse_loops = true;
    if (!(cfg::get().use_scaffolder && cfg::get().pe_params.param_set.scaffolder_options.on)) {
        name = "final_contigs";
        traverse_loops = false;
    }
    path_extend::ResolveRepeatsPe(
            gp, long_reads_libs, cfg::get().output_dir, name, traverse_loops,
            boost::optional<std::string>("final_contigs"));
}

void RepeatResolution::run(conj_graph_pack &gp, const char*) {
    OutputContigs(gp.g, cfg::get().output_dir + "simplified_contigs", true, cfg::get().use_unipaths,
                  cfg::get().simp.tec.plausibility_length);
    OutputContigs(gp.g, cfg::get().output_dir + "before_rr", true);
    if (cfg::get().developer_mode) {
        FillPos(gp, gp.genome, "ref0");
        FillPos(gp, !gp.genome, "ref1");
//        gp.ClearQuality();
//        gp.FillQuality();
    }

    bool no_valid_libs = true;
    for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
        if (cfg::get().ds.reads[i].data().mean_insert_size != 0.0 ||
                cfg::get().ds.reads[i].is_pacbio_alignable()) {

            no_valid_libs = false;
            break;
        }
    }
    bool use_single_reads = cfg::get().always_single_reads_rr || (no_valid_libs && cfg::get().single_reads_rr);

    if (cfg::get().rr_enable && no_valid_libs && !use_single_reads)
        WARN("Insert size was not estimated for any of the paired libraries, repeat resolution module will not run.");

    if ((no_valid_libs ||
            cfg::get().rm == debruijn_graph::resolving_mode::rm_none) &&
            !use_single_reads) {
        OutputContigs(gp.g, cfg::get().output_dir + "final_contigs", true);
        return;
    }

    // Repeat resolving begins
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

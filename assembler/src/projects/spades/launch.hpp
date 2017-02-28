//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pipeline/config_struct.hpp"

#include "pipeline/graph_pack.hpp"
#include "stages/construction.hpp"
#include "pipeline/genomic_info_filler.hpp"
#include "gap_closer.hpp"
#include "stages/simplification.hpp"
#include "mismatch_correction.hpp"
#include "pair_info_count.hpp"
#include "second_phase_setup.hpp"
#include "repeat_resolving.hpp"
#include "distance_estimation.hpp"
#include "hybrid_aligning.hpp"
#include "chromosome_removal.hpp"
#include "series_analysis.hpp"
#include "pipeline/stage.hpp"
#include "contig_output_stage.hpp"

namespace spades {

inline bool MetaCompatibleLibraries() {
    const auto& libs = cfg::get().ds.reads;
    if (libs[0].type() != io::LibraryType::PairedEnd)
        return false;
    if (libs.lib_count() > 2)
        return false;
    if (libs.lib_count() == 2 &&
        libs[1].type() != io::LibraryType::TSLReads &&
        libs[1].type() != io::LibraryType::PacBioReads && libs[1].type() != io::LibraryType::NanoporeReads)
            return false;
    return true;
}

inline bool HybridLibrariesPresent() {
    for (size_t lib_id = 0; lib_id < cfg::get().ds.reads.lib_count(); ++lib_id) 
        if (cfg::get().ds.reads[lib_id].is_hybrid_lib()) 
            return true;
    return false;
}

void assemble_genome() {
    INFO("SPAdes started");
    if (cfg::get().mode == debruijn_graph::config::pipeline_type::meta && !MetaCompatibleLibraries()) {
        ERROR("Sorry, current version of metaSPAdes can work either with single library (paired-end only) "
                      "or in paired-end + TSLR mode.");
        exit(239);
    }

    INFO("Starting from stage: " << cfg::get().entry_point);

    bool two_step_rr = cfg::get().two_step_rr && cfg::get().rr_enable;
    INFO("Two-step RR enabled: " << two_step_rr);

    StageManager SPAdes({cfg::get().developer_mode,
                         cfg::get().load_from,
                         cfg::get().output_saves});

    size_t read_index_cnt = cfg::get().ds.reads.lib_count();
    if (two_step_rr)
        read_index_cnt++;

    debruijn_graph::conj_graph_pack conj_gp(cfg::get().K,
                                            cfg::get().tmp_dir,
                                            read_index_cnt,
                                            cfg::get().ds.reference_genome,
                                            cfg::get().flanking_range,
                                            cfg::get().pos.max_mapping_gap,
                                            cfg::get().pos.max_gap_diff);
    if (cfg::get().need_mapping) {
        INFO("Will need read mapping, kmer mapper will be attached");
        conj_gp.kmer_mapper.Attach();
    }
    // Build the pipeline
    SPAdes.add(new debruijn_graph::Construction())
          .add(new debruijn_graph::GenomicInfoFiller());
    if (cfg::get().gap_closer_enable && cfg::get().gc.before_simplify)
        SPAdes.add(new debruijn_graph::GapClosing("early_gapcloser"));

    SPAdes.add(new debruijn_graph::Simplification(two_step_rr));

    if (cfg::get().gap_closer_enable && cfg::get().gc.after_simplify)
        SPAdes.add(new debruijn_graph::GapClosing("late_gapcloser"));

    SPAdes.add(new debruijn_graph::SimplificationCleanup());
    //currently cannot be used with two step rr
    if (cfg::get().correct_mismatches && !cfg::get().two_step_rr)
        SPAdes.add(new debruijn_graph::MismatchCorrection());
    if (cfg::get().rr_enable) {
        if (two_step_rr) {
            string prelim_prefix = "preliminary_";
            if (cfg::get().use_intermediate_contigs)
                SPAdes.add(new debruijn_graph::PairInfoCount(true))
                      .add(new debruijn_graph::DistanceEstimation(true))
                      .add(new debruijn_graph::RepeatResolution(true))
                      .add(new debruijn_graph::ContigOutput(true, true, prelim_prefix))
                      .add(new debruijn_graph::SecondPhaseSetup(prelim_prefix));

            SPAdes.add(new debruijn_graph::Simplification());
        }

        if (!cfg::get().series_analysis.empty())
            SPAdes.add(new debruijn_graph::SeriesAnalysis());

        if (cfg::get().pd) {
            SPAdes.add(new debruijn_graph::ChromosomeRemoval());
        }

        if (HybridLibrariesPresent()) {
            SPAdes.add(new debruijn_graph::HybridLibrariesAligning());
        }

        //No graph modification allowed after HybridLibrariesAligning stage!

        SPAdes.add(new debruijn_graph::ContigOutput(false, false, "pre_pe_"))
              .add(new debruijn_graph::PairInfoCount())
              .add(new debruijn_graph::DistanceEstimation())
              .add(new debruijn_graph::RepeatResolution());
    }

    SPAdes.add(new debruijn_graph::ContigOutput());

    SPAdes.run(conj_gp, cfg::get().entry_point.c_str());

    // For informing spades.py about estimated params
    debruijn_graph::config::write_lib_data(path::append_path(cfg::get().output_dir, "final"));

    INFO("SPAdes finished");
}

}

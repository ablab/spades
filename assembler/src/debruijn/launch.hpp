//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "standard.hpp"

#include "config_struct.hpp"

#include "graph_pack.hpp"
#include "construction.hpp"
#include "genomic_info_filler.hpp"
#include "gap_closer.hpp"
#include "simplification.hpp"
#include "mismatch_correction.hpp"
#include "pair_info_count.hpp"
#include "repeat_resolving.hpp"
#include "distance_estimation.hpp"
#include "pacbio_aligning.hpp"
#include "stage.hpp"

namespace spades {

void assemble_genome() {
    INFO("SPAdes started");
    INFO("Starting from stage: " << cfg::get().entry_point);

    StageManager SPAdes({ cfg::get().developer_mode,
                          cfg::get().load_from,
                          cfg::get().output_saves});

    debruijn_graph::conj_graph_pack conj_gp(cfg::get().K,
                                            cfg::get().tmp_dir,
                                            cfg::get().ds.reads.lib_count(),
                                            cfg::get().ds.reference_genome,
                                            cfg::get().flanking_range,
                                            cfg::get().pos.max_mapping_gap,
                                            cfg::get().pos.max_gap_diff);
    if (!cfg::get().developer_mode) {
        conj_gp.edge_pos.Detach();
        conj_gp.paired_indices.Detach();
        conj_gp.clustered_indices.Detach();
        conj_gp.scaffolding_indices.Detach();
        conj_gp.element_finder.Detach();
        if (!cfg::get().gap_closer_enable && !cfg::get().rr_enable)
            conj_gp.kmer_mapper.Detach();
    }

    // Build the pipeline
    SPAdes.add(new debruijn_graph::Construction());
    SPAdes.add(new debruijn_graph::GenomicInfoFiller());
    if (cfg::get().gap_closer_enable && cfg::get().gc.before_simplify)
        SPAdes.add(new debruijn_graph::GapClosing("early_gapcloser"));
    SPAdes.add(new debruijn_graph::Simplification());
    if (cfg::get().gap_closer_enable && cfg::get().gc.after_simplify)
        SPAdes.add(new debruijn_graph::GapClosing("late_gapcloser"));
    SPAdes.add(new debruijn_graph::SimplificationCleanup());
    if (cfg::get().correct_mismatches)
        SPAdes.add(new debruijn_graph::MismatchCorrection());
    if (cfg::get().rr_enable) {
        bool run_pacbio = false;
        for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
            if (cfg::get().ds.reads[i].is_pacbio_alignable()) {
                run_pacbio = true;
                break;
            }
        }
        if (run_pacbio) {
            SPAdes.add(new debruijn_graph::PacBioAligning());
        }
        SPAdes.add(new debruijn_graph::PairInfoCount());
        SPAdes.add(new debruijn_graph::DistanceEstimation());
        SPAdes.add(new debruijn_graph::RepeatResolution());
    } else {
        SPAdes.add(new debruijn_graph::ContigOutput());
    }

    SPAdes.run(conj_gp, cfg::get().entry_point.c_str());

    // For informing spades.py about estimated params
    debruijn_graph::write_lib_data(path::append_path(cfg::get().output_dir, "final"));

    INFO("SPAdes finished");
}

}

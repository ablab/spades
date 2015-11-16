//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

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
#include "second_phase_setup.hpp"
#include "repeat_resolving.hpp"
#include "distance_estimation.hpp"
#include "pacbio_aligning.hpp"
#include "stage.hpp"

namespace spades {

void assemble_genome() {
    INFO("SPAdes started");
    if (cfg::get().ds.meta && cfg::get().ds.reads.lib_count() != 1) {
        ERROR("Sorry, current version of metaSPAdes can work with single library only (paired-end only).");
        exit(239);
    }

    INFO("Starting from stage: " << cfg::get().entry_point);

    bool two_step_rr = cfg::get().two_step_rr && cfg::get().rr_enable && cfg::get().ds.meta;
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
    SPAdes.add(new debruijn_graph::Construction());
    SPAdes.add(new debruijn_graph::GenomicInfoFiller());
    if (cfg::get().gap_closer_enable && cfg::get().gc.before_simplify)
        SPAdes.add(new debruijn_graph::GapClosing("early_gapcloser"));

    SPAdes.add(new debruijn_graph::Simplification(two_step_rr));

    if (cfg::get().gap_closer_enable && cfg::get().gc.after_simplify)
        SPAdes.add(new debruijn_graph::GapClosing("late_gapcloser"));

    SPAdes.add(new debruijn_graph::SimplificationCleanup());
    //currently cannot be used with two step rr
    if (cfg::get().correct_mismatches && !two_step_rr)
        SPAdes.add(new debruijn_graph::MismatchCorrection());
    if (cfg::get().rr_enable) {
        if (two_step_rr) {
            SPAdes.add(new debruijn_graph::PairInfoCount(true));
            SPAdes.add(new debruijn_graph::DistanceEstimation(true));
            SPAdes.add(new debruijn_graph::RepeatResolution(true));
            SPAdes.add(new debruijn_graph::SecondPhaseSetup());
            SPAdes.add(new debruijn_graph::Simplification());
        }

        //begin pacbio
        bool run_pacbio = false;
        for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
            if (cfg::get().ds.reads[i].is_pacbio_alignable()) {
                run_pacbio = true;
                break;
            }
        }
        if (run_pacbio) {
			//currently not integrated with two step rr process
			VERIFY(!two_step_rr);
            SPAdes.add(new debruijn_graph::PacBioAligning());
        }
        //end pacbio

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

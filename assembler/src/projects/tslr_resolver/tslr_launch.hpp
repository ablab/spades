#pragma once

#include <barcode_map_construction.hpp>
#include <barcode_repeat_resolution.hpp>
#include <projects/spades/pair_info_count.hpp>
#include <projects/spades/gap_closer.hpp>
#include <common/stages/construction.hpp>
#include <common/pipeline/genomic_info_filler.hpp>
#include <common/stages/simplification.hpp>
#include <projects/spades/hybrid_aligning.hpp>
#include <projects/spades/mismatch_correction.hpp>
#include "projects/spades/distance_estimation.hpp"

namespace spades {

    void run_tslr_resolver() {
        INFO("Starting from stage " << cfg::get().entry_point.c_str());

        debruijn_graph::conj_graph_pack conj_gp(cfg::get().K,
                                                cfg::get().tmp_dir,
                                                cfg::get().ds.reads.lib_count(),
                                                cfg::get().ds.reference_genome,
                                                cfg::get().flanking_range,
                                                cfg::get().pos.max_mapping_gap,
                                                cfg::get().pos.max_gap_diff);


        StageManager barcode_resolver_pipeline({cfg::get().developer_mode,
                              cfg::get().load_from,
                              cfg::get().output_saves});

        conj_gp.kmer_mapper.Attach();

        bool two_step_rr = cfg::get().rr_enable && cfg::get().two_step_rr;

        barcode_resolver_pipeline.add(new debruijn_graph::Construction())
                                 .add(new debruijn_graph::GenomicInfoFiller());

        if (cfg::get().gap_closer_enable && cfg::get().gc.before_simplify)
            barcode_resolver_pipeline.add(new debruijn_graph::GapClosing("early_gapcloser"));

        barcode_resolver_pipeline.add(new debruijn_graph::Simplification(true /*preliminary*/));

        if (cfg::get().gap_closer_enable && cfg::get().gc.after_simplify)
            barcode_resolver_pipeline.add(new debruijn_graph::GapClosing("late_gapcloser"));

        barcode_resolver_pipeline.add(new debruijn_graph::SimplificationCleanup());

        if (cfg::get().correct_mismatches && !cfg::get().two_step_rr)
            barcode_resolver_pipeline.add(new debruijn_graph::MismatchCorrection());
        if (cfg::get().rr_enable) {

            if (two_step_rr) {
//                if (cfg::get().use_intermediate_contigs)
//                    barcode_resolver_pipeline.add(new debruijn_graph::PairInfoCount(true))
//                                             .add(new debruijn_graph::DistanceEstimation(true))
//                                             .add(new debruijn_graph::RepeatResolution(true))
//                                             .add(new debruijn_graph::SecondPhaseSetup());
                barcode_resolver_pipeline.add(new debruijn_graph::Simplification());

            }
            barcode_resolver_pipeline.add(new debruijn_graph::HybridLibrariesAligning());

            //No graph modification allowed after HybridLibrariesAligning stage!

            barcode_resolver_pipeline.add(new BarcodeMapConstructionStage(cfg::get().K))
//                    .add(new debruijn_graph::PairInfoCount())
//                    .add(new debruijn_graph::DistanceEstimation())
                    .add(new TslrResolverStage(cfg::get().K, cfg::get().output_dir + "resolver_output.fasta"));


        }
//        else {
//            barcode_resolver_pipeline.add(new debruijn_graph::ContigOutput());
//        }

        INFO("Output directory: " << cfg::get().output_dir);

        barcode_resolver_pipeline.run(conj_gp, cfg::get().entry_point.c_str());

        INFO("Read cloud resolver finished.");
    }
} //spades

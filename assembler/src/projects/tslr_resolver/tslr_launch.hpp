#pragma once

#include <barcode_map_construction.hpp>
#include <barcode_repeat_resolution.hpp>
#include <projects/spades/pair_info_count.hpp>
#include <projects/spades/gap_closer.hpp>
#include <common/stages/construction.hpp>
#include <common/pipeline/genomic_info_filler.hpp>
#include <common/stages/simplification.hpp>
#include <projects/spades/hybrid_aligning.hpp>
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
        barcode_resolver_pipeline.add(new debruijn_graph::Construction())
                .add(new debruijn_graph::GenomicInfoFiller());

        if (!cfg::get().ts_res.ideal_reads) {
            barcode_resolver_pipeline.add(new debruijn_graph::GapClosing("early_gapcloser"))
                   .add(new debruijn_graph::Simplification)
                   .add(new debruijn_graph::GapClosing("late_gapcloser"))
                   .add(new debruijn_graph::SimplificationCleanup());
        }


        barcode_resolver_pipeline.add(new debruijn_graph::HybridLibrariesAligning());

        barcode_resolver_pipeline.add(new BarcodeMapConstructionStage(cfg::get().K))
                .add(new debruijn_graph::PairInfoCount())
                .add(new debruijn_graph::DistanceEstimation())
                .add(new TslrResolverStage(cfg::get().K, cfg::get().output_dir + "resolver_output.fasta"));
        INFO("Output directory: " << cfg::get().output_dir);

        barcode_resolver_pipeline.run(conj_gp, cfg::get().entry_point.c_str());
        INFO("TSLR resolver finished.");
    }
} //spades

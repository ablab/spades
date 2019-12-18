//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pipeline/config_struct.hpp"

#include "pipeline/graph_pack.hpp"
#include "stages/read_conversion.hpp"
#include "stages/construction.hpp"
#include "pipeline/genomic_info_filler.hpp"
#include "gap_closer.hpp"
#include "stages/simplification.hpp"
#include "stages/ss_edge_split.hpp"
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
#include "molecule_extraction_stage.hpp"

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

inline bool RNA10XLibraries() {
    for (size_t lib_id = 0; lib_id < cfg::get().ds.reads.lib_count(); ++lib_id)
        if (cfg::get().ds.reads[lib_id].type() == io::LibraryType::RNA10x)
            return true;
    return false;
}

inline std::string GetContigName(std::string contig_id, size_t cov) {
    std::string res = std::to_string(cov);
    while (res.length() < 4) {
        res = "_" + res;
    }
    return contig_id + res;
}

inline void AddMetaplasmidStages(StageManager &SPAdes) {
    size_t cov = cfg::get().pd->additive_step;
    size_t add = cfg::get().pd->additive_step;
    double multiplier = cfg::get().pd->relative_step;
    size_t max_cov = 600;
    SPAdes.add<debruijn_graph::ContigOutput>(true, true, GetContigName(cfg::get().co.contigs_name, 0));
    while (cov < max_cov) {
        SPAdes.add<debruijn_graph::ChromosomeRemoval>(cov);
        SPAdes.add<debruijn_graph::RepeatResolution>();
        SPAdes.add<debruijn_graph::ContigOutput>(true, true, GetContigName(cfg::get().co.contigs_name, cov));
        cov = std::max(cov + add, size_t((double) cov*multiplier));
    }

}

void assemble_genome() {
    INFO("SPAdes started");
    if (cfg::get().mode == debruijn_graph::config::pipeline_type::meta && !MetaCompatibleLibraries()) {
        ERROR("Sorry, current version of metaSPAdes can work either with single library (paired-end only) "
                      "or in paired-end + (TSLR or PacBio or Nanopore) mode.");
        exit(239);
    }

    INFO("Starting from stage: " << cfg::get().entry_point);

    StageManager SPAdes(SavesPolicy(cfg::get().checkpoints,
                                    cfg::get().output_saves, cfg::get().load_from));

    bool two_step_rr = cfg::get().two_step_rr && cfg::get().rr_enable;
    INFO("Two-step RR enabled: " << two_step_rr);

    debruijn_graph::conj_graph_pack conj_gp(cfg::get().K,
                                            cfg::get().tmp_dir,
                                            two_step_rr ? cfg::get().ds.reads.lib_count() + 1
                                                        : cfg::get().ds.reads.lib_count(),
                                            cfg::get().ds.reference_genome,
                                            cfg::get().flanking_range,
                                            cfg::get().pos.max_mapping_gap,
                                            cfg::get().pos.max_gap_diff);
    if (cfg::get().need_mapping) {
        INFO("Will need read mapping, kmer mapper will be attached");
        conj_gp.kmer_mapper.Attach();
    }

    // Build the pipeline
    SPAdes.add<ReadConversion>();
    SPAdes.add<debruijn_graph::Construction>();

    if (! debruijn_graph::config::PipelineHelper::IsMetagenomicPipeline(cfg::get().mode))
        SPAdes.add<debruijn_graph::GenomicInfoFiller>();

    VERIFY(!cfg::get().gc.before_raw_simplify || !cfg::get().gc.before_simplify);

    if (cfg::get().gap_closer_enable &&
        cfg::get().gc.before_raw_simplify)
        SPAdes.add<debruijn_graph::GapClosing>("early_gapcloser");

    //Using two_step_rr is hacky here. Fix soon!
    SPAdes.add<debruijn_graph::RawSimplification>(two_step_rr);

    if (cfg::get().gap_closer_enable &&
        cfg::get().gc.before_simplify)
        SPAdes.add<debruijn_graph::GapClosing>("early_gapcloser");

    if (two_step_rr) {
        SPAdes.add<debruijn_graph::Simplification>(true);
        if (cfg::get().gap_closer_enable && cfg::get().gc.after_simplify)
            SPAdes.add<debruijn_graph::GapClosing>("prelim_gapcloser");
        if (cfg::get().use_intermediate_contigs) {
            SPAdes.add<debruijn_graph::PairInfoCount>(true)
                    .add<debruijn_graph::DistanceEstimation>(true)
                    .add<debruijn_graph::RepeatResolution>(true)
                    .add<debruijn_graph::ContigOutput>(true)
                    .add<debruijn_graph::SecondPhaseSetup>();
        }
    }

    SPAdes.add<debruijn_graph::Simplification>();

    if (cfg::get().gap_closer_enable && cfg::get().gc.after_simplify)
        SPAdes.add<debruijn_graph::GapClosing>("late_gapcloser");

    SPAdes.add<debruijn_graph::SimplificationCleanup>();
    if (cfg::get().correct_mismatches)
        SPAdes.add<debruijn_graph::MismatchCorrection>();

    if (cfg::get().ss_coverage_splitter.enabled)
        SPAdes.add<debruijn_graph::SSEdgeSplit>();

    if (cfg::get().rr_enable) {
        if (!cfg::get().series_analysis.empty())
            SPAdes.add<debruijn_graph::SeriesAnalysis>();
//Not metaplasmid!
        if (cfg::get().mode == debruijn_graph::config::pipeline_type::plasmid)
            SPAdes.add<debruijn_graph::ChromosomeRemoval>();

        if (HybridLibrariesPresent()) {
            SPAdes.add<debruijn_graph::HybridLibrariesAligning>();
        }

        //No graph modification allowed after HybridLibrariesAligning stage!

        SPAdes.add<debruijn_graph::ContigOutput>(cfg::get().main_iteration, false)
               .add<debruijn_graph::PairInfoCount>();

        if (cfg::get().mode == debruijn_graph::config::pipeline_type::rna10x)
            SPAdes.add<debruijn_graph::MoleculeExtractionStage>();
        else {
            SPAdes.add<debruijn_graph::DistanceEstimation>()
                    .add<debruijn_graph::RepeatResolution>();
        }

        if (cfg::get().mode == debruijn_graph::config::pipeline_type::metaplasmid) {
            AddMetaplasmidStages(SPAdes);
        }

    } else {
        SPAdes.add<debruijn_graph::ContigOutput>(false);
    }

    SPAdes.add<debruijn_graph::ContigOutput>(cfg::get().main_iteration);

    SPAdes.run(conj_gp, cfg::get().entry_point.c_str());

    // For informing spades.py about estimated params
    debruijn_graph::config::write_lib_data(fs::append_path(cfg::get().output_dir, "final"));

    INFO("SPAdes finished");
}

}

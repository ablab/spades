//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "gap_closer.hpp"
#include "mismatch_correction.hpp"
#include "pair_info_count.hpp"
#include "pipeline/library.hpp"
#include "second_phase_setup.hpp"
#include "repeat_resolving.hpp"
#include "distance_estimation.hpp"
#include "hybrid_aligning.hpp"
#include "chromosome_removal.hpp"
#include "series_analysis.hpp"
#include "pipeline/stage.hpp"
#include "contig_output_stage.hpp"
#include "extract_domains.hpp"
#include "domain_graph_construction.hpp"
#include "restricted_edges_filling.hpp"

#include "modules/alignment/kmer_mapper.hpp"

#include "stages/genomic_info_filler.hpp"
#include "stages/read_conversion.hpp"
#include "stages/construction.hpp"
#include "stages/simplification.hpp"
#include "stages/ss_edge_split.hpp"

#include "pipeline/config_struct.hpp"
#include "pipeline/graph_pack.hpp"

#include "load_graph.hpp"

namespace spades {

static bool MetaCompatibleLibraries() {
    const auto& libs = cfg::get().ds.reads;
    if (libs.lib_count() > 2)
        return false;

    size_t paired_end_libs = 0, long_read_libs = 0;
    for (const auto &lib : libs) {
        auto type = lib.type();
        paired_end_libs += (type == io::LibraryType::PairedEnd);
        long_read_libs +=
            (type == io::LibraryType::TSLReads ||
             type == io::LibraryType::PacBioReads ||
             type == io::LibraryType::NanoporeReads);
    }

    return (paired_end_libs == 1 && long_read_libs <= 1);
}

static bool HybridLibrariesPresent() {
    for (const auto &lib : cfg::get().ds.reads)
        if (lib.is_hybrid_lib())
            return true;

    return false;
}

static bool AssemblyGraphPresent() {
    for (const auto &lib : cfg::get().ds.reads)
        if (lib.is_assembly_graph())
            return true;

    return false;
}

static std::string GetContigName(std::string contig_id, size_t cov) {
    std::string res = std::to_string(cov);
    while (res.length() < 4) {
        res = "_" + res;
    }
    return contig_id + res;
}

static debruijn_graph::ContigOutput::OutputList GetMetaplasmidOutput(size_t cov) {
    return {{debruijn_graph::ContigOutput::Kind::PlasmidContigs,
             GetContigName(cfg::get().co.contigs_name, cov) }};
}

static void AddMetaplasmidStages(StageManager &SPAdes) {
    size_t cov = cfg::get().pd->additive_step;
    size_t add = cfg::get().pd->additive_step;
    double multiplier = cfg::get().pd->relative_step;
    size_t max_cov = 600;
    SPAdes.add<debruijn_graph::ContigOutput>(GetMetaplasmidOutput(0));
    while (cov < max_cov) {
        SPAdes.add<debruijn_graph::ChromosomeRemoval>(cov);
        SPAdes.add<debruijn_graph::RepeatResolution>();
        SPAdes.add<debruijn_graph::ContigOutput>(GetMetaplasmidOutput(cov));
        cov = std::max(cov + add, size_t((double) cov*multiplier));
    }
}

static debruijn_graph::ContigOutput::OutputList GetPreliminaryStageOutput() {
    using namespace debruijn_graph;

    return {
        {ContigOutput::Kind::GFAGraph, "strain_graph"},
        {ContigOutput::Kind::FinalContigs, cfg::get().co.contigs_name}
    };
}

static debruijn_graph::ContigOutput::OutputList GetNonFinalStageOutput() {
    return { { debruijn_graph::ContigOutput::Kind::BinaryContigs, "simplified_contigs"} };
}

static debruijn_graph::ContigOutput::OutputList GetBeforeRROutput() {
    using namespace debruijn_graph;

    return {
        { ContigOutput::Kind::GFAGraph, "assembly_graph_after_simplification"},
        { ContigOutput::Kind::EdgeSequences, "before_rr"}
    };
}

static debruijn_graph::ContigOutput::OutputList GetFinalStageOutput() {
    using namespace debruijn_graph;

    return {
        { ContigOutput::Kind::EdgeSequences, "before_rr" },
        { ContigOutput::Kind::GFAGraph, "assembly_graph_with_scaffolds" },
        { ContigOutput::Kind::FASTGGraph, "assembly_graph" },
        { ContigOutput::Kind::FinalContigs, cfg::get().co.contigs_name },
        { ContigOutput::Kind::Scaffolds, cfg::get().co.scaffolds_name }
    };
}

static void AddPreliminarySimplificationStages(StageManager &SPAdes) {
    using namespace debruijn_graph::config;
    pipeline_type mode = cfg::get().mode;

    SPAdes.add<debruijn_graph::Simplification>(true);
    if (cfg::get().gap_closer_enable && cfg::get().gc.after_simplify)
        SPAdes.add<debruijn_graph::GapClosing>("prelim_gapcloser");

    if (cfg::get().use_intermediate_contigs) {
        SPAdes.add<debruijn_graph::PairInfoCount>(true)
              .add<debruijn_graph::DistanceEstimation>(true)
              .add<debruijn_graph::RepeatResolution>(true);

        if (cfg::get().hm)
            SPAdes.add<debruijn_graph::ExtractDomains>();

        SPAdes.add<debruijn_graph::ContigOutput>(GetPreliminaryStageOutput())
              .add<debruijn_graph::SecondPhaseSetup>();
        if (cfg::get().hm)
            SPAdes.add<debruijn_graph::RestrictedEdgesFilling>();
    }
}

static void AddSimplificationStages(StageManager &SPAdes) {
    VERIFY(!cfg::get().gc.before_raw_simplify || !cfg::get().gc.before_simplify);
    bool two_step_rr = cfg::get().two_step_rr && cfg::get().rr_enable;

    if (cfg::get().gap_closer_enable &&
        cfg::get().gc.before_raw_simplify)
        SPAdes.add<debruijn_graph::GapClosing>("early_gapcloser");

    // Using two_step_rr is hacky here. Fix soon!
    SPAdes.add<debruijn_graph::RawSimplification>(two_step_rr);

    if (cfg::get().gap_closer_enable &&
        cfg::get().gc.before_simplify)
        SPAdes.add<debruijn_graph::GapClosing>("early_gapcloser");

    if (two_step_rr)
        AddPreliminarySimplificationStages(SPAdes);

    SPAdes.add<debruijn_graph::Simplification>();

    if (cfg::get().gap_closer_enable && cfg::get().gc.after_simplify)
        SPAdes.add<debruijn_graph::GapClosing>("late_gapcloser");

    SPAdes.add<debruijn_graph::SimplificationCleanup>();

    if (cfg::get().correct_mismatches)
        SPAdes.add<debruijn_graph::MismatchCorrection>();

    if (cfg::get().ss_coverage_splitter.enabled)
        SPAdes.add<debruijn_graph::SSEdgeSplit>();
}

static void AddConstructionStages(StageManager &SPAdes) {
    using namespace debruijn_graph::config;
    pipeline_type mode = cfg::get().mode;

    SPAdes.add<debruijn_graph::Construction>();
    if (!PipelineHelper::IsMetagenomicPipeline(mode))
        SPAdes.add<debruijn_graph::GenomicInfoFiller>();
}

static void AddRepeatResolutionStages(StageManager &SPAdes) {
    using namespace debruijn_graph::config;

    if (!cfg::get().series_analysis.empty())
        SPAdes.add<debruijn_graph::SeriesAnalysis>();

    SPAdes.add<debruijn_graph::PairInfoCount>()
          .add<debruijn_graph::DistanceEstimation>()
          .add<debruijn_graph::RepeatResolution>();
}

void assemble_genome() {
    using namespace debruijn_graph::config;
    pipeline_type mode = cfg::get().mode;

    INFO("SPAdes started");

    // Perform various sanity checks
    if (mode == pipeline_type::meta && !MetaCompatibleLibraries()) {
        FATAL_ERROR("Sorry, current version of metaSPAdes can work either with single library (paired-end only) "
                    "or in hybrid paired-end + (TSLR or PacBio or Nanopore) mode.");
    } else if (AssemblyGraphPresent() &&
               (mode != pipeline_type::metaextrachromosomal &&
                !cfg::get().hm)) {
        // Disallow generic assembly graph inputs for now
        FATAL_ERROR("Assembly graph inputs are supported only for plasmid / metaextrachromosomal and / bgc modes!");
    }

    INFO("Starting from stage: " << cfg::get().entry_point);

    StageManager SPAdes(SavesPolicy(cfg::get().checkpoints,
                                    cfg::get().output_saves, cfg::get().load_from));

    bool two_step_rr = cfg::get().two_step_rr && cfg::get().rr_enable;
    INFO("Two-step repeat resolution " << (two_step_rr ? "enabled" : "disabled"));

    debruijn_graph::GraphPack conj_gp(cfg::get().K,
                                            cfg::get().tmp_dir,
                                            two_step_rr ? cfg::get().ds.reads.lib_count() + 1
                                                        : cfg::get().ds.reads.lib_count(),
                                            cfg::get().ds.reference_genome,
                                            cfg::get().flanking_range,
                                            cfg::get().pos.max_mapping_gap,
                                            cfg::get().pos.max_gap_diff);
    if (cfg::get().need_mapping) {
        INFO("Will need read mapping, kmer mapper will be attached");
        conj_gp.get_mutable<debruijn_graph::KmerMapper<debruijn_graph::Graph>>().Attach();
    }

    // Build the pipeline
    SPAdes.add<ReadConversion>();

    if (!AssemblyGraphPresent()) {
        AddConstructionStages(SPAdes);

        AddSimplificationStages(SPAdes);

        SPAdes.add<debruijn_graph::ContigOutput>(cfg::get().main_iteration ?
                                                 GetBeforeRROutput() : GetNonFinalStageOutput());
    } else {
        SPAdes.add<debruijn_graph::LoadGraph>();
    }
    
    if (cfg::get().main_iteration) {
        // Not metaextrachromosomal!
        if (mode == pipeline_type::plasmid)
            SPAdes.add<debruijn_graph::ChromosomeRemoval>();

        if (HybridLibrariesPresent())
            SPAdes.add<debruijn_graph::HybridLibrariesAligning>();

        // No graph modification allowed after HybridLibrariesAligning stage!

        if (cfg::get().rr_enable)
            AddRepeatResolutionStages(SPAdes);

        if (mode == pipeline_type::metaextrachromosomal)
            AddMetaplasmidStages(SPAdes);
        else
            SPAdes.add<debruijn_graph::ContigOutput>(GetFinalStageOutput());

        if (cfg::get().hm)
            SPAdes.add<debruijn_graph::DomainGraphConstruction>();
    }

    SPAdes.run(conj_gp, cfg::get().entry_point.c_str());

    // For informing spades.py about estimated params
    write_lib_data(fs::append_path(cfg::get().output_dir, "final"));

    INFO("SPAdes finished");
}

}

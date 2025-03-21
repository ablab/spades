//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "pair_info_count.hpp"

#include "alignment/bwa_sequence_mapper.hpp"
#include "alignment/long_read_mapper.hpp"
#include "alignment/rna/ss_coverage_filler.hpp"
#include "assembly_graph/core/basic_graph_stats.hpp"
#include "io/dataset_support/read_converter.hpp"
#include "paired_info/pair_info_filler.hpp"
#include "paired_info/paired_info_utils.hpp"
#include "pipeline/graph_pack_helpers.h"
#include "pipeline/sequence_mapper_gp_api.hpp"

namespace debruijn_graph {

namespace {

using paired_info::SequencingLib;

std::shared_ptr<SequenceMapper<Graph>> ChooseProperMapper(const graph_pack::GraphPack& gp,
                                                          const SequencingLib& library) {
    const auto &graph = gp.get<Graph>();

    if (library.type() == io::LibraryType::MatePairs) {
        INFO("Mapping mate-pairs using BWA-mem mapper");
        return std::make_shared<alignment::BWAReadMapper<Graph>>(graph);
    }

    if (library.data().unmerged_read_length < gp.k() && library.type() == io::LibraryType::PairedEnd) {
        INFO("Mapping PE reads shorter than K with BWA-mem mapper");
        return std::make_shared<alignment::BWAReadMapper<Graph>>(graph);
    }

    INFO("Selecting usual mapper");
    return MapperInstance(gp);
}

bool HasGoodRRLibs() {
    for (const auto &lib : cfg::get().ds.reads) {
        if (lib.is_contig_lib())
            continue;

        if (lib.is_paired() &&
            lib.data().mean_insert_size == 0.0)
            continue;

        if (lib.is_repeat_resolvable())
            return true;
    }

    return false;
}

bool HasOnlyMP() {
    for (const auto &lib : cfg::get().ds.reads) {
        if (lib.type() == io::LibraryType::PathExtendContigs)
            continue;

        if (lib.type() != io::LibraryType::MatePairs &&
            lib.type() != io::LibraryType::HQMatePairs)
            return false;
    }

    return true;
}

bool ShouldObtainLibCoverage() {
    return cfg::get().calculate_coverage_for_each_lib;
}

//todo improve logic
bool ShouldObtainSingleReadsPaths(size_t ilib) {
    using config::single_read_resolving_mode;
    using config::PipelineHelper;
    switch (cfg::get().single_reads_rr) {
        case single_read_resolving_mode::all:
            return true;
        case single_read_resolving_mode::only_single_libs:
            // Map when no PacBio/paried libs or only mate-pairs or single lib itself
            if (!HasGoodRRLibs() || HasOnlyMP() ||
                cfg::get().ds.reads[ilib].type() == io::LibraryType::SingleReads) {
                if (!PipelineHelper::IsMetagenomicPipeline(cfg::get().mode) ||
                    cfg::get().mode == config::pipeline_type::rnaviral) {
                    return true;
                } else {
                    WARN("Single reads are not used in metagenomic mode");
                }
            }
            break;
        case single_read_resolving_mode::none:
            break;
        default:
            CHECK_FATAL_ERROR(false, "Invalid mode value");
    }
    return false;
}

size_t ProcessSingleReads(graph_pack::GraphPack &gp, size_t ilib,
                          bool use_binary = true, bool map_paired = false) {
    //FIXME make const
    auto& reads = cfg::get_writable().ds.reads[ilib];
    const auto &graph = gp.get<Graph>();

    SequenceMapperNotifier notifier;

    auto &single_long_reads = gp.get_mutable<LongReadContainer<Graph>>()[ilib];
    auto& trusted_paths = gp.get_mutable<path_extend::TrustedPathsContainer>()[ilib];
    LongReadMapper read_mapper(graph, single_long_reads, trusted_paths, reads.type());

    if (ShouldObtainSingleReadsPaths(ilib) || reads.is_contig_lib()) {
        //FIXME pretty awful, would be much better if listeners were shared ptrs
        notifier.Subscribe(&read_mapper);
        cfg::get_writable().ds.reads[ilib].data().single_reads_mapped = true;
    }

    SSCoverageFiller ss_coverage_filler(graph, gp.get_mutable<SSCoverageContainer>()[ilib],
                                        !cfg::get().ss.ss_enabled);
    if (cfg::get().calculate_coverage_for_each_lib) {
        INFO("Will calculate lib coverage as well");
        map_paired = true;
        notifier.Subscribe(&ss_coverage_filler);
    }

    auto mapper_ptr = ChooseProperMapper(gp, reads);
    if (use_binary) {
        auto single_streams = single_binary_readers(reads, false, map_paired);
        notifier.ProcessLibrary(single_streams, *mapper_ptr);
    } else {
        auto single_streams = single_easy_readers(reads, false,
                                                  map_paired, /*handle Ns*/false);
        notifier.ProcessLibrary(single_streams, *mapper_ptr);
    }

    return single_long_reads.size();
}

} // namespace

void PairInfoCount::run(graph_pack::GraphPack &gp, const char *) {
    InitRRIndices(gp);
    EnsureBasicMapping(gp);

    const auto &graph = gp.get<Graph>();

    //TODO implement better universal logic
    size_t edge_length_threshold = cfg::get().min_edge_length_for_is_count;
    if (!debruijn_graph::config::PipelineHelper::IsMetagenomicPipeline(cfg::get().mode))
        edge_length_threshold = std::max(edge_length_threshold, Nx(graph, 50));

    INFO("Min edge length for estimation: " << edge_length_threshold);
    for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
        auto &lib = cfg::get_writable().ds.reads[i];
        if (lib.is_hybrid_lib()) {
            INFO("Library #" << i << " was mapped earlier on hybrid aligning stage, skipping");
            continue;
        } else if (lib.is_contig_lib()) {
            INFO("Mapping contigs library #" << i);
            ProcessSingleReads(gp, i, false);
        } else {
            if (lib.is_paired()) {
                INFO("Estimating insert size for library #" << i);
                const auto &lib_data = lib.data();
                size_t rl = lib_data.unmerged_read_length;
                size_t k = cfg::get().K;

                size_t edgepairs = 0;
                if (!paired_info::CollectLibInformation(graph, *ChooseProperMapper(gp, lib),
                                                        edgepairs, lib, edge_length_threshold)) {
                    cfg::get_writable().ds.reads[i].data().mean_insert_size = 0.0;
                    WARN("Unable to estimate insert size for paired library #" << i);
                    if (rl > 0 && rl <= k) {
                        WARN("Maximum read length (" << rl << ") should be greater than K (" << k << ")");
                    } else if (rl <= k * 11 / 10) {
                        WARN("Maximum read length (" << rl << ") is probably too close to K (" << k << ")");
                    } else {
                        WARN("None of paired reads aligned properly. Please, check orientation of your read pairs.");
                    }
                    continue;
                }

                INFO("  Insert size = " << lib_data.mean_insert_size <<
                     ", deviation = " << lib_data.insert_size_deviation <<
                     ", left quantile = " << lib_data.insert_size_left_quantile <<
                     ", right quantile = " << lib_data.insert_size_right_quantile <<
                     ", read length = " << lib_data.unmerged_read_length);

                if (lib_data.mean_insert_size < 1.1 * (double) rl)
                    WARN("Estimated mean insert size " << lib_data.mean_insert_size
                         << " is very small compared to read length " << rl);

                std::unique_ptr<paired_info::PairedInfoFilter> filter;
                unsigned filter_threshold = cfg::get().de.raw_filter_threshold;

                // Only filter paired-end libraries
                if (filter_threshold && lib.type() == io::LibraryType::PairedEnd) {
                    INFO("Filtering data for library #" << i);
                    filter = paired_info::FillEdgePairFilter(graph, *ChooseProperMapper(gp, lib), lib, edgepairs);
                }

                INFO("Mapping library #" << i);
                if (lib.data().mean_insert_size != 0.0) {
                    INFO("Mapping paired reads (takes a while) ");
                    using Indices = omnigraph::de::UnclusteredPairedInfoIndicesT<Graph>;

                    unsigned round_thr = 0;
                    // Do not round if filtering is disabled
                    if (filter)
                        round_thr = unsigned(std::min(cfg::get().de.max_distance_coeff * lib.data().insert_size_deviation * cfg::get().de.rounding_coeff,
                                                      cfg::get().de.rounding_thr));

                    paired_info::FillPairedIndex(graph, *ChooseProperMapper(gp, lib),
                                                 lib, gp.get_mutable<Indices>()[i],
                                                 std::move(filter), filter_threshold, round_thr);
                }
            }

            if (ShouldObtainSingleReadsPaths(i) || ShouldObtainLibCoverage()) {
                cfg::get_writable().use_single_reads |= ShouldObtainSingleReadsPaths(i);
                INFO("Mapping single reads of library #" << i);
                size_t n = ProcessSingleReads(gp, i, /*use_binary*/true, /*map_paired*/true);
                INFO("Total paths obtained from single reads: " << n);
            }
        }
    }
}

} // namespace debruijn_graph

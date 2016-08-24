//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include <paired_info/is_counter.hpp>
#include "io/dataset_support/read_converter.hpp"

#include "pair_info_count.hpp"
#include "assembly_graph/graph_alignment/short_read_mapper.hpp"
#include "assembly_graph/graph_alignment/long_read_mapper.hpp"
#include "paired_info/pair_info_filler.hpp"
#include "algorithms/path_extend/split_graph_pair_info.hpp"
#include "paired_info/bwa_pair_info_filler.hpp"

namespace debruijn_graph {


bool RefineInsertSizeForLib(conj_graph_pack &gp, size_t ilib, size_t edge_length_threshold) {

    INFO("Estimating insert size (takes a while)");
    InsertSizeCounter hist_counter(gp, edge_length_threshold);
    SequenceMapperNotifier notifier(gp);
    notifier.Subscribe(ilib, &hist_counter);

    auto& reads = cfg::get_writable().ds.reads[ilib];
    auto paired_streams = paired_binary_readers(reads, false);

    VERIFY(reads.data().read_length != 0);
    notifier.ProcessLibrary(paired_streams, ilib, *ChooseProperMapper(gp, reads));

    INFO(hist_counter.mapped() << " paired reads (" <<
         ((double) hist_counter.mapped() * 100.0 / (double) hist_counter.total()) <<
         "% of all) aligned to long edges");
    if (hist_counter.negative() > 3 * hist_counter.mapped())
        WARN("Too much reads aligned with negative insert size. Is the library orientation set properly?");
    if (hist_counter.mapped() == 0)
        return false;

    std::map<size_t, size_t> percentiles;
    hist_counter.FindMean(reads.data().mean_insert_size, reads.data().insert_size_deviation, percentiles);
    hist_counter.FindMedian(reads.data().median_insert_size, reads.data().insert_size_mad,
                            reads.data().insert_size_distribution);
    if (reads.data().median_insert_size < gp.k_value + 2) {
        return false;
    }

    std::tie(reads.data().insert_size_left_quantile,
             reads.data().insert_size_right_quantile) = omnigraph::GetISInterval(0.8,
                                                                                 reads.data().insert_size_distribution);

    return !reads.data().insert_size_distribution.empty();
}

void ProcessSingleReads(conj_graph_pack &gp,
                        size_t ilib,
                        bool use_binary,
                        bool map_paired,
                        SequenceMapperListener* mapping_listener_ptr) {
    //FIXME make const
    auto& reads = cfg::get_writable().ds.reads[ilib];

    SequenceMapperNotifier notifier(gp);
    //FIXME pretty awful, would be much better if listeners were shared ptrs
    LongReadMapper read_mapper(gp.g, gp.single_long_reads[ilib],
                               ChooseProperReadPathExtractor(gp.g, reads.type()));

    if (mapping_listener_ptr != nullptr) {
        notifier.Subscribe(ilib, mapping_listener_ptr);
    }

    notifier.Subscribe(ilib, &read_mapper);

    auto mapper_ptr = ChooseProperMapper(gp, reads);
    if (use_binary) {
        auto single_streams = single_binary_readers(reads, false, map_paired);
        notifier.ProcessLibrary(single_streams, ilib, *mapper_ptr);
    } else {
        auto single_streams = single_easy_readers(reads, false,
                                                  map_paired, /*handle Ns*/false);
        notifier.ProcessLibrary(single_streams, ilib, *mapper_ptr);
    }
    cfg::get_writable().ds.reads[ilib].data().single_reads_mapped = true;
}

void ProcessPairedReads(conj_graph_pack &gp, size_t ilib) {
    auto& reads = cfg::get_writable().ds.reads[ilib];
    bool calculate_threshold = (reads.type() == io::LibraryType::PairedEnd);
    SequenceMapperNotifier notifier(gp);
    INFO("Left insert size qauntile " << reads.data().insert_size_left_quantile <<
         ", right insert size quantile " << reads.data().insert_size_right_quantile);

    path_extend::SplitGraphPairInfo split_graph(
            gp, (size_t) reads.data().median_insert_size,
            (size_t) reads.data().insert_size_deviation,
            (size_t) reads.data().insert_size_left_quantile,
            (size_t) reads.data().insert_size_right_quantile,
            reads.data().read_length, gp.g.k(),
            cfg::get().pe_params.param_set.split_edge_length,
            reads.data().insert_size_distribution);
    if (calculate_threshold) {
        notifier.Subscribe(ilib, &split_graph);
    }

    LatePairedIndexFiller pif(gp.g, PairedReadCountWeight, gp.paired_indices[ilib]);
    notifier.Subscribe(ilib, &pif);

    auto paired_streams = paired_binary_readers(reads, false, (size_t) reads.data().mean_insert_size);
    notifier.ProcessLibrary(paired_streams, ilib, *ChooseProperMapper(gp, reads));
    cfg::get_writable().ds.reads[ilib].data().pi_threshold = split_graph.GetThreshold();
}

bool HasGoodRRLibs() {
    for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
        const auto &lib = cfg::get().ds.reads[i];
        if (lib.is_contig_lib())
            continue;
        if (lib.is_paired() &&
            lib.data().mean_insert_size == 0.0) {
            continue;
        }
        if (lib.is_repeat_resolvable()) {
            return true;
        }
    }
    return false;
}

bool HasOnlyMP() {
    for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
        if (cfg::get().ds.reads[i].type() == io::LibraryType::PathExtendContigs)
            continue;
        if (cfg::get().ds.reads[i].type() != io::LibraryType::MatePairs &&
            cfg::get().ds.reads[i].type() != io::LibraryType::HQMatePairs) {
            return false;
        }
    }
    return true;
}

//todo improve logic
bool ShouldMapSingleReads(size_t ilib) {
    using config::single_read_resolving_mode;
    switch (cfg::get().single_reads_rr) {
        case single_read_resolving_mode::all: {
            return true;
        }
        case single_read_resolving_mode::only_single_libs: {
            //Map when no PacBio/paried libs or only mate-pairs or single lib itself
            if (!HasGoodRRLibs() || HasOnlyMP() ||
                cfg::get().ds.reads[ilib].type() == io::LibraryType::SingleReads) {
                if (cfg::get().mode != debruijn_graph::config::pipeline_type::meta) {
                    return true;
                } else {
                    WARN("Single reads are not used in metagenomic mode");
                }
            }
            break;
        }
        case single_read_resolving_mode::none: {
            break;
        }
        default:
            VERIFY_MSG(false, "Invalid mode value");
    }
    return false;
}

void PairInfoCount::run(conj_graph_pack &gp, const char *) {
    gp.InitRRIndices();
    gp.EnsureBasicMapping();

    //fixme implement better universal logic
    size_t edge_length_threshold = cfg::get().mode == config::pipeline_type::meta ? 1000 : stats::Nx(gp.g, 50);
    INFO("Min edge length for estimation: " << edge_length_threshold);
    bwa_pair_info::BWAPairInfoFiller bwa_counter(gp.g,
                                                 cfg::get().bwa.path_to_bwa,
                                                 path::append_path(cfg::get().output_dir, "bwa_count"),
                                                 cfg::get().max_threads, !cfg::get().bwa.debug);

    for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
        const auto &lib = cfg::get().ds.reads[i];

        if (cfg::get().bwa.bwa_enable && lib.is_bwa_alignable()) {
            //Run insert size estimation and pair index filler together to save disc space (removes SAM file right after processing the lib)
            bwa_counter.ProcessLib(i, cfg::get_writable().ds.reads[i], gp.paired_indices[i],
                                   edge_length_threshold, cfg::get().bwa.min_contig_len);
        } else if (lib.is_paired()) {
            INFO("Estimating insert size for library #" << i);
            const auto &lib_data = lib.data();
            size_t rl = lib_data.read_length;
            size_t k = cfg::get().K;
            bool insert_size_refined = RefineInsertSizeForLib(gp, i, edge_length_threshold);

            if (!insert_size_refined) {
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
            } else {
                INFO("  Insert size = " << lib_data.mean_insert_size <<
                     ", deviation = " << lib_data.insert_size_deviation <<
                     ", left quantile = " << lib_data.insert_size_left_quantile <<
                     ", right quantile = " << lib_data.insert_size_right_quantile <<
                     ", read length = " << lib_data.read_length);

                if (lib_data.mean_insert_size < 1.1 * (double) rl) {
                    WARN("Estimated mean insert size " << lib_data.mean_insert_size
                         << " is very small compared to read length " << rl);
                }
            }
        }
    }

    for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
        const auto &lib = cfg::get().ds.reads[i];
        if (lib.is_hybrid_lib()) {
            INFO("Library #" << i << " was mapped earlier on hybrid aligning stage, skipping");
            continue;
        } else if (lib.is_contig_lib()) {
            INFO("Mapping contigs library #" << i);
            ProcessSingleReads(gp, i, /*use_binary*/false);
        } else if (cfg::get().bwa.bwa_enable && lib.is_bwa_alignable()) {
            INFO("Library #" << i << " was mapped by BWA, skipping");
            continue;
        } else {
            INFO("Mapping library #" << i);
            bool map_single_reads = ShouldMapSingleReads(i);
            cfg::get_writable().use_single_reads |= map_single_reads;

            if (lib.is_paired() && lib.data().mean_insert_size != 0.0) {
                INFO("Mapping paired reads (takes a while) ");
                ProcessPairedReads(gp, i);
            }
            if (map_single_reads) {
                INFO("Mapping single reads (takes a while) ");
                ProcessSingleReads(gp, i, /*use_binary*/true, /*map_paired*/true);
                INFO("Total paths obtained from single reads: " << gp.single_long_reads[i].size());
            }

        }
    }

    SensitiveReadMapper<Graph>::EraseIndices();
}

}

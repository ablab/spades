//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "standard.hpp"
#include "dataset_readers.hpp"
#include "read_converter.hpp"

#include "de/paired_info.hpp"

#include "utils.hpp"
#include "stats/debruijn_stats.hpp"

#include "is_counter.hpp"
#include "pair_info_count.hpp"
#include "sequence_mapper.hpp"
#include "short_read_mapper.hpp"
#include "long_read_mapper.hpp"
#include "pair_info_filler.hpp"
#include "stats/debruijn_stats.hpp"
#include "path_extend/split_graph_pair_info.hpp"

namespace debruijn_graph {
    typedef io::SequencingLibrary<debruijn_config::DataSetData> SequencingLib;


bool RefineInsertSizeForLib(conj_graph_pack& gp, size_t ilib, size_t edge_length_threshold) {

  INFO("Estimating insert size (takes a while)");
  InsertSizeCounter hist_counter(gp, edge_length_threshold, /* ignore negative */ true);
  SequenceMapperNotifier notifier(gp, false);
  notifier.Subscribe(ilib, &hist_counter);

  VERIFY(cfg::get().use_multithreading);
  SequencingLib& reads = cfg::get_writable().ds.reads[ilib];
  VERIFY(reads.data().read_length != 0);
  auto paired_streams = paired_binary_readers(reads, false, (size_t) reads.data().mean_insert_size);
  notifier.ProcessLibrary(paired_streams, ilib, *ChooseProperMapper(gp, reads), paired_streams.size());

  INFO(hist_counter.mapped() << " paired reads (" << ((double) hist_counter.mapped() * 100.0 / (double) hist_counter.total()) << "% of all) aligned to long edges");
  if (hist_counter.negative() > 3 * hist_counter.mapped())
      WARN("Too much reads aligned with negative insert size. Does the library orientation set properly?");
  if (hist_counter.mapped() == 0)
    return false;

  std::map<size_t, size_t> percentiles;
  hist_counter.FindMean(reads.data().mean_insert_size, reads.data().insert_size_deviation, percentiles);
  hist_counter.FindMedian(reads.data().median_insert_size, reads.data().insert_size_mad, reads.data().insert_size_distribution);

  std::tie(reads.data().insert_size_left_quantile, reads.data().insert_size_right_quantile) = omnigraph::GetISInterval(0.8, reads.data().insert_size_distribution);

  return !reads.data().insert_size_distribution.empty();
}

void ProcessSingleReads(conj_graph_pack& gp, size_t ilib) {
    const SequencingLib& reads = cfg::get().ds.reads[ilib];
    SequenceMapperNotifier notifier(gp);
    SimpleLongReadMapper read_mapper(gp, gp.single_long_reads[ilib]);
    notifier.Subscribe(ilib, &read_mapper);

    VERIFY(cfg::get().use_multithreading);
    auto single_streams = single_binary_readers(reads, true, false);
    notifier.ProcessLibrary(single_streams, ilib, *ChooseProperMapper(gp, reads), single_streams.size());
}

void ProcessPairedReads(conj_graph_pack& gp, size_t ilib, bool map_single_reads) {
    const SequencingLib& reads = cfg::get().ds.reads[ilib];
    bool calculate_threshold = (reads.type() == io::LibraryType::PairedEnd);
    SequenceMapperNotifier notifier(gp);
    INFO("Left insert size qauntile " << reads.data().insert_size_left_quantile << ", right insert size quantile " << reads.data().insert_size_right_quantile);

    SimpleLongReadMapper read_mapper(gp, gp.single_long_reads[ilib]);
    if (map_single_reads) {
        notifier.Subscribe(ilib, &read_mapper);
    }

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

    VERIFY(cfg::get().use_multithreading);
    auto paired_streams = paired_binary_readers(reads, true, (size_t) reads.data().mean_insert_size);
    notifier.ProcessLibrary(paired_streams, ilib, *ChooseProperMapper(gp, reads), paired_streams.size());
    cfg::get_writable().ds.reads[ilib].data().pi_threshold = split_graph.GetThreshold();

    if (map_single_reads) {
        ProcessSingleReads(gp, ilib);
    }
}

static bool HasGoodRRLibs() {
    static bool has_good_rr_reads = false;
    for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
        if (cfg::get().ds.reads[i].is_paired() &&
                (cfg::get().ds.reads[i].data().mean_insert_size == 0.0 ||
                cfg::get().ds.reads[i].data().mean_insert_size < 1.1 * (double) cfg::get().ds.reads[i].data().read_length)) {
            continue;
        }
        if (cfg::get().ds.reads[i].is_repeat_resolvable()) {
            has_good_rr_reads = true;
            break;
        }
    }
    return has_good_rr_reads;
}

static bool HasOnlyMP() {
    for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
        if (cfg::get().ds.reads[i].type() != io::LibraryType::MatePairs && cfg::get().ds.reads[i].type() != io::LibraryType::HQMatePairs) {
            return false;
        }
    }
    return true;
}

bool ShouldMapSingleReads(size_t ilib) {
    switch (cfg::get().single_reads_rr) {
        case sr_none: {
            return false;
        }
        case sr_all: {
            return true;
        }
        case sr_only_single_libs: {
            //Map when no PacBio/paried libs or only mate-pairs or single lib itself
            return !HasGoodRRLibs() || HasOnlyMP() || (cfg::get().ds.reads[ilib].type() == io::LibraryType::SingleReads);
        }
    }
    return false;
}

void PairInfoCount::run(conj_graph_pack &gp, const char*) {
    if (!cfg::get().developer_mode) {
        gp.paired_indices.Attach();
        gp.paired_indices.Init();
    }

    size_t edge_length_threshold = stats::Nx(gp.g, 50);
    INFO("Graph N50: " << edge_length_threshold);
    for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
        INFO("Estimating insert size for library #" << i);
        if (cfg::get().ds.reads[i].is_paired()) {

            bool insert_size_refined = RefineInsertSizeForLib(gp, i, edge_length_threshold);

            if (!insert_size_refined) {
                cfg::get_writable().ds.reads[i].data().mean_insert_size = 0.0;
                WARN("Unable to estimate insert size for paired library #" << i);
                if (cfg::get().ds.reads[i].data().read_length > 0 && cfg::get().ds.reads[i].data().read_length <= cfg::get().K) {
                    WARN("Maximum read length (" << cfg::get().ds.reads[i].data().read_length << ") should be greater than K (" << cfg::get().K << ")");
                }
                else if (cfg::get().ds.reads[i].data().read_length <= cfg::get().K * 11 / 10) {
                    WARN("Maximum read length (" << cfg::get().ds.reads[i].data().read_length << ") is probably too close to K (" << cfg::get().K << ")");
                } else {
                    WARN("None of paired reads aligned properly. Please, check orientation of your read pairs.");
                }
                continue;
            } else {
                INFO("  Estimated insert size for paired library #" << i);
                INFO("  Insert size = " << cfg::get().ds.reads[i].data().mean_insert_size <<
                        ", deviation = " << cfg::get().ds.reads[i].data().insert_size_deviation <<
                        ", left quantile = " << cfg::get().ds.reads[i].data().insert_size_left_quantile <<
                        ", right quantile = " << cfg::get().ds.reads[i].data().insert_size_right_quantile <<
                        ", read length = " << cfg::get().ds.reads[i].data().read_length);
            }
        }
    }



    for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
        bool map_single_reads = ShouldMapSingleReads(i);
        cfg::get_writable().use_single_reads |= map_single_reads;


        INFO("Mapping library #" << i);
        if (cfg::get().ds.reads[i].is_paired() && cfg::get().ds.reads[i].data().mean_insert_size != 0.0) {
            INFO("Mapping paired reads (takes a while) ");
            ProcessPairedReads(gp, i, map_single_reads);
        }
        else if (map_single_reads) {
            INFO("Mapping single reads (takes a while) ");
            ProcessSingleReads(gp, i);
        }

        if (map_single_reads) {
            INFO("Total paths obtained from single reads: " << gp.single_long_reads[i].size());
        }
    }

    SensitiveReadMapper<Graph>::EraseIndices();
}

}

//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "standard.hpp"
#include "dataset_readers.hpp"
#include "read_converter.hpp"

#include "de/insert_size_refiner.hpp"
#include "de/paired_info.hpp"

#include "utils.hpp"

#include "pair_info_count.hpp"
#include "long_read_mapper.hpp"
#include "pair_info_filler.hpp"
#include "path_extend/split_graph_pair_info.hpp"

namespace debruijn_graph {

typedef io::ReadStreamVector<io::SequencePairedReadStream> MultiStreamType;
typedef io::ReadStreamVector<io::PairedReadStream> SingleStreamType;


void ProcessSingleReads(conj_graph_pack& gp, size_t ilib) {
    const io::SequencingLibrary<debruijn_config::DataSetData>& reads = cfg::get().ds.reads[ilib];
    SequenceMapperNotifier notifier(gp);
    SimpleLongReadMapper read_mapper(gp, gp.single_long_reads[ilib]);
    notifier.Subscribe(ilib, &read_mapper);

    if (cfg::get().use_multithreading) {
        auto single_streams = single_binary_readers(reads, true, false);
        notifier.ProcessLibrary(*single_streams, ilib, single_streams->size());
    }
    else {
        auto single_stream = single_easy_reader(reads, true, false);
        single_stream.release();
        SingleStreamType single_streams(single_stream.get());
        notifier.ProcessLibrary(single_streams, ilib, single_streams.size());
    }
}


void ProcessPairedReads(conj_graph_pack& gp, size_t ilib) {
    SequenceMapperNotifier notifier(gp);
    const io::SequencingLibrary<debruijn_config::DataSetData>& reads = cfg::get().ds.reads[ilib];

    path_extend::SplitGraphPairInfo split_graph(
            gp, (size_t) reads.data().mean_insert_size, reads.data().read_length,
            (size_t) reads.data().insert_size_deviation, gp.g.k(),
            cfg::get().pe_params.param_set.split_edge_length);

    LatePairedIndexFiller pif(gp.g, PairedReadCountWeight, gp.paired_indices[ilib]);

    SimpleLongReadMapper read_mapper(gp, gp.single_long_reads[ilib]);
    if (cfg::get().long_single_mode) {
        notifier.Subscribe(ilib, &read_mapper);
    }
    notifier.Subscribe(ilib, &split_graph);
    notifier.Subscribe(ilib, &pif);

    if (cfg::get().use_multithreading) {
        auto paired_streams = paired_binary_readers(reads, true, (size_t) reads.data().mean_insert_size);
        notifier.ProcessLibrary(*paired_streams, ilib, paired_streams->size());
        cfg::get_writable().ds.reads[ilib].data().pi_threshold = split_graph.GetThreshold();
    }
    else {
        auto paired_stream = paired_easy_reader(reads, true, (size_t) reads.data().mean_insert_size);
        SingleStreamType paired_streams(paired_stream.get());
        notifier.ProcessLibrary(paired_streams, ilib, paired_streams.size());
        cfg::get_writable().ds.reads[ilib].data().pi_threshold = split_graph.GetThreshold();
    }

    if (cfg::get().long_single_mode) {
        ProcessSingleReads(gp, ilib);
    }
}

template<class graph_pack, class PairedRead, class ConfigType>
bool RefineInsertSize(const graph_pack& gp,
                      io::ReadStreamVector<io::IReader<PairedRead> >& streams,
                      ConfigType& config,
                      size_t edge_length_threshold) {
  size_t rl;
  double mean;
  double delta;
  double median;
  double mad;
  std::map<size_t, size_t> percentiles;
  std::map<int, size_t> hist;
  // calling default method
  refine_insert_size(streams, gp, edge_length_threshold, rl, mean, delta, median, mad, percentiles, hist);

  if (hist.size() == 0) {
    config.paired_mode = false;
    WARN("Failed to estimate the insert size of paired reads, because none of the paired reads aligned to long edges.");
    WARN("Paired reads will not be used.");
    return false;
  }

  config.ds.set_IS(mean);
  config.ds.set_is_var(delta);
  config.ds.set_median(median);
  config.ds.set_mad(mad);
  config.ds.set_hist(hist);
  INFO("Mean Insert Size = " << mean);
  INFO("Insert Size stddev= " << delta);
  INFO("Median Insert Size = " << median);
  INFO("Insert Size MAD = " << mad);
  DEBUG("Delta_Mad = " << 1.4826 * mad);

  return true;
}

template<class graph_pack, class PairedRead, class DataSet>
bool RefineInsertSizeForLib(const graph_pack& gp,
                      io::ReadStreamVector<io::IReader<PairedRead> >& streams,
                      DataSet& data,
                      size_t edge_length_threshold) {

  std::map<size_t, size_t> percentiles;
  // calling default method
  data.read_length = 0;
  refine_insert_size(streams, gp, edge_length_threshold,
          data.read_length,
          data.mean_insert_size,
          data.insert_size_deviation,
          data.median_insert_size,
          data.insert_size_mad,
          percentiles,
          data.insert_size_distribution);

  if (data.insert_size_distribution.size() == 0) {
    return false;
  }

  return true;
}

void PairInfoCount::run(conj_graph_pack &gp) {
    if (!cfg::get().developer_mode) {
        gp.paired_indices.Attach();
        gp.paired_indices.Init();
    }

    size_t edge_length_threshold = Nx(gp.g, 50);
    for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {

        if (cfg::get().ds.reads[i].data().read_length > 0 && cfg::get().ds.reads[i].data().read_length <= cfg::get().K) {
            WARN("Unable to estimate insert size for paired library #" << i);
            WARN("Maximum read length (" << cfg::get().ds.reads[i].data().read_length << ") should be greater than K (" << cfg::get().K << ")");
            //TODO: run short read aligner in this case
            continue;
        }

        if (cfg::get().ds.reads[i].type() == io::LibraryType::PairedEnd ||
            cfg::get().ds.reads[i].type() == io::LibraryType::MatePairs) {

            bool insert_size_refined;
            if (cfg::get().use_multithreading) {
                auto streams = paired_binary_readers(cfg::get().ds.reads[i], false, 0);
                insert_size_refined = RefineInsertSizeForLib(gp, *streams, cfg::get_writable().ds.reads[i].data(), edge_length_threshold);
            } else {
                std::auto_ptr<PairedReadStream> stream = paired_easy_reader(cfg::get().ds.reads[i], false, 0);
                SingleStreamType streams(stream.get());
                streams.release();
                insert_size_refined = RefineInsertSizeForLib(gp, streams, cfg::get_writable().ds.reads[i].data(), edge_length_threshold);
            }

            if (!insert_size_refined) {
                cfg::get_writable().ds.reads[i].data().mean_insert_size = 0.0;
                WARN("Unable to estimate insert size for paired library #" << i);
                if (cfg::get().ds.reads[i].data().read_length <= cfg::get().K * 11 / 10) {
                    WARN("Maximum read length (" << cfg::get().ds.reads[i].data().read_length << ") is probably too close to K (" << cfg::get().K << ")");
                } else {
                    WARN("None of paired reads aligned properly. Please, check orientation of your read pairs.");
                }
                continue;
            } else {
                INFO("  Estimated insert size for paired library #" << i);
                INFO("  Insert size = " << cfg::get().ds.reads[i].data().mean_insert_size <<
                        ", deviation = " << cfg::get().ds.reads[i].data().insert_size_deviation <<
                        ", read length = " << cfg::get().ds.reads[i].data().read_length);
            }
            ProcessPairedReads(gp, i);
        }

        if (cfg::get().ds.reads[i].type() == io::LibraryType::SingleReads) {
            ProcessSingleReads(gp, i);
        }
    }
}

}

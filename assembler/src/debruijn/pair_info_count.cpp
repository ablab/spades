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

template<class graph_pack, class PairedReadType>
void RefineInsertSizeParallel(const graph_pack& gp,
        io::ReadStreamList<PairedReadType>& streams,
        InsertSizeHistogramCounter<graph_pack>& counter,
        size_t& rl) {

  debruijn_graph::MapperFactory<graph_pack> mapper_factory(gp);
  auto mapper = mapper_factory.GetSequenceMapper(rl == 0 ? gp.k_value + 1 : rl);

  size_t nthreads = streams.size();
  std::vector<size_t> rls(nthreads, 0);
  counter.Init(nthreads);

#pragma omp parallel for num_threads(nthreads)
  for (size_t i = 0; i < nthreads; ++i) {
    PairedReadType r;
    auto& stream = streams[i];
    stream.reset();

    while (!stream.eof()) {
      stream >> r;
      Sequence sequence_left = r.first().sequence();
      Sequence sequence_right = r.second().sequence();

      if (sequence_left.size() > rls[i]) {
          rls[i] = sequence_left.size();
      }
      if (sequence_right.size() > rls[i]) {
          rls[i] = sequence_right.size();
      }
      auto pos_left = mapper->GetLastKmerPos(sequence_left);
      auto pos_right = mapper->GetFirstKmerPos(sequence_right);
      counter.ProcessPairedRead(i, r, pos_left, pos_right);
    }
  }
  counter.Finalize();

  if (rl == 0) {
    rl = rls[0];
    for (size_t i = 1; i < nthreads; ++i) {
      if (rl < rls[i]) {
        rl = rls[i];
      }
    }
  }
}


template<class graph_pack, class PairedReadType, class DataSet>
bool RefineInsertSizeForLib(const graph_pack& gp,
                      io::ReadStreamList<PairedReadType>& streams,
                      DataSet& data,
                      size_t edge_length_threshold) {

  INFO("Estimating insert size (takes a while)");
  InsertSizeHistogramCounter<graph_pack> hist_counter(gp, edge_length_threshold, /* ignore negative */ true);
  RefineInsertSizeParallel<graph_pack, PairedReadType>(gp, streams, hist_counter, data.read_length);

  INFO(hist_counter.mapped() << " paired reads (" << ((double) hist_counter.mapped() * 100.0 / (double) hist_counter.total()) << "% of all) aligned to long edges");
  if (hist_counter.negative() > 3 * hist_counter.mapped())
      WARN("Too much reads aligned with negative insert size. Does the library orientation set properly?");
  if (hist_counter.mapped() == 0)
    return false;

  std::map<size_t, size_t> percentiles;
  hist_counter.FindMean(data.mean_insert_size,  data.insert_size_deviation, percentiles);
  hist_counter.FindMedian(data.median_insert_size, data.insert_size_mad, data.insert_size_distribution);

  omnigraph::ISInterval(0.8, data.mean_insert_size,
          data.insert_size_distribution,
          data.insert_size_left_quantile,
          data.insert_size_right_quantile);

  if (data.insert_size_distribution.size() == 0) {
    return false;
  }
  return true;
}

void ProcessSingleReads(conj_graph_pack& gp, size_t ilib) {
    const io::SequencingLibrary<debruijn_config::DataSetData>& reads = cfg::get().ds.reads[ilib];
    SequenceMapperNotifier notifier(gp);
    SimpleLongReadMapper read_mapper(gp, gp.single_long_reads[ilib]);
    notifier.Subscribe(ilib, &read_mapper);

    if (cfg::get().use_multithreading) {
        auto single_streams = single_binary_readers(reads, true, false);
        notifier.ProcessLibrary(single_streams, ilib, reads.data().read_length, single_streams.size());
    }
    else {
        io::SingleStreams single_streams(single_easy_reader(reads, true, false));
        notifier.ProcessLibrary(single_streams, ilib, reads.data().read_length, single_streams.size());
    }
}


void ProcessPairedReads(conj_graph_pack& gp, size_t ilib) {
    SequenceMapperNotifier notifier(gp);
    const io::SequencingLibrary<debruijn_config::DataSetData>& reads = cfg::get().ds.reads[ilib];
    INFO("left qauntile " << reads.data().insert_size_left_quantile << " right " << reads.data().insert_size_right_quantile);
    path_extend::SplitGraphPairInfo split_graph(
            gp, (size_t) reads.data().median_insert_size,
            (size_t) reads.data().insert_size_deviation,
            (size_t) reads.data().insert_size_left_quantile,
            (size_t) reads.data().insert_size_right_quantile,
            reads.data().read_length, gp.g.k(),
            cfg::get().pe_params.param_set.split_edge_length,
            reads.data().insert_size_distribution);

    LatePairedIndexFiller pif(gp.g, PairedReadCountWeight, gp.paired_indices[ilib]);

    SimpleLongReadMapper read_mapper(gp, gp.single_long_reads[ilib]);
    if (cfg::get().long_single_mode) {
        notifier.Subscribe(ilib, &read_mapper);
    }
    notifier.Subscribe(ilib, &split_graph);
    notifier.Subscribe(ilib, &pif);

    if (cfg::get().use_multithreading) {
        auto paired_streams = paired_binary_readers(reads, true, (size_t) reads.data().mean_insert_size);
        notifier.ProcessLibrary(paired_streams, ilib, reads.data().read_length, paired_streams.size());
        cfg::get_writable().ds.reads[ilib].data().pi_threshold = split_graph.GetThreshold();
    }
    else {
        io::PairedStreams paired_streams(paired_easy_reader(reads, true, (size_t) reads.data().mean_insert_size));
        notifier.ProcessLibrary(paired_streams, ilib, reads.data().read_length, paired_streams.size());
        cfg::get_writable().ds.reads[ilib].data().pi_threshold = split_graph.GetThreshold();
    }

    if (cfg::get().long_single_mode) {
        ProcessSingleReads(gp, ilib);
    }
}


void PairInfoCount::run(conj_graph_pack &gp, const char*) {
    if (!cfg::get().developer_mode) {
        gp.paired_indices.Attach();
        gp.paired_indices.Init();
    }

    size_t edge_length_threshold = Nx(gp.g, 50);
    for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
        INFO("Processing library #" << i);
        if (cfg::get().ds.reads[i].type() == io::LibraryType::PairedEnd ||
            cfg::get().ds.reads[i].type() == io::LibraryType::MatePairs) {

            bool insert_size_refined;
            if (cfg::get().use_multithreading) {
                auto streams = paired_binary_readers(cfg::get().ds.reads[i], false, 0);
                insert_size_refined = RefineInsertSizeForLib(gp, streams, cfg::get_writable().ds.reads[i].data(), edge_length_threshold);
            } else {
                io::PairedStreams streams(paired_easy_reader(cfg::get().ds.reads[i], false, 0));
                insert_size_refined = RefineInsertSizeForLib(gp, streams, cfg::get_writable().ds.reads[i].data(), edge_length_threshold);
            }

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
            INFO("Mapping paired reads (takes a while) ");
            ProcessPairedReads(gp, i);
        }

        if (cfg::get().ds.reads[i].type() == io::LibraryType::SingleReads && cfg::get().long_single_mode) {
            INFO("Mapping single reads (takes a while) ");
            ProcessSingleReads(gp, i);
        }

        if (cfg::get().long_single_mode) {
            INFO("Total paths obtained from single reads: " << gp.single_long_reads[i].size());
        }
    }
}

}

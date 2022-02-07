//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "paired_info_utils.hpp"
#include "is_counter.hpp"

#include "modules/alignment/sequence_mapper_notifier.hpp"
#include "paired_info/pair_info_filler.hpp"
#include "library/library.hpp"
#include "library/library_data.hpp"

#include "io/dataset_support/read_converter.hpp"
#include "utils/parallel/openmp_wrapper.h"

#include "adt/hll.hpp"

#define XXH_INLINE_ALL
#include "xxh/xxhash.h"

namespace paired_info {

using namespace debruijn_graph;

using EdgePairCounter = hll::hll_with_hasher<std::pair<EdgeId, EdgeId>>;

class EdgePairCounterFiller : public SequenceMapperListener {
    static uint64_t EdgePairHash(const std::pair<EdgeId, EdgeId> &e) {
        uint64_t h1 = e.first.hash();

        return XXH3_64bits_withSeed(&h1, sizeof(h1), e.second.hash());
    }

  public:
    EdgePairCounterFiller(size_t thread_num)
            : counter_(EdgePairHash) {
        buf_.reserve(thread_num);
        for (unsigned i = 0; i < thread_num; ++i)
          buf_.emplace_back(EdgePairHash);
    }

    void MergeBuffer(size_t i) override {
        counter_.merge(buf_[i]);
        buf_[i].clear();
    }

    void ProcessPairedRead(size_t idx,
                           const io::PairedRead&,
                           const MappingPath<EdgeId>& read1,
                           const MappingPath<EdgeId>& read2) override {
        ProcessPairedRead(buf_[idx], read1, read2);
    }
    void ProcessPairedRead(size_t idx,
                           const io::PairedReadSeq&,
                           const MappingPath<EdgeId>& read1,
                           const MappingPath<EdgeId>& read2) override {
        ProcessPairedRead(buf_[idx], read1, read2);
    }

    double cardinality() const {
        return counter_.cardinality();
    }
  private:
    void ProcessPairedRead(EdgePairCounter &buf,
                           const MappingPath<EdgeId>& path1,
                           const MappingPath<EdgeId>& path2) {
        for (size_t i = 0; i < path1.size(); ++i) {
            EdgeId edge1 = path1.edge_at(i);
            for (size_t j = 0; j < path2.size(); ++j) {
                EdgeId edge2 = path2.edge_at(j);
                buf.add({edge1, edge2});
            }
        }
    }

    std::vector<EdgePairCounter> buf_;
    EdgePairCounter counter_;
};

bool CollectLibInformation(const Graph &graph,
                           const SequenceMapperNotifier::SequenceMapperT &mapper,
                           size_t &edgepairs, SequencingLib &reads,
                           size_t edge_length_threshold) {
    INFO("Estimating insert size (takes a while)");
    InsertSizeCounter hist_counter(graph, edge_length_threshold);
    EdgePairCounterFiller pcounter(omp_get_max_threads());

    SequenceMapperNotifier notifier;
    notifier.Subscribe(&hist_counter);
    notifier.Subscribe(&pcounter);

    auto &data = reads.data();
    auto paired_streams = paired_binary_readers(reads, /*followed by rc*/false, /*insert_size*/0,
                                                /*include_merged*/true);

    notifier.ProcessLibrary(paired_streams, mapper);
    //Check read length after lib processing since mate pairs a not used until this step
    VERIFY(reads.data().unmerged_read_length != 0);

    edgepairs = size_t(pcounter.cardinality());
    INFO("Edge pairs: " << edgepairs);

    INFO(hist_counter.mapped() << " paired reads (" <<
         ((double) hist_counter.mapped() * 100.0 / (double) hist_counter.total()) <<
         "% of all) aligned to long edges");
    if (hist_counter.negative() > 3 * hist_counter.mapped())
        WARN("Too much reads aligned with negative insert size. Is the library orientation set properly?");
    if (hist_counter.mapped() == 0)
        return false;


    std::map<size_t, size_t> percentiles;
    hist_counter.FindMean(data.mean_insert_size, data.insert_size_deviation, percentiles);
    hist_counter.FindMedian(data.median_insert_size, data.insert_size_mad,
                            data.insert_size_distribution);
    if (data.median_insert_size < double(graph.k() + 2))
        return false;

    std::tie(data.insert_size_left_quantile,
             data.insert_size_right_quantile) = omnigraph::GetISInterval(0.8,
                                                                         data.insert_size_distribution);

    return !data.insert_size_distribution.empty();
}

void FillPairedIndex(const Graph &graph,
                     const SequenceMapperNotifier::SequenceMapperT &mapper,
                     SequencingLib &reads,
                     PairedIndex &index,
                     std::unique_ptr<PairedInfoFilter> filter, unsigned filter_threshold,
                     unsigned round_thr, bool use_binary) {
    const auto &data = reads.data();

    SequenceMapperNotifier notifier;
    INFO("Left insert size quantile " << data.insert_size_left_quantile <<
         ", right insert size quantile " << data.insert_size_right_quantile <<
         ", filtering threshold " << filter_threshold <<
         ", rounding threshold " << round_thr);

    LatePairedIndexFiller::WeightF weight;
    if (filter) {
        weight = [&](const std::pair<EdgeId, EdgeId> &ep,
                     const MappingRange&, const MappingRange&) {
            return (filter->lookup(ep) > filter_threshold ? 1. : 0.);
        };
    } else {
        weight = [&](const std::pair<EdgeId, EdgeId> &,
                     const MappingRange&, const MappingRange&) {
            return 1.;
        };
    }

    LatePairedIndexFiller pif(graph, weight, round_thr, index);
    notifier.Subscribe(&pif);

    if (use_binary) {
        auto paired_streams = paired_binary_readers(reads, /*followed by rc*/false, (size_t) data.mean_insert_size,
                                                    /*include merged*/true);
        notifier.ProcessLibrary(paired_streams, mapper);
    } else {
        auto paired_streams = paired_easy_readers(reads, /*followed by rc*/false,
                                                  (size_t)data.mean_insert_size, /*use_orientation*/false);
        notifier.ProcessLibrary(paired_streams, mapper);
    }
}

class DEFilter : public SequenceMapperListener {
  public:
    DEFilter(paired_info::PairedInfoFilter &filter, const Graph &g)
            : bf_(filter), g_(g) {}

    void ProcessPairedRead(size_t,
                           const io::PairedRead&,
                           const MappingPath<EdgeId>& read1,
                           const MappingPath<EdgeId>& read2) override {
        ProcessPairedRead(read1, read2);
    }
    void ProcessPairedRead(size_t,
                           const io::PairedReadSeq&,
                           const MappingPath<EdgeId>& read1,
                           const MappingPath<EdgeId>& read2) override {
        ProcessPairedRead(read1, read2);
    }
  private:
    void ProcessPairedRead(const MappingPath<EdgeId>& path1,
                           const MappingPath<EdgeId>& path2) {
        for (size_t i = 0; i < path1.size(); ++i) {
            EdgeId edge1 = path1.edge_at(i);
            for (size_t j = 0; j < path2.size(); ++j) {
                EdgeId edge2 = path2.edge_at(j);
                bf_.add({edge1, edge2});
                bf_.add({g_.conjugate(edge2), g_.conjugate(edge1)});
            }
        }
    }

    paired_info::PairedInfoFilter &bf_;
    const Graph &g_;
};

std::unique_ptr<PairedInfoFilter> FillEdgePairFilter(const Graph &graph,
                                                     const SequenceMapperNotifier::SequenceMapperT &mapper,
                                                     SequencingLib &reads,
                                                     size_t edgepairs) {
    auto filter = std::make_unique<paired_info::PairedInfoFilter>(
        [](const std::pair<EdgeId, EdgeId> &e, uint64_t seed) {
            uint64_t h1 = e.first.hash();
            return XXH3_64bits_withSeed(&h1, sizeof(h1), (e.second.hash() * seed) ^ seed);
        },
        12 * edgepairs);

    SequenceMapperNotifier notifier;
    DEFilter filter_counter(*filter, graph);
    notifier.Subscribe(&filter_counter);

    VERIFY(reads.data().unmerged_read_length != 0);
    auto stream = paired_binary_readers(reads, /*followed by rc*/false, 0, /*include merged*/true);
    notifier.ProcessLibrary(stream, mapper);

    return filter;
}

}


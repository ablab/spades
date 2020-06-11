//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "pair_info_count.hpp"

#include "assembly_graph/core/basic_graph_stats.hpp"
#include "paired_info/is_counter.hpp"
#include "paired_info/pair_info_filler.hpp"

#include "modules/alignment/long_read_mapper.hpp"
#include "modules/alignment/bwa_sequence_mapper.hpp"
#include "modules/alignment/rna/ss_coverage_filler.hpp"

#include "io/dataset_support/read_converter.hpp"

#include "adt/bf.hpp"
#include "adt/hll.hpp"

#define XXH_INLINE_ALL
#include "xxh/xxhash.h"

namespace debruijn_graph {

namespace {

using SequencingLib = io::SequencingLibrary<config::LibraryData>;
using PairedInfoFilter = bf::counting_bloom_filter<std::pair<EdgeId, EdgeId>, 2>;
using EdgePairCounter = hll::hll_with_hasher<std::pair<EdgeId, EdgeId>>;

std::shared_ptr<SequenceMapper<Graph>> ChooseProperMapper(const GraphPack& gp,
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

class DEFilter : public SequenceMapperListener {
  public:
    DEFilter(PairedInfoFilter &filter, const Graph &g)
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

    PairedInfoFilter &bf_;
    const Graph &g_;
};

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

bool CollectLibInformation(const GraphPack &gp,
                           size_t &edgepairs,
                           size_t ilib, size_t edge_length_threshold) {
    INFO("Estimating insert size (takes a while)");
    InsertSizeCounter hist_counter(gp.get<Graph>(), edge_length_threshold);
    EdgePairCounterFiller pcounter(cfg::get().max_threads);

    SequenceMapperNotifier notifier(gp, cfg::get_writable().ds.reads.lib_count());
    notifier.Subscribe(ilib, &hist_counter);
    notifier.Subscribe(ilib, &pcounter);

    SequencingLib &reads = cfg::get_writable().ds.reads[ilib];
    auto &data = reads.data();
    auto paired_streams = paired_binary_readers(reads, /*followed by rc*/false, /*insert_size*/0,
                                                /*include_merged*/true);

    notifier.ProcessLibrary(paired_streams, ilib, *ChooseProperMapper(gp, reads));
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
    if (data.median_insert_size < gp.k() + 2)
        return false;

    std::tie(data.insert_size_left_quantile,
             data.insert_size_right_quantile) = omnigraph::GetISInterval(0.8,
                                                                         data.insert_size_distribution);

    return !data.insert_size_distribution.empty();
}

size_t ProcessSingleReads(GraphPack &gp, size_t ilib,
                                 bool use_binary = true, bool map_paired = false) {
    //FIXME make const
    auto& reads = cfg::get_writable().ds.reads[ilib];
    const auto &graph = gp.get<Graph>();

    SequenceMapperNotifier notifier(gp, cfg::get_writable().ds.reads.lib_count());

    auto &single_long_reads = gp.get_mutable<LongReadContainer<Graph>>()[ilib];
    auto& trusted_paths = gp.get_mutable<path_extend::TrustedPathsContainer>()[ilib];
    LongReadMapper read_mapper(graph, single_long_reads, trusted_paths, reads.type());

    if (ShouldObtainSingleReadsPaths(ilib) || reads.is_contig_lib()) {
        //FIXME pretty awful, would be much better if listeners were shared ptrs
        notifier.Subscribe(ilib, &read_mapper);
        cfg::get_writable().ds.reads[ilib].data().single_reads_mapped = true;
    }

    SSCoverageFiller ss_coverage_filler(graph, gp.get_mutable<SSCoverageContainer>()[ilib],
                                        !cfg::get().ss.ss_enabled);
    if (cfg::get().calculate_coverage_for_each_lib) {
        INFO("Will calculate lib coverage as well");
        map_paired = true;
        notifier.Subscribe(ilib, &ss_coverage_filler);
    }

    auto mapper_ptr = ChooseProperMapper(gp, reads);
    if (use_binary) {
        auto single_streams = single_binary_readers(reads, false, map_paired);
        notifier.ProcessLibrary(single_streams, ilib, *mapper_ptr);
    } else {
        auto single_streams = single_easy_readers(reads, false,
                                                  map_paired, /*handle Ns*/false);
        notifier.ProcessLibrary(single_streams, ilib, *mapper_ptr);
    }

    return single_long_reads.size();
}


void ProcessPairedReads(GraphPack &gp,
                               std::unique_ptr<PairedInfoFilter> filter,
                               unsigned filter_threshold,
                               size_t ilib) {
    SequencingLib &reads = cfg::get_writable().ds.reads[ilib];
    const auto &data = reads.data();

    unsigned round_thr = 0;
    // Do not round if filtering is disabled
    if (filter)
        round_thr = unsigned(std::min(cfg::get().de.max_distance_coeff * data.insert_size_deviation * cfg::get().de.rounding_coeff,
                                      cfg::get().de.rounding_thr));

    SequenceMapperNotifier notifier(gp, cfg::get_writable().ds.reads.lib_count());
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

    using Indices = omnigraph::de::UnclusteredPairedInfoIndicesT<Graph>;
    LatePairedIndexFiller pif(gp.get<Graph>(), weight, round_thr, gp.get_mutable<Indices>()[ilib]);
    notifier.Subscribe(ilib, &pif);

    auto paired_streams = paired_binary_readers(reads, /*followed by rc*/false, (size_t) data.mean_insert_size,
                                                /*include merged*/true);
    notifier.ProcessLibrary(paired_streams, ilib, *ChooseProperMapper(gp, reads));
}

} // namespace

void PairInfoCount::run(GraphPack &gp, const char *) {
    gp.InitRRIndices();
    gp.EnsureBasicMapping();

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
                if (!CollectLibInformation(gp, edgepairs, i, edge_length_threshold)) {
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

                std::unique_ptr<PairedInfoFilter> filter;
                unsigned filter_threshold = cfg::get().de.raw_filter_threshold;

                // Only filter paired-end libraries
                if (filter_threshold && lib.type() == io::LibraryType::PairedEnd) {
                    filter.reset(new PairedInfoFilter([](const std::pair<EdgeId, EdgeId> &e, uint64_t seed) {
                                uint64_t h1 = e.first.hash();
                                return XXH3_64bits_withSeed(&h1, sizeof(h1), (e.second.hash() * seed) ^ seed);
                            },
                            12 * edgepairs));

                    INFO("Filtering data for library #" << i);
                    {
                        SequenceMapperNotifier notifier(gp, cfg::get_writable().ds.reads.lib_count());
                        DEFilter filter_counter(*filter, graph);
                        notifier.Subscribe(i, &filter_counter);

                        VERIFY(lib.data().unmerged_read_length != 0);
                        auto reads = paired_binary_readers(lib, /*followed by rc*/false, 0, /*include merged*/true);
                        notifier.ProcessLibrary(reads, i, *ChooseProperMapper(gp, lib));
                    }
                }

                INFO("Mapping library #" << i);
                if (lib.data().mean_insert_size != 0.0) {
                    INFO("Mapping paired reads (takes a while) ");
                    ProcessPairedReads(gp, std::move(filter), filter_threshold, i);
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

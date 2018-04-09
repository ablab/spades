//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include <paired_info/is_counter.hpp>
#include "io/dataset_support/read_converter.hpp"

#include "pair_info_count.hpp"
#include "modules/alignment/long_read_mapper.hpp"
#include "modules/alignment/bwa_sequence_mapper.hpp"
#include "paired_info/pair_info_filler.hpp"
#include "modules/alignment/rna/ss_coverage_filler.hpp"


#include "adt/bf.hpp"
#include "adt/hll.hpp"

namespace debruijn_graph {

typedef io::SequencingLibrary<config::LibraryData> SequencingLib;
using PairedInfoFilter = bf::counting_bloom_filter<std::pair<EdgeId, EdgeId>, 2>;
using EdgePairCounter = hll::hll_with_hasher<std::pair<EdgeId, EdgeId>>;

std::shared_ptr<SequenceMapper<Graph>> ChooseProperMapper(const conj_graph_pack& gp,
                                                          const SequencingLib& library) {
    if (library.type() == io::LibraryType::MatePairs) {
        INFO("Mapping mate-pairs using BWA-mem mapper");
        return std::make_shared<alignment::BWAReadMapper<Graph>>(gp.g);
    }

    if (library.data().unmerged_read_length < gp.k_value && library.type() == io::LibraryType::PairedEnd) {
        INFO("Mapping PE reads shorter than K with BWA-mem mapper");
        return std::make_shared<alignment::BWAReadMapper<Graph>>(gp.g);
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
            std::pair<EdgeId, MappingRange> mapping_edge_1 = path1[i];
            for (size_t j = 0; j < path2.size(); ++j) {
                std::pair<EdgeId, MappingRange> mapping_edge_2 = path2[j];
                bf_.add({mapping_edge_1.first, mapping_edge_2.first});
                bf_.add({g_.conjugate(mapping_edge_2.first), g_.conjugate(mapping_edge_1.first)});
            }
        }
    }

    PairedInfoFilter &bf_;
    const Graph &g_;
};

class EdgePairCounterFiller : public SequenceMapperListener {
    static uint64_t EdgePairHash(const std::pair<EdgeId, EdgeId> &e) {
        uint64_t h1 = e.first.hash();
        return CityHash64WithSeeds((const char*)&h1, sizeof(h1), e.second.hash(), 0x0BADF00D);
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

    std::pair<double, bool> cardinality() const {
        return counter_.cardinality();
    }
  private:
    void ProcessPairedRead(EdgePairCounter &buf,
                           const MappingPath<EdgeId>& path1,
                           const MappingPath<EdgeId>& path2) {
        for (size_t i = 0; i < path1.size(); ++i) {
            std::pair<EdgeId, MappingRange> mapping_edge_1 = path1[i];
            for (size_t j = 0; j < path2.size(); ++j) {
                std::pair<EdgeId, MappingRange> mapping_edge_2 = path2[j];
                buf.add({mapping_edge_1.first, mapping_edge_2.first});
            }
        }
    }

    std::vector<EdgePairCounter> buf_;
    EdgePairCounter counter_;
};

static bool HasGoodRRLibs() {
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

static bool HasOnlyMP() {
    for (const auto &lib : cfg::get().ds.reads) {
        if (lib.type() == io::LibraryType::PathExtendContigs)
            continue;

        if (lib.type() != io::LibraryType::MatePairs &&
            lib.type() != io::LibraryType::HQMatePairs)
            return false;
    }

    return true;
}

static bool ShouldObtainLibCoverage() {
    return cfg::get().calculate_coverage_for_each_lib;
}

//todo improve logic
static bool ShouldObtainSingleReadsPaths(size_t ilib) {
    using config::single_read_resolving_mode;
    switch (cfg::get().single_reads_rr) {
        case single_read_resolving_mode::all:
            return true;
        case single_read_resolving_mode::only_single_libs:
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
        case single_read_resolving_mode::none:
            break;
        default:
            VERIFY_MSG(false, "Invalid mode value");
    }
    return false;
}

static bool CollectLibInformation(const conj_graph_pack &gp,
                                  size_t &edgepairs,
                                  size_t ilib, size_t edge_length_threshold) {
    INFO("Estimating insert size (takes a while)");
    InsertSizeCounter hist_counter(gp, edge_length_threshold);
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

    auto pres = pcounter.cardinality();
    edgepairs = (!pres.second ? 64ull * 1024 * 1024 : size_t(pres.first));
    INFO("Edge pairs: " << edgepairs << (!pres.second ? " (rough upper limit)" : ""));

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
    if (data.median_insert_size < gp.k_value + 2)
        return false;

    std::tie(data.insert_size_left_quantile,
             data.insert_size_right_quantile) = omnigraph::GetISInterval(0.8,
                                                                         data.insert_size_distribution);

    return !data.insert_size_distribution.empty();
}

// FIXME: This needs to be static
static void ProcessSingleReads(conj_graph_pack &gp,
                        size_t ilib,
                        bool use_binary = true,
                        bool map_paired = false) {
    //FIXME make const
    auto& reads = cfg::get_writable().ds.reads[ilib];

    SequenceMapperNotifier notifier(gp, cfg::get_writable().ds.reads.lib_count());

    LongReadMapper read_mapper(gp.g, gp.single_long_reads[ilib],
                               ChooseProperReadPathExtractor(gp.g, reads.type()));

    if (ShouldObtainSingleReadsPaths(ilib) || reads.is_contig_lib()) {
        //FIXME pretty awful, would be much better if listeners were shared ptrs
        notifier.Subscribe(ilib, &read_mapper);
        cfg::get_writable().ds.reads[ilib].data().single_reads_mapped = true;
    }

    SSCoverageFiller ss_coverage_filler(gp.g, gp.ss_coverage[ilib], !cfg::get().ss.ss_enabled);
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
}


static void ProcessPairedReads(conj_graph_pack &gp,
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

    LatePairedIndexFiller pif(gp.g,
                              weight, round_thr,
                              gp.paired_indices[ilib]);
    notifier.Subscribe(ilib, &pif);

    auto paired_streams = paired_binary_readers(reads, /*followed by rc*/false, (size_t) data.mean_insert_size,
                                                /*include merged*/true);
    notifier.ProcessLibrary(paired_streams, ilib, *ChooseProperMapper(gp, reads));
}

void PairInfoCount::run(conj_graph_pack &gp, const char *) {
    gp.InitRRIndices();
    gp.EnsureBasicMapping();

    //TODO implement better universal logic
    size_t edge_length_threshold = cfg::get().mode == config::pipeline_type::meta ? 900 : stats::Nx(gp.g, 50);
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
                                return CityHash64WithSeeds((const char*)&h1, sizeof(h1), e.second.hash(), seed);
                            },
                            12 * edgepairs));

                    INFO("Filtering data for library #" << i);
                    {
                        SequenceMapperNotifier notifier(gp, cfg::get_writable().ds.reads.lib_count());
                        DEFilter filter_counter(*filter, gp.g);
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
                ProcessSingleReads(gp, i, /*use_binary*/true, /*map_paired*/true);
                INFO("Total paths obtained from single reads: " << gp.single_long_reads[i].size());
            }
        }
    }
}

}

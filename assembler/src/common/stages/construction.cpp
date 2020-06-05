//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "construction.hpp"

#include "assembly_graph/construction/early_simplification.hpp"
#include "modules/alignment/edge_index.hpp"
#include "modules/graph_construction.hpp"

#include "pipeline/graph_pack.hpp"
#include "pipeline/genomic_info.hpp"

#include "io/dataset_support/dataset_readers.hpp"
#include "io/dataset_support/read_converter.hpp"
#include "io/reads/coverage_filtering_read_wrapper.hpp"
#include "io/reads/multifile_reader.hpp"

#include "utils/filesystem/temporary.hpp"
#include "utils/ph_map/coverage_hash_map_builder.hpp"


namespace debruijn_graph {

struct ConstructionStorage {
    using CoverageMap = utils::PerfectHashMap<RtSeq, uint32_t, utils::slim_kmer_index_traits<RtSeq>, utils::DefaultStoring>;

    ConstructionStorage(unsigned k)
            : ext_index(k) {}

    utils::DeBruijnExtensionIndex<> ext_index;

    std::unique_ptr<qf::cqf> cqf;
    std::unique_ptr<kmers::KMerDiskStorage<RtSeq>> kmers;
    std::unique_ptr<CoverageMap> coverage_map;
    config::debruijn_config::construction params;
    io::ReadStreamList<io::SingleReadSeq> read_streams;
    io::ReadStreamList<io::SingleReadSeq> contigs_streams;
    fs::TmpDir workdir;
};

bool add_trusted_contigs(io::DataSet<config::LibraryData> &libraries,
                       io::ReadStreamList<io::SingleReadSeq> &trusted_list) {
    std::vector<size_t> trusted_contigs;
    for (size_t i = 0; i < libraries.lib_count(); ++i) {
        auto& lib = libraries[i];
        if (lib.type() != io::LibraryType::TrustedContigs)
            continue;
        trusted_contigs.push_back(i);
    }

    if (!trusted_contigs.empty()) {
        trusted_list = io::single_binary_readers_for_libs(libraries, trusted_contigs, true, false);
    }
    return !trusted_contigs.empty();
}

void merge_read_streams(io::ReadStreamList<io::SingleReadSeq> &streams1,
                        io::ReadStreamList<io::SingleReadSeq> &streams2) {
    for (size_t i = 0; i < streams2.size(); ++i) {
        if (i < streams1.size()) {
            streams1[i] = io::MultifileWrap<io::SingleReadSeq>(std::move(streams1[i]), std::move(streams2[i]));
        } else {
            streams1.push_back(std::move(streams2[i]));
        }
    }
}

io::ReadStreamList<io::SingleReadSeq> temp_merge_read_streams(io::ReadStreamList<io::SingleReadSeq> &streams1,
                                                              io::ReadStreamList<io::SingleReadSeq> &streams2) {
    io::ReadStreamList<io::SingleReadSeq> merge_stream_list;

    for (size_t i = 0; i < std::max(streams1.size(), streams2.size()); ++i) {
        if (i < streams1.size() && i < streams2.size()) {
            merge_stream_list.push_back(io::ScopedMultifileWrap<io::SingleReadSeq>(streams1[i], streams2[i]));
        } else if (i < streams1.size()) {
            merge_stream_list.push_back(io::ScopedMultifileWrap<io::SingleReadSeq>(streams1[i]));
        } else {
            merge_stream_list.push_back(io::ScopedMultifileWrap<io::SingleReadSeq>(streams2[i]));
        }
    }

    return merge_stream_list;
}




void add_additional_contigs_to_lib(std::string path_to_additional_contigs_dir, size_t max_threads,
                                   io::ReadStreamList<io::SingleReadSeq> &trusted_list) {
    io::SequencingLibraryT seq_lib;
    seq_lib.set_type(io::LibraryType::TrustedContigs);
    seq_lib.set_orientation(io::LibraryOrientation::Undefined);
    seq_lib.data().lib_index = size_t(-1);
    auto& bin_info = seq_lib.data().binary_reads_info;
    bin_info.single_read_prefix = path_to_additional_contigs_dir + "/contigs";
    bin_info.bin_reads_info_file = path_to_additional_contigs_dir + "/contigs_info";
    bin_info.binary_converted = true;
    bin_info.chunk_num = max_threads;

    io::ReadStreamList<io::SingleReadSeq> lib_streams = io::single_binary_readers(seq_lib, true, false);
    merge_read_streams(trusted_list, lib_streams);
}

void Construction::init(debruijn_graph::GraphPack &gp, const char *) {
    init_storage(unsigned(gp.k()));

    auto& dataset = cfg::get_writable().ds;

    // Has to be separate stream for not counting it in coverage
    if (add_trusted_contigs(dataset.reads, storage().contigs_streams))
        INFO("Trusted contigs will be used in graph construction");

    if (cfg::get().use_additional_contigs) {
        INFO("Contigs from previous K will be used: " << cfg::get().additional_contigs);
        add_additional_contigs_to_lib(cfg::get().additional_contigs, cfg::get().max_threads, storage().contigs_streams);
    }

    // FIXME: indices here are awful
    std::vector<size_t> libs_for_construction;
    for (size_t i = 0; i < dataset.reads.lib_count(); ++i) {
        if (dataset.reads[i].is_graph_constructable()) {
            libs_for_construction.push_back(i);
        }
    }

    storage().params = cfg::get().con;
    storage().workdir = fs::tmp::make_temp_dir(gp.workdir(), "construction");
    //FIXME needs to be changed if we move to hash only filtering
    storage().read_streams = io::single_binary_readers_for_libs(dataset.reads, libs_for_construction);

    //Updating dataset stats
    VERIFY(dataset.RL == 0 && dataset.aRL == 0.);
    size_t merged_max_len = 0;
    uint64_t total_nucls = 0;
    size_t read_count = 0;
    for (size_t lib_id : libs_for_construction) {
        auto lib_data = dataset.reads[lib_id].data();
        if (lib_data.unmerged_read_length == 0) {
            FATAL_ERROR("Failed to determine read length for library #" << lib_data.lib_index << ". "
                    "Check that not only merged reads are present.");
        }
        dataset.no_merge_RL = std::max(dataset.no_merge_RL, lib_data.unmerged_read_length);
        merged_max_len = std::max(merged_max_len, lib_data.merged_read_length);
        total_nucls += dataset.reads[lib_id].data().total_nucls;
        read_count += dataset.reads[lib_id].data().read_count;
    }

    dataset.RL = std::max(dataset.no_merge_RL, merged_max_len);
    INFO("Max read length " << dataset.RL);

    if (merged_max_len > 0)
        INFO("Max read length without merged " << dataset.no_merge_RL);

    dataset.aRL = double(total_nucls) / double(read_count);
    INFO("Average read length " << dataset.aRL);
}

void Construction::fini(debruijn_graph::GraphPack &) {
    reset_storage();
}

Construction::~Construction() {}

namespace {

class CoverageFilter: public Construction::Phase {
  public:
    CoverageFilter()
            : Construction::Phase("k-mer multiplicity estimation", "cqf_filter") { }
    virtual ~CoverageFilter() = default;

    void run(debruijn_graph::GraphPack &, const char*) override {
        auto &read_streams = storage().read_streams;
        const auto &index = storage().ext_index;
        using storing_type = decltype(storage().ext_index)::storing_type;

        VERIFY_MSG(read_streams.size(), "No input streams specified");

        unsigned rthr = storage().params.read_cov_threshold;

        using KmerFilter = utils::StoringTypeFilter<storing_type>;

        unsigned kplusone = index.k() + 1;
        rolling_hash::SymmetricCyclicHash<rolling_hash::NDNASeqHash> hasher(kplusone);

        INFO("Estimating k-mers cardinality");
        size_t kmers = EstimateCardinalityUpperBound(kplusone, read_streams, hasher, KmerFilter());

        // Create main CQF using # of slots derived from estimated # of k-mers
        storage().cqf.reset(new qf::cqf(kmers));

        INFO("Building k-mer coverage histogram");
        FillCoverageHistogram(*storage().cqf, kplusone, hasher, read_streams, rthr, KmerFilter());

        // Replace input streams with wrapper ones
        storage().read_streams = io::CovFilteringWrap(std::move(read_streams), kplusone, hasher, *storage().cqf, rthr);
    }

    void load(debruijn_graph::GraphPack&,
              const std::string &,
              const char*) override {
        VERIFY_MSG(false, "implement me");
    }

    void save(const debruijn_graph::GraphPack&,
              const std::string &,
              const char*) const override {
        // VERIFY_MSG(false, "implement me");
    }

};


class KMerCounting : public Construction::Phase {
    typedef rolling_hash::SymmetricCyclicHash<> SeqHasher;
public:
    KMerCounting()
            : Construction::Phase("k+1-mer counting", "kpomer_counting") { }

    virtual ~KMerCounting() = default;

    void run(debruijn_graph::GraphPack &, const char*) override {
        auto &read_streams = storage().read_streams;
        auto &contigs_streams = storage().contigs_streams;
        const auto &index = storage().ext_index;
        size_t buffer_size = storage().params.read_buffer_size;
        using storing_type = decltype(storage().ext_index)::storing_type;

        VERIFY_MSG(read_streams.size(), "No input streams specified");


        io::ReadStreamList<io::SingleReadSeq> merge_streams = temp_merge_read_streams(read_streams, contigs_streams);

        unsigned nthreads = (unsigned)merge_streams.size();
        using Splitter =  utils::DeBruijnReadKMerSplitter<io::SingleReadSeq,
                                                          utils::StoringTypeFilter<storing_type>>;

        kmers::KMerDiskCounter<RtSeq>
                counter(storage().workdir,
                        Splitter(storage().workdir, index.k() + 1, merge_streams, buffer_size));
        auto kmers = counter.Count(10 * nthreads, nthreads);
        storage().kmers.reset(new kmers::KMerDiskStorage<RtSeq>(std::move(kmers)));
    }

    void load(debruijn_graph::GraphPack&,
              const std::string &,
              const char*) override {
        VERIFY_MSG(false, "implement me");
    }

    void save(const debruijn_graph::GraphPack&,
              const std::string &,
              const char*) const override {
        // VERIFY_MSG(false, "implement me");
    }
};

class ExtensionIndexBuilder : public Construction::Phase {
public:
    ExtensionIndexBuilder()
            : Construction::Phase("Extension index construction", "extension_index_construction") { }

    virtual ~ExtensionIndexBuilder() = default;

    void run(debruijn_graph::GraphPack &, const char*) override {
        // FIXME: We just need files here, not the full counter. Implement refererence counting scheme!
        utils::DeBruijnExtensionIndexBuilder().BuildExtensionIndexFromKPOMers(storage().workdir,
                                                                              storage().ext_index,
                                                                              *storage().kmers,
                                                                              unsigned(storage().read_streams.size()),
                                                                              storage().params.read_buffer_size);
    }

    void load(debruijn_graph::GraphPack&,
              const std::string &,
              const char*) override {
        VERIFY_MSG(false, "implement me");
    }

    void save(const debruijn_graph::GraphPack&,
              const std::string &,
              const char*) const override {
        // VERIFY_MSG(false, "implement me");
    }
};


class EarlyTipClipper : public Construction::Phase {
public:
    EarlyTipClipper()
            : Construction::Phase("Early tip clipping", "early_tip_clipper") { }

    virtual ~EarlyTipClipper() = default;

    void run(debruijn_graph::GraphPack &gp, const char*) override {
        if (!storage().params.early_tc.length_bound) {
            INFO("Early tip clipper length bound set as (RL - K)");
            storage().params.early_tc.length_bound.reset(cfg::get().ds.RL - gp.k());
        }
        EarlyTipClipperProcessor(storage().ext_index, *storage().params.early_tc.length_bound).ClipTips();
    }

    void load(debruijn_graph::GraphPack&,
              const std::string &,
              const char*) override {
        VERIFY_MSG(false, "implement me");
    }

    void save(const debruijn_graph::GraphPack&,
              const std::string &,
              const char*) const override {
        // VERIFY_MSG(false, "implement me");
    }
};

class EarlyATClipper : public Construction::Phase {
public:
    EarlyATClipper()
            : Construction::Phase("Early A/T remover", "early_at_remover") { }

    virtual ~EarlyATClipper() = default;

    void run(debruijn_graph::GraphPack &, const char*) override {
        EarlyLowComplexityClipperProcessor at_processor(storage().ext_index, 0.8, 10, 200);
        at_processor.RemoveATEdges();
        at_processor.RemoveATTips();
    }

    void load(debruijn_graph::GraphPack&,
              const std::string &,
              const char*) override {
        VERIFY_MSG(false, "implement me");
    }

    void save(const debruijn_graph::GraphPack&,
              const std::string &,
              const char*) const override {
        // VERIFY_MSG(false, "implement me");
    }
};

class GraphCondenser : public Construction::Phase {
public:
    GraphCondenser()
            : Construction::Phase("Condensing graph", "graph_condensing") { }

    virtual ~GraphCondenser() = default;

    void run(debruijn_graph::GraphPack &gp, const char*) override {
        auto &index = gp.get_mutable<EdgeIndex<Graph>>();
        if (index.IsAttached())
            index.Detach();
        DeBruijnGraphExtentionConstructor<Graph>(gp.get_mutable<Graph>(), storage().ext_index).ConstructGraph(storage().params.keep_perfect_loops);
    }

    void load(debruijn_graph::GraphPack&,
              const std::string &,
              const char*) override {
        VERIFY_MSG(false, "implement me");
    }

    void save(const debruijn_graph::GraphPack&,
              const std::string &,
              const char*) const override {
        //FIXME why commented here and others
        // VERIFY_MSG(false, "implement me");
    }
};

//FIXME unused?
class EdgeIndexFiller : public Construction::Phase {
public:
    EdgeIndexFiller()
            : Construction::Phase("Edge index filling", "initial_edge_index_filling") { }

    virtual ~EdgeIndexFiller() = default;

    void run(debruijn_graph::GraphPack &gp, const char*) override {
        auto &index = gp.get_mutable<EdgeIndex<Graph>>();
        index.Refill();
        index.Attach();
    }

    void load(debruijn_graph::GraphPack&,
              const std::string &,
              const char*) override {
        VERIFY_MSG(false, "implement me");
    }

    void save(const debruijn_graph::GraphPack&,
              const std::string &,
              const char*) const override {
        // VERIFY_MSG(false, "implement me");
    }
};

class PHMCoverageFiller : public Construction::Phase {
public:
    PHMCoverageFiller()
            : Construction::Phase("Filling coverage indices (PHM)", "coverage_filling_phm") {}
    virtual ~PHMCoverageFiller() = default;

    void run(debruijn_graph::GraphPack &gp, const char *) override {
        storage().coverage_map.reset(new ConstructionStorage::CoverageMap(storage().kmers->k()));
        auto &coverage_map = *storage().coverage_map;

        utils::CoverageHashMapBuilder().BuildIndex(coverage_map,
                                                   *storage().kmers,
                                                   storage().read_streams);
        /*
        INFO("Checking the PHM");

        auto &index = gp.index.inner_index();
        for (auto I = index.value_cbegin(), E = index.value_cend();
             I != E; ++I) {
            const auto& edge_info = *I;

            Sequence sk = gp.g.EdgeNucls(edge_info.edge_id).Subseq(edge_info.offset, edge_info.offset + index.k());
            auto kwh = coverage_map.ConstructKWH(sk.start<RtSeq>(index.k()));

            uint32_t cov = coverage_map.get_value(kwh, utils::InvertableStoring::trivial_inverter<uint32_t>());
            if (edge_info.count != cov)
                INFO("" << kwh << ":" << edge_info.count << ":" << cov);
        } */

        INFO("Filling coverage and flanking coverage from PHM");
        FillCoverageAndFlankingFromPHM(coverage_map,
                                       gp.get_mutable<Graph>(), gp.get_mutable<omnigraph::FlankingCoverage<Graph>>());

        std::vector<size_t> hist;
        size_t maxcov = 0;
        size_t kmer_per_record = 1;
        if (EdgeIndex<Graph>::IsInvertable())
            kmer_per_record = 2;

        for (auto I = coverage_map.value_cbegin(), E = coverage_map.value_cend(); I != E; ++I) {
            size_t ccov = *I;
            if (!ccov)
                continue;
            maxcov = std::max(ccov, maxcov);
            if (maxcov > hist.size())
                hist.resize(maxcov, 0);
            hist[ccov - 1] += kmer_per_record;
        }

        gp.get_mutable<GenomicInfo>().set_cov_histogram(hist);
    }

    void load(debruijn_graph::GraphPack&,
              const std::string &,
              const char*) override {
        VERIFY_MSG(false, "implement me");
    }

    void save(const debruijn_graph::GraphPack&,
              const std::string &,
              const char*) const override {
        // VERIFY_MSG(false, "implement me");
    }

};

} // namespace

Construction::Construction()
        : spades::CompositeStageDeferred<ConstructionStorage>("de Bruijn graph construction", "construction") {
    if (cfg::get().con.read_cov_threshold)
        add<CoverageFilter>();

    add<KMerCounting>();

    add<ExtensionIndexBuilder>();
    if (config::PipelineHelper::IsRNAPipeline(cfg::get().mode))
        add<EarlyATClipper>();
    if (cfg::get().con.early_tc.enable && !cfg::get().gap_closer_enable)
        add<EarlyTipClipper>();
    add<GraphCondenser>();
    add<PHMCoverageFiller>();
}


} //namespace debruijn_graph

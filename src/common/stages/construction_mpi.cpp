//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "construction_mpi.hpp"

#include "alignment/edge_index.hpp"
#include "assembly_graph/construction/early_simplification.hpp"
#include "assembly_graph/construction/debruijn_graph_constructor_mpi.hpp"
#include "io/dataset_support/dataset_readers.hpp"
#include "io/dataset_support/read_converter.hpp"
#include "io/reads/coverage_filtering_read_wrapper.hpp"
#include "io/reads/multifile_reader.hpp"
#include "kmer_index/ph_map/coverage_hash_map_builder.hpp"
#include "kmer_index/extension_index/kmer_extension_index_builder_mpi.hpp"
#include "modules/graph_construction.hpp"
#include "pipeline/genomic_info.hpp"
#include "pipeline/graph_pack.hpp"
#include "pipeline/mpi_stage.hpp"
#include "pipeline/partask_mpi.hpp"
#include "utils/filesystem/temporary.hpp"

namespace debruijn_graph {

struct ConstructionStorage {
    using CoverageMap = kmers::PerfectHashMap<RtSeq, uint32_t, kmers::slim_kmer_index_traits<RtSeq>, kmers::DefaultStoring>;

    ConstructionStorage(unsigned k)
            : ext_index(k) {}

    kmers::DeBruijnExtensionIndex<> ext_index;

    std::unique_ptr<qf::cqf> cqf;
    std::unique_ptr<kmers::KMerDiskStorage<RtSeq>> kmers;
    std::unique_ptr<CoverageMap> coverage_map;
    config::debruijn_config::construction params;
    io::ReadStreamList<io::SingleReadSeq> read_streams;
    io::ReadStreamList<io::SingleReadSeq> contigs_streams;
    fs::TmpDir workdir;
};

static bool add_trusted_contigs(io::DataSet<config::LibraryData> &libraries,
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

static void merge_read_streams(io::ReadStreamList<io::SingleReadSeq> &streams1,
                        io::ReadStreamList<io::SingleReadSeq> &streams2) {
    for (size_t i = 0; i < streams2.size(); ++i) {
        if (i < streams1.size()) {
            streams1[i] = io::MultifileWrap<io::SingleReadSeq>(std::move(streams1[i]), std::move(streams2[i]));
        } else {
            streams1.push_back(std::move(streams2[i]));
        }
    }
}

static io::ReadStreamList<io::SingleReadSeq> temp_merge_read_streams(io::ReadStreamList<io::SingleReadSeq> &streams1,
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




static void add_additional_contigs_to_lib(std::filesystem::path path_to_additional_contigs_dir, size_t max_threads,
                                   io::ReadStreamList<io::SingleReadSeq> &trusted_list) {
    io::SequencingLibraryT seq_lib;
    seq_lib.set_type(io::LibraryType::TrustedContigs);
    seq_lib.set_orientation(io::LibraryOrientation::Undefined);
    seq_lib.data().lib_index = size_t(-1);
    auto& bin_info = seq_lib.data().binary_reads_info;
    bin_info.single_read_prefix = path_to_additional_contigs_dir / "contigs";
    bin_info.bin_reads_info_file = path_to_additional_contigs_dir / "contigs_info";
    bin_info.binary_converted = true;
    bin_info.chunk_num = max_threads;

    io::ReadStreamList<io::SingleReadSeq> lib_streams = io::single_binary_readers(seq_lib, true, false);
    merge_read_streams(trusted_list, lib_streams);
}

void ConstructionMPI::init(graph_pack::GraphPack &gp, const char *) {
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
    size_t num_readers = partask::overall_num_threads();
    storage().read_streams = io::single_binary_readers_for_libs(dataset.reads, libs_for_construction, true, true, num_readers);
    INFO("Overall number of readers (actual): " << storage().read_streams.size());

    // Updating dataset stats
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

void ConstructionMPI::fini(graph_pack::GraphPack &) {
    reset_storage();
}

ConstructionMPI::~ConstructionMPI() {}

namespace {

class CoverageFilter: public ConstructionMPI::Phase {
  public:
    CoverageFilter()
            : ConstructionMPI::Phase("k-mer multiplicity estimation", "cqf_filter") { }
    virtual ~CoverageFilter() = default;

    void run(graph_pack::GraphPack &, const char*) override {
        auto &read_streams = storage().read_streams;
        const auto &index = storage().ext_index;
        using storing_type = decltype(storage().ext_index)::storing_type;

        VERIFY_MSG(read_streams.size(), "No input streams specified");

        unsigned rthr = storage().params.read_cov_threshold;

        using KmerFilter = kmers::StoringTypeFilter<storing_type>;

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

    void load(graph_pack::GraphPack&,
              const std::filesystem::path &,
              const char*) override {
        VERIFY_MSG(false, "implement me");
    }

void save(const graph_pack::GraphPack&,
              const std::filesystem::path &,
              const char*) const override {
        // VERIFY_MSG(false, "implement me");
    }

};


template <typename Seq, typename storing_type>
class ReadKMerCountingTask {
 private:
    ReadKMerCountingTask() = default;
    using ReadStreams = io::ReadStreamList<io::SingleReadSeq>;
 public:
    ReadKMerCountingTask(const std::string &dir, unsigned k, size_t buffer_size,
                         unsigned num_buckets, unsigned num_threads)
        : dir_{dir}, k_{k}, buffer_size_{buffer_size}, num_buckets_{num_buckets}, num_threads_{num_threads} {};
    ReadKMerCountingTask(std::istream &is) { deserialize(is); }
    std::ostream &serialize(std::ostream &os) const {
        io::binary::BinWrite(os, dir_, k_, buffer_size_, num_buckets_, num_threads_);
        return os;
    }

    std::istream &deserialize(std::istream &is) {
        io::binary::BinRead(is, dir_, k_, buffer_size_, num_buckets_, num_threads_);
        return is;
    }

    template<typename... Args>
    auto make_splitter(size_t, ReadStreams &read_streams, Args &&...) {
        return partask::make_seq_along_generator(read_streams);
    }

    void process(std::istream &is, std::ostream &os, ReadStreams &read_streams) {
        auto chunks = partask::get_seq(is);

        ReadStreams streams = partask::create_empty_stream_list<ReadStreams>(chunks.size());
        partask::swap_streams(read_streams, streams, chunks);
        streams.reset();

        if (streams.size() == 0) {
            return;
        }

        auto workdir = fs::tmp::acquire_temp_dir(dir_);
        workdir->release();
        kmers::DeBruijnReadKMerSplitter<io::SingleReadSeq, kmers::StoringTypeFilter<storing_type>> splitter(
            workdir, k_, streams, buffer_size_);

        kmers::KMerDiskCounter<RtSeq> counter(workdir, splitter);

        auto kmerstorage = counter.CountAll(num_buckets_, num_threads_, /* merge */ false);
        INFO("k-mers counted successfully");

        if (kmerstorage.total_kmers() == 0) {
            WARN("No kmers were extracted from reads. Check the read lengths and k-mer length settings");
        }

        kmerstorage.BinWrite(os);
        kmerstorage.release_all();

        partask::swap_streams(read_streams, streams, chunks);
    }

    auto merge(const std::vector<std::istream *> &piss, ReadStreams & /*read_streams*/) {
        auto workdir = fs::tmp::acquire_temp_dir(dir_);
        workdir->release();


        std::vector<kmers::KMerDiskStorage<Seq>> storages;
        for (size_t i = 0; i < piss.size(); ++i) {
            auto &is = *piss[i];
            kmers::KMerDiskStorage<Seq> kmerstorage(workdir, k_, typename kmer::KMerSegmentPolicy<Seq>(num_buckets_));
            kmerstorage.BinRead(is);
            storages.push_back(std::move(kmerstorage));
        }

        return storages;
    }

 private:
    std::string dir_;
    unsigned k_;
    size_t buffer_size_;
    unsigned num_buckets_;
    unsigned num_threads_;
};


class KMerCounting : public ConstructionMPI::Phase {
    typedef rolling_hash::SymmetricCyclicHash<> SeqHasher;
public:
    KMerCounting()
            : ConstructionMPI::Phase("k+1-mer counting", "kpomer_counting") { }

    virtual ~KMerCounting() = default;

    bool distributed() const override { return true; }

    void run(graph_pack::GraphPack &, const char*) override {
        sync();  // TODO change for syncfs (it requires opened file descriptor from the fs)
        auto &read_streams = storage().read_streams;
        auto &contigs_streams = storage().contigs_streams;
        const auto &index = storage().ext_index;
        size_t buffer_size = storage().params.read_buffer_size;
        using storing_type = decltype(storage().ext_index)::storing_type;

        VERIFY_MSG(read_streams.size(), "No input streams specified");

        io::ReadStreamList<io::SingleReadSeq> merge_streams = temp_merge_read_streams(read_streams, contigs_streams);

        unsigned k = index.k();
        unsigned nthreads = cfg::get().max_threads;
        unsigned num_buckets = 10 * nthreads;
        using Seq = RtSeq;

        kmers::KMerDiskStorage<Seq> kmerfiles2(storage().workdir, k + 1, typename kmer::KMerSegmentPolicy<Seq>(num_buckets));


        VERIFY(partask::all_equal(num_buckets));
        VERIFY(partask::all_equal(k));
        VERIFY(partask::all_equal(nthreads));

        partask::TaskRegistry treg;
        auto merge_kmer_files = treg.add<kmers::MergeKMerFilesTask<Seq>>();
        auto kmercount = treg.add<ReadKMerCountingTask<Seq, storing_type>>(std::ref(merge_streams));
        treg.listen();
        INFO("Listening started " << (partask::master() ? "(master)" : "(worker)"));

        if (partask::master()) {
            INFO("Start KMer Counting")
            auto unmerged_kmerfiles = kmercount(storage().workdir->dir(), k + 1, buffer_size, num_buckets, nthreads);

            std::vector<std::string> outputfiles;
            for (size_t i = 0; i < kmerfiles2.num_buckets(); ++i) {
                outputfiles.push_back(kmerfiles2.create(i)->file());
            }

            INFO("Start Merge results from different nodes")
            merge_kmer_files(std::move(unmerged_kmerfiles), outputfiles);
        }
        treg.stop_listening();
        partask::broadcast(kmerfiles2);

        size_t kmers = kmerfiles2.total_kmers();

        if (!kmers) {
            FATAL_ERROR("No kmers were extracted from reads. Check the read lengths and k-mer length settings");
        }
        INFO(kmers << " k+1-mers (k=" << k << ") were extracted");

        storage().kmers.reset(new kmers::KMerDiskStorage<RtSeq>(std::move(kmerfiles2)));
    }

    void load(graph_pack::GraphPack&,
              const std::filesystem::path &,
              const char*) override {
        VERIFY_MSG(false, "implement me");
    }

    void save(const graph_pack::GraphPack&,
              const std::filesystem::path &,
              const char*) const override {
        // VERIFY_MSG(false, "implement me");
    }
};

class ExtensionIndexBuilder : public ConstructionMPI::Phase {
public:
    ExtensionIndexBuilder()
            : ConstructionMPI::Phase("Extension index construction", "extension_index_construction") { }

    virtual ~ExtensionIndexBuilder() = default;

    bool distributed() const override { return true; }

    void run(graph_pack::GraphPack &, const char*) override {
        // FIXME: We just need files here, not the full counter. Implement refererence counting scheme!
        if (!storage().kmers)
            storage().kmers.reset(new kmers::KMerDiskStorage<RtSeq>());
        partask::broadcast(*storage().kmers);

        kmers::DeBruijnExtensionIndexBuilderMPI().BuildExtensionIndexFromKPOMersMPI(storage().workdir,
                                                                                    storage().ext_index,
                                                                                    *storage().kmers,
                                                                                    cfg::get().max_threads,
                                                                                    storage().params.read_buffer_size);
    }

    void load(graph_pack::GraphPack&,
              const std::filesystem::path &,
              const char*) override {
        VERIFY_MSG(false, "implement me");
    }

    void save(const graph_pack::GraphPack&,
              const std::filesystem::path &,
              const char*) const override {
        // VERIFY_MSG(false, "implement me");
    }
};


template <typename Index>
class TipClippingTask {
    TipClippingTask() = default;
 public:
    TipClippingTask(size_t length_bound) : length_bound_{length_bound} {}
    TipClippingTask(std::istream &is) { deserialize(is); }
    std::ostream &serialize(std::ostream &os) const {
        io::binary::BinWrite(os, length_bound_);
        return os;
    }

    std::istream &deserialize(std::istream &is) {
        io::binary::BinRead(is, length_bound_);
        return is;
    }

    auto make_splitter(size_t size, Index &) {
        return partask::make_seq_plus_n_generator(size);
    }

    void process(std::istream &is, std::ostream &os, Index &index) {
        size_t n = 0;
        std::vector<size_t> chunks = partask::get_seq_plus_n(is, n);

        INFO("Job got, " << chunks.size() << "/" << n << "chunks");
        auto iters = index.kmer_begin(n);
        std::vector<typename Index::kmer_iterator> local_iters;
        for (size_t i : chunks) {
            if (i < iters.size()) {
                local_iters.push_back(std::move(iters[i]));
            }
        }
        size_t kpo_mers_removed = EarlyTipClipperProcessor(index, length_bound_).ClipTips(local_iters);  // TODO support empty input

        INFO("K+1-mers removed: " << kpo_mers_removed);
        partask::allreduce(index.raw_data(), index.raw_size(), MPI_BAND);
        io::binary::BinWrite(os, kpo_mers_removed);
    }

    size_t merge(const std::vector<std::istream *> &piss, Index&) {
        size_t kpo_mers_removed = 0;
        for (auto &pis : piss) {
            kpo_mers_removed += io::binary::BinRead<size_t>(*pis);
        }
        return kpo_mers_removed;
    }

 private:
    size_t length_bound_;
};


class EarlyTipClipper : public ConstructionMPI::Phase {
public:
    EarlyTipClipper()
            : ConstructionMPI::Phase("Early tip clipping (MPI)", "early_tip_clipper_mpi") { }

    virtual ~EarlyTipClipper() = default;

    bool distributed() const override { return true; }

    void run(graph_pack::GraphPack &gp, const char*) override {
        partask::TaskRegistry treg;
        auto &index = storage().ext_index;
        using Index = std::remove_reference_t<decltype(index)>;
        VERIFY(partask::all_equal(index.size()));

        auto clip_tips = treg.add<TipClippingTask<Index>>(std::ref(index));
        treg.listen();

        if (partask::master()) {
            if (!storage().params.early_tc.length_bound) {
                INFO("Early tip clipper length bound set as (RL - K)");
            storage().params.early_tc.length_bound = cfg::get().ds.RL - gp.k();
            }

            size_t length_bound = *storage().params.early_tc.length_bound;
            size_t kpo_mers_removed = clip_tips(length_bound);
            INFO(kpo_mers_removed << " short edges ("  << (index.k() + 1) << "-mers) were removed by early tip clipper");
        }

        treg.stop_listening();
    }

    void load(graph_pack::GraphPack&,
              const std::filesystem::path &,
              const char*) override {
        VERIFY_MSG(false, "implement me");
    }

    void save(const graph_pack::GraphPack&,
              const std::filesystem::path &,
              const char*) const override {
        // VERIFY_MSG(false, "implement me");
    }
};

class EarlyATClipper : public ConstructionMPI::Phase {
public:
    EarlyATClipper()
            : ConstructionMPI::Phase("Early A/T remover", "early_at_remover") { }

    virtual ~EarlyATClipper() = default;

    bool distributed() const override { return true; }

    void run(graph_pack::GraphPack &, const char*) override {
        EarlyLowComplexityClipperProcessor at_processor(storage().ext_index, 0.8, 10, 200);
        at_processor.RemoveATEdges();
        at_processor.RemoveATTips();
    }

    void load(graph_pack::GraphPack&,
              const std::filesystem::path &,
              const char*) override {
        VERIFY_MSG(false, "implement me");
    }

    void save(const graph_pack::GraphPack&,
              const std::filesystem::path &,
              const char*) const override {
        // VERIFY_MSG(false, "implement me");
    }
};

class GraphCondenser : public ConstructionMPI::Phase {
public:
    GraphCondenser()
            : ConstructionMPI::Phase("Condensing graph (MPI)", "graph_condensing_mpi") { }

    virtual ~GraphCondenser() = default;

    bool distributed() const override { return true; }

    void run(graph_pack::GraphPack &gp, const char*) override {
        auto &index = gp.get_mutable<EdgeIndex<Graph>>();
        if (index.IsAttached())
            index.Detach();

        partask::TaskRegistry treg;
        using GraphT = std::decay_t<decltype(gp.get_mutable<Graph>())>;
        auto condense = treg.add<DeBruijnGraphExtentionConstructorTask<GraphT>>(std::ref(gp.get_mutable<Graph>()), std::ref(storage().ext_index));
        treg.listen();
        if (partask::master()) {
            condense(storage().params.keep_perfect_loops);
        }
        treg.stop_listening();
        INFO("Graph synced, edges " << gp.get<Graph>().e_size() << ", vertices " << gp.get<Graph>().size());
        VERIFY(partask::all_equal(gp.get<Graph>().e_size()));
        VERIFY(partask::all_equal(gp.get<Graph>().size()));
    }

    void load(graph_pack::GraphPack&,
              const std::filesystem::path &,
              const char*) override {
        VERIFY_MSG(false, "implement me");
    }

    void save(const graph_pack::GraphPack&,
              const std::filesystem::path &,
              const char*) const override {
        //FIXME why commented here and others
        // VERIFY_MSG(false, "implement me");
    }
};

template <typename Index>
class CollectKMerCoverageTask {
    using ReadStreams = io::ReadStreamList<io::SingleReadSeq>;

public:
    CollectKMerCoverageTask() = default;
    CollectKMerCoverageTask(std::istream &is) { deserialize(is); }
    std::ostream &serialize(std::ostream &os) const { return os; }
    std::istream &deserialize(std::istream &is) { return is; }

    auto make_splitter(size_t, ReadStreams &read_streams, Index &) {
        return partask::make_seq_along_generator(read_streams);
    }

    void process(std::istream &is, std::ostream &, ReadStreams &read_streams, Index &index) {
        auto chunks = partask::get_seq(is);
        if (!chunks.size())
            return;

        INFO("Selected streams: " << chunks);
        partask::execute_on_subset(read_streams, chunks,
                                   [&](ReadStreams& local_streams) {
                                       # pragma omp parallel for
                                       for (size_t i = 0; i < local_streams.size(); ++i)
                                           kmers::CoverageHashMapBuilder().FillCoverageFromStream(local_streams[i], index);
                                   });
    }

    void sync(ReadStreams & /*read_streams*/, Index &index) {
        auto &values = index.values();
        partask::allreduce(values.data(), values.size(), MPI_SUM);
    }
};


class PHMCoverageFiller : public ConstructionMPI::Phase {
public:
    PHMCoverageFiller()
            : ConstructionMPI::Phase("Filling coverage indices (PHM)", "coverage_filling_phm") {}
    virtual ~PHMCoverageFiller() = default;

    bool distributed() const override { return true; }

    void run(graph_pack::GraphPack &gp, const char *) override {
        if (!storage().kmers)
            storage().kmers.reset(new kmers::KMerDiskStorage<RtSeq>());
        partask::broadcast(*storage().kmers);

        storage().coverage_map.reset(new ConstructionStorage::CoverageMap(storage().kmers->k()));
        auto &coverage_map = *storage().coverage_map;
        auto &streams = storage().read_streams;

        unsigned nthreads = cfg::get().max_threads;
        kmers::PerfectHashMapBuilder().BuildIndex(coverage_map, *storage().kmers, nthreads);

        INFO("Collecting k-mer coverage information from reads, this takes a while.");
        {
            partask::TaskRegistry treg;
            auto fill_kmer_coverage = treg.add<CollectKMerCoverageTask<ConstructionStorage::CoverageMap>>(std::ref(streams), std::ref(coverage_map));
            treg.listen();
            if (partask::master())
                fill_kmer_coverage();
            treg.stop_listening();
        }

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

    void load(graph_pack::GraphPack&,
              const std::filesystem::path &,
              const char*) override {
        VERIFY_MSG(false, "implement me");
    }

    void save(const graph_pack::GraphPack&,
              const std::filesystem::path &,
              const char*) const override {
        // VERIFY_MSG(false, "implement me");
    }

};

} // namespace

ConstructionMPI::ConstructionMPI()
        : spades::MPICompositeStageDeferred<ConstructionStorage>("de Bruijn graph construction", "construction") {
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

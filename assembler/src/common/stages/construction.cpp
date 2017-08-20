//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "io/reads/vector_reader.hpp"
#include "io/dataset_support/dataset_readers.hpp"
#include "pipeline/graph_pack.hpp"
#include "io/dataset_support/read_converter.hpp"

#include "modules/graph_construction.hpp"
#include "assembly_graph/stats/picture_dump.hpp"
#include "construction.hpp"

namespace debruijn_graph {

struct ConstructionStorage {
    ConstructionStorage(unsigned k)
            : ext_index(k) {}

    utils::DeBruijnExtensionIndex<> ext_index;
    std::unique_ptr<utils::KMerDiskCounter<RtSeq>> counter;
    config::debruijn_config::construction params;
    io::ReadStreamList<io::SingleReadSeq> read_streams;
    io::SingleStreamPtr contigs_stream;
    utils::ReadStatistics read_stats;
    std::string workdir;
};

void ConstructionNew::init(debruijn_graph::conj_graph_pack &gp, const char *) {
    init_storage(unsigned(gp.g.k()));

    // Has to be separate stream for not counting it in coverage
    io::ReadStreamList<io::SingleRead> trusted_contigs;
    if (cfg::get().use_additional_contigs) {
        DEBUG("Contigs from previous K will be used: " << cfg::get().additional_contigs);
        trusted_contigs.push_back(io::EasyStream(cfg::get().additional_contigs, true));
    }

    bool trusted_contigs_exist = false;
    for (const auto& lib : cfg::get().ds.reads) {
        if (lib.type() != io::LibraryType::TrustedContigs)
            continue;

        for (const auto& read : lib.single_reads()) {
            trusted_contigs.push_back(io::EasyStream(read, true));
            trusted_contigs_exist = true;
        }
    }

    if (trusted_contigs_exist)
        INFO("Trusted contigs will be used in graph construction");
    storage().contigs_stream = MultifileWrap(trusted_contigs);

    // FIXME: indices here are awful
    auto& dataset = cfg::get_writable().ds;
    std::vector<size_t> libs_for_construction;
    for (size_t i = 0; i < dataset.reads.lib_count(); ++i)
        if (dataset.reads[i].is_graph_contructable())
            libs_for_construction.push_back(i);

    storage().read_streams = io::single_binary_readers_for_libs(dataset, libs_for_construction, true, true);
    storage().params = cfg::get().con;
    storage().workdir = fs::make_temp_dir(gp.workdir, "construction");
};

ConstructionNew::~ConstructionNew() {}

class KMerCounting : public ConstructionNew::Phase {
public:
    KMerCounting()
            : ConstructionNew::Phase("k+1-mer counting", "kpomer_counting") { }

    virtual ~KMerCounting() = default;

    void run(debruijn_graph::conj_graph_pack &, const char*) override {
        auto &read_streams = storage().read_streams;
        auto &contigs_stream = storage().contigs_stream;
        const auto &index = storage().ext_index;

        VERIFY_MSG(read_streams.size(), "No input streams specified");

        unsigned nthreads = (unsigned)read_streams.size();
        utils::DeBruijnReadKMerSplitter<io::SingleReadSeq,
                                        utils::StoringTypeFilter<decltype(storage().ext_index)::storing_type>>
                splitter(storage().workdir, index.k() + 1, 0xDEADBEEF,
                         read_streams,(contigs_stream == 0) ? 0 : &(*contigs_stream),
                         storage().params.read_buffer_size);
        storage().counter.reset(new utils::KMerDiskCounter<RtSeq>(storage().workdir, splitter));
        storage().counter->CountAll(nthreads, nthreads, /* merge */false);

        auto stats = splitter.stats();
        size_t rl = stats.max_read_length_;
        if (!cfg::get().ds.RL()) {
            INFO("Figured out: read length = " << rl);
            cfg::get_writable().ds.set_RL(rl);
            cfg::get_writable().ds.set_aRL((double) stats.bases_ / (double) stats.reads_);
        } else if (cfg::get().ds.RL() != rl)
            WARN("In datasets.info, wrong RL is specified: " << cfg::get().ds.RL() << ", not " << rl);
        storage().read_stats = stats;
    }

    void load(debruijn_graph::conj_graph_pack&,
              const std::string &,
              const char*) override {
        VERIFY_MSG(false, "implement me");
    }

    void save(const debruijn_graph::conj_graph_pack&,
              const std::string &,
              const char*) const override {
        VERIFY_MSG(false, "implement me");
    }
};


class ExtensionIndexBuilder : public ConstructionNew::Phase {
public:
    ExtensionIndexBuilder()
            : ConstructionNew::Phase("Extension index construction", "extension_index_construction") { }

    virtual ~ExtensionIndexBuilder() = default;

    void run(debruijn_graph::conj_graph_pack &, const char*) override {
        // FIXME: We just need files here, not the full counter. Implement refererence counting scheme!
        utils::DeBruijnExtensionIndexBuilder().BuildExtensionIndexFromKPOMers(storage().workdir,
                                                                              storage().ext_index,
                                                                              *storage().counter,
                                                                              unsigned(storage().read_streams.size()),
                                                                              storage().params.read_buffer_size);
    }

    void load(debruijn_graph::conj_graph_pack&,
              const std::string &,
              const char*) override {
        VERIFY_MSG(false, "implement me");
    }

    void save(const debruijn_graph::conj_graph_pack&,
              const std::string &,
              const char*) const override {
        VERIFY_MSG(false, "implement me");
    }
};


class EarlyTipClipper : public ConstructionNew::Phase {
public:
    EarlyTipClipper()
            : ConstructionNew::Phase("Early tip clipping", "early_tip_clipper") { }

    virtual ~EarlyTipClipper() = default;

    void run(debruijn_graph::conj_graph_pack &gp, const char*) override {
        size_t length_bound = storage().read_stats.max_read_length_ - gp.g.k();
        if (storage().params.early_tc.length_bound)
            length_bound = storage().params.early_tc.length_bound.get();
        AlternativeEarlyTipClipper(storage().ext_index, length_bound).ClipTips();
    }

    void load(debruijn_graph::conj_graph_pack&,
              const std::string &,
              const char*) override {
        VERIFY_MSG(false, "implement me");
    }

    void save(const debruijn_graph::conj_graph_pack&,
              const std::string &,
              const char*) const override {
        VERIFY_MSG(false, "implement me");
    }
};

class GraphCondenser : public ConstructionNew::Phase {
public:
    GraphCondenser()
            : ConstructionNew::Phase("Condensing graph", "graph_condensing") { }

    virtual ~GraphCondenser() = default;

    void run(debruijn_graph::conj_graph_pack &gp, const char*) override {
        VERIFY(!gp.index.IsAttached());
        DeBruijnGraphExtentionConstructor<Graph>(gp.g, storage().ext_index).ConstructGraph(storage().params.keep_perfect_loops);
    }

    void load(debruijn_graph::conj_graph_pack&,
              const std::string &,
              const char*) override {
        VERIFY_MSG(false, "implement me");
    }

    void save(const debruijn_graph::conj_graph_pack&,
              const std::string &,
              const char*) const override {
        VERIFY_MSG(false, "implement me");
    }
};

class EdgeIndexFiller : public ConstructionNew::Phase {
public:
    EdgeIndexFiller()
            : ConstructionNew::Phase("Edge index filling", "initial_edge_index_filling") { }

    virtual ~EdgeIndexFiller() = default;

    void run(debruijn_graph::conj_graph_pack &gp, const char*) override {
        gp.index.Refill();
        gp.index.Attach();
    }

    void load(debruijn_graph::conj_graph_pack&,
              const std::string &,
              const char*) override {
        VERIFY_MSG(false, "implement me");
    }

    void save(const debruijn_graph::conj_graph_pack&,
              const std::string &,
              const char*) const override {
        VERIFY_MSG(false, "implement me");
    }
};

class CoverageFiller : public ConstructionNew::Phase {
public:
    CoverageFiller()
            : ConstructionNew::Phase("Filling coverage indices", "coverage_filling") { }

    virtual ~CoverageFiller() = default;

    void run(debruijn_graph::conj_graph_pack &gp, const char*) override {
        typedef typename decltype(gp.index)::InnerIndex InnerIndex;
        typedef typename EdgeIndexHelper<InnerIndex>::CoverageAndGraphPositionFillingIndexBuilderT IndexBuilder;
        INFO("Filling coverage index");
        IndexBuilder().ParallelFillCoverage(gp.index.inner_index(), storage().read_streams);
        INFO("Filling coverage and flanking coverage from index");
        FillCoverageAndFlanking(gp.index.inner_index(), gp.g, gp.flanking_cov);
    }

    void load(debruijn_graph::conj_graph_pack&,
              const std::string &,
              const char*) override {
        VERIFY_MSG(false, "implement me");
    }

    void save(const debruijn_graph::conj_graph_pack&,
              const std::string &,
              const char*) const override {
        VERIFY_MSG(false, "implement me");
    }
};

class CQFCoverageFiller : public ConstructionNew::Phase {
public:
    CQFCoverageFiller()
            : ConstructionNew::Phase("Filling coverage indices (CQF)", "coverage_filling_cqf") {}
    virtual ~CQFCoverageFiller() = default;

    void run(debruijn_graph::conj_graph_pack &gp, const char *) override {
        INFO("Filling coverage and flanking coverage from index");
        FillCoverageAndFlanking(gp.index.inner_index(), gp.g, gp.flanking_cov);
    }

    void load(debruijn_graph::conj_graph_pack&,
              const std::string &,
              const char*) override {
        VERIFY_MSG(false, "implement me");
    }

    void save(const debruijn_graph::conj_graph_pack&,
              const std::string &,
              const char*) const override {
        VERIFY_MSG(false, "implement me");
    }

};


ConstructionNew::ConstructionNew()
        : spades::CompositeStageDeferred<ConstructionStorage>("de Bruijn graph construction", "construction") {
    add<KMerCounting>();
    add<ExtensionIndexBuilder>();
    if (cfg::get().con.early_tc.enable && !cfg::get().gap_closer_enable)
        add<EarlyTipClipper>();
    add<GraphCondenser>();
    add<EdgeIndexFiller>();
    add<CoverageFiller>();
}

template<class Read>
void construct_graph(io::ReadStreamList<Read>& streams,
                     conj_graph_pack& gp, io::SingleStreamPtr contigs_stream = io::SingleStreamPtr()) {
    config::debruijn_config::construction params = cfg::get().con;
    params.early_tc.enable &= !cfg::get().gap_closer_enable;

    std::string workdir = fs::make_temp_dir(gp.workdir, "construction");
    utils::ReadStatistics stats = ConstructGraphWithCoverage(params, workdir, streams, gp.g,
                                                             gp.index, gp.flanking_cov, contigs_stream);
    size_t rl = stats.max_read_length_;

    if (!cfg::get().ds.RL()) {
        INFO("Figured out: read length = " << rl);
        cfg::get_writable().ds.set_RL(rl);
        cfg::get_writable().ds.set_aRL((double) stats.bases_ / (double) stats.reads_);
    } else if (cfg::get().ds.RL() != rl)
        WARN("In datasets.info, wrong RL is specified: " << cfg::get().ds.RL() << ", not " << rl);
}

void Construction::run(conj_graph_pack &gp, const char*) {
    // Has to be separate stream for not counting it in coverage
    io::ReadStreamList<io::SingleRead> trusted_contigs;
    if (cfg::get().use_additional_contigs) {
        DEBUG("Contigs from previous K will be used: " << cfg::get().additional_contigs);
        trusted_contigs.push_back(io::EasyStream(cfg::get().additional_contigs, true));
    }

    bool trusted_contigs_exist = false;
    for (const auto& lib : cfg::get().ds.reads) {
        if (lib.type() != io::LibraryType::TrustedContigs)
            continue;

        for (const auto& read : lib.single_reads()) {
            trusted_contigs.push_back(io::EasyStream(read, true));
            trusted_contigs_exist = true;
        }
    }

    if (trusted_contigs_exist)
        INFO("Trusted contigs will be used in graph construction");
    auto contigs_stream = MultifileWrap(trusted_contigs);

    auto& dataset = cfg::get_writable().ds;
    std::vector<size_t> libs_for_construction;
    for (size_t i = 0; i < dataset.reads.lib_count(); ++i)
        if (dataset.reads[i].is_graph_contructable())
            libs_for_construction.push_back(i);

    auto streams = io::single_binary_readers_for_libs(dataset, libs_for_construction, true, true);
    construct_graph<io::SingleReadSeq>(streams, gp, contigs_stream);
}

} //namespace debruijn_graph

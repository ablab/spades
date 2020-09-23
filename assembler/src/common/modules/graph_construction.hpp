//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once


#include "assembly_graph/core/graph.hpp"

#include "assembly_graph/construction/debruijn_graph_constructor.hpp"
#include "assembly_graph/construction/early_simplification.hpp"
#include "assembly_graph/graph_support/coverage_filling.hpp"

#include "modules/alignment/edge_index.hpp"

// FIXME: layering violation
#include "pipeline/config_struct.hpp"
#include "utils/extension_index/kmer_extension_index_builder.hpp"
#include "utils/ph_map/coverage_hash_map_builder.hpp"
#include "utils/perf/perfcounter.hpp"

#include "io/reads/io_helper.hpp"

namespace debruijn_graph {

template<class ExtensionIndex>
void EarlyClipTips(const config::debruijn_config::construction& params, ExtensionIndex& ext) {
    if (!params.early_tc.enable)
        return;

    VERIFY(params.early_tc.length_bound);
    EarlyTipClipperProcessor(ext, *params.early_tc.length_bound).ClipTips();
}

using KMerFiles = kmers::KMerDiskStorage<RtSeq>;

template<class Read>
KMerFiles ConstructGraphUsingExtensionIndex(const config::debruijn_config::construction &params,
                                            fs::TmpDir workdir,
                                            io::ReadStreamList<Read>& streams, Graph& g) {
    unsigned k = unsigned(g.k());
    INFO("Constructing DeBruijn graph for k=" << k);

    TRACE("Filling indices");
    VERIFY_MSG(streams.size(), "No input streams specified");

    TRACE("... in parallel");
    utils::DeBruijnExtensionIndex<> ext(k);

    KMerFiles kmers = utils::DeBruijnExtensionIndexBuilder().BuildExtensionIndexFromStream(workdir, ext, streams, params.read_buffer_size);

    EarlyClipTips(params, ext);

    INFO("Condensing graph");
    DeBruijnGraphExtentionConstructor<Graph> g_c(g, ext);
    g_c.ConstructGraph(params.keep_perfect_loops);

    return kmers;
}

//FIXME these methods are tested, but not used!
template<class Streams>
KMerFiles ConstructGraph(const config::debruijn_config::construction &params,
                         fs::TmpDir workdir, Streams& streams, Graph& g) {
    return ConstructGraphUsingExtensionIndex(params, workdir, streams, g);
}

template<class Streams>
KMerFiles ConstructGraphWithIndex(const config::debruijn_config::construction &params,
                                  fs::TmpDir workdir, Streams& streams, Graph& g,
                                  EdgeIndex<Graph>& index) {
    VERIFY(!index.IsAttached());

    KMerFiles kmers = ConstructGraph(params, workdir, streams, g);

    INFO("Building index with from graph")
    //todo pass buffer size
    index.Refill();
    index.Attach();

    return kmers;
}

template<class Streams>
void ConstructGraphWithCoverage(const config::debruijn_config::construction &params,
                                fs::TmpDir workdir, Streams &streams, Graph &g,
                                EdgeIndex<Graph> &index, omnigraph::FlankingCoverage<Graph> &flanking_cov) {
    KMerFiles kmers = ConstructGraphWithIndex(params, workdir, streams, g, index);

    INFO("Filling coverage index");
    using CoverageMap = utils::PerfectHashMap<RtSeq, uint32_t, utils::slim_kmer_index_traits<RtSeq>, utils::DefaultStoring>;
    CoverageMap coverage_map(unsigned(g.k() + 1));

    utils::CoverageHashMapBuilder().BuildIndex(coverage_map, kmers, streams);

    INFO("Filling coverage and flanking coverage from PHM");
    FillCoverageAndFlankingFromPHM(coverage_map, g, flanking_cov);
}

}

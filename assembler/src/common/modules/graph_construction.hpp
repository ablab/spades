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
#include "assembly_graph/index/edge_index_builders.hpp"

// FIXME: layering violation
#include "pipeline/config_struct.hpp"
#include "utils/extension_index/kmer_extension_index_builder.hpp"
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

template<class Graph, class Read, class Index>
void ConstructGraphUsingExtensionIndex(const config::debruijn_config::construction &params,
                                       fs::TmpDir workdir,
                                       io::ReadStreamList<Read>& streams, Graph& g,
                                       Index& index) {
    unsigned k = unsigned(g.k());
    INFO("Constructing DeBruijn graph for k=" << k);

    TRACE("Filling indices");
    VERIFY_MSG(streams.size(), "No input streams specified");

    TRACE("... in parallel");
    utils::DeBruijnExtensionIndex<> ext(k);

    utils::DeBruijnExtensionIndexBuilder().BuildExtensionIndexFromStream(workdir, ext, streams, params.read_buffer_size);

    EarlyClipTips(params, ext);

    INFO("Condensing graph");
    VERIFY(!index.IsAttached());
    DeBruijnGraphExtentionConstructor<Graph> g_c(g, ext);
    g_c.ConstructGraph(params.keep_perfect_loops);

    INFO("Building index with from graph")
    //todo pass buffer size
    index.Refill();
    index.Attach();
}

//FIXME these methods are tested, but not used!
template<class Graph, class Index, class Streams>
void ConstructGraph(const config::debruijn_config::construction &params,
                    fs::TmpDir workdir, Streams& streams, Graph& g,
                    Index& index) {
    ConstructGraphUsingExtensionIndex(params, workdir, streams, g, index);
}

//FIXME these methods are tested, but not used!
template<class Graph, class Index, class Streams>
void ConstructGraphWithCoverage(const config::debruijn_config::construction &params,
                                fs::TmpDir workdir, Streams &streams, Graph &g,
                                Index &index, omnigraph::FlankingCoverage<Graph> &flanking_cov) {
    ConstructGraph(params, workdir, streams, g, index);

    typedef typename Index::InnerIndex InnerIndex;
    typedef typename EdgeIndexHelper<InnerIndex>::CoverageAndGraphPositionFillingIndexBuilderT IndexBuilder;
    INFO("Filling coverage index")
    IndexBuilder().ParallelFillCoverage(index.inner_index(), streams);
    INFO("Filling coverage and flanking coverage from index");
    FillCoverageAndFlanking(index.inner_index(), g, flanking_cov);
}

}

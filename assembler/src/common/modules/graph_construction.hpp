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

#include "modules/alignment/edge_index.hpp"

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

template<class Read>
void ConstructGraphUsingExtensionIndex(const config::debruijn_config::construction &params,
                                       fs::TmpDir workdir,
                                       io::ReadStreamList<Read>& streams, Graph& g) {
    unsigned k = unsigned(g.k());
    INFO("Constructing DeBruijn graph for k=" << k);

    TRACE("Filling indices");
    VERIFY_MSG(streams.size(), "No input streams specified");

    TRACE("... in parallel");
    utils::DeBruijnExtensionIndex<> ext(k);

    utils::DeBruijnExtensionIndexBuilder().BuildExtensionIndexFromStream(workdir, ext, streams, params.read_buffer_size);

    EarlyClipTips(params, ext);

    INFO("Condensing graph");
    DeBruijnGraphExtentionConstructor<Graph> g_c(g, ext);
    g_c.ConstructGraph(params.keep_perfect_loops);
}

//FIXME these methods are tested, but not used!
template<class Streams>
void ConstructGraph(const config::debruijn_config::construction &params,
                    fs::TmpDir workdir, Streams& streams, Graph& g) {
    ConstructGraphUsingExtensionIndex(params, workdir, streams, g);
}

template<class Streams>
void ConstructGraphWithIndex(const config::debruijn_config::construction &params,
                             fs::TmpDir workdir, Streams& streams, Graph& g,
                             EdgeIndex<Graph>& index) {
    VERIFY(!index.IsAttached());

    ConstructGraph(params, workdir, streams, g);

    INFO("Building index with from graph")
    //todo pass buffer size
    index.Refill();
    index.Attach();
}

//FIXME these methods are tested, but not used!
template<class Streams>
void ConstructGraphWithCoverage(const config::debruijn_config::construction &params,
                                fs::TmpDir workdir, Streams &streams, Graph &g,
                                EdgeIndex<Graph> &index, omnigraph::FlankingCoverage<Graph> &flanking_cov) {
    ConstructGraphWithIndex(params, workdir, streams, g, index);

    typedef typename EdgeIndex<Graph>::InnerIndex InnerIndex;
    typedef typename EdgeIndexHelper<InnerIndex>::CoverageAndGraphPositionFillingIndexBuilderT IndexBuilder;
    INFO("Filling coverage index")
    IndexBuilder().ParallelFillCoverage(index.inner_index(), streams);
    INFO("Filling coverage and flanking coverage from index");
    FillCoverageAndFlanking(index.inner_index(), g, flanking_cov);
}

}

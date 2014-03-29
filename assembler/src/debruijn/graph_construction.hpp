//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * graph_construction.hpp
 *
 *  Created on: Aug 12, 2011
 *      Author: sergey
 */
#pragma once

#include "openmp_wrapper.h"

#include "io/io_helper.hpp"
#include "omni/edges_position_handler.hpp"

#include "debruijn_graph_constructor.hpp"
#include "indices/edge_index_builders.hpp"
#include "debruijn_graph.hpp"
#include "graph_pack.hpp"
#include "utils.hpp"
#include "perfcounter.hpp"
#include "early_simplification.hpp"

#include "read_converter.hpp"
#include "detail_coverage.hpp"
#include "indices/storing_traits.hpp"

namespace debruijn_graph {

template<class StoringType>
struct CoverageCollector {
};

template<>
struct CoverageCollector<SimpleStoring> {
    template<class Info>
    static void CollectCoverage(Info edge_info) {
        edge_info.edge_id->IncCoverage(edge_info.count);
    }
};

template<>
struct CoverageCollector<InvertableStoring> {
    template<class Info>
    static void CollectCoverage(Info edge_info) {
        edge_info.edge_id->IncCoverage(edge_info.count);
        edge_info.edge_id->conjugate()->IncCoverage(edge_info.count);
    }
};


template<class Index>
void FillCoverageFromIndex(const Index &index) {
    for (auto I = index.value_cbegin(), E = index.value_cend();
            I != E; ++I) {
        const auto& edge_info = *I;
        VERIFY(edge_info.offset != -1u);
//      VERIFY(edge_info.edge_id.get() != NULL);
        if(edge_info.offset != -1u) {
            CoverageCollector<typename Index::storing_type>::CollectCoverage(edge_info);
        }
    }
    DEBUG("Coverage counted");
}

template<class Graph, class Readers, class Index>
size_t ConstructGraphUsingOldIndex(Readers& streams, Graph& g,
		Index& index, io::SingleStreamPtr contigs_stream = io::SingleStreamPtr()) {
	INFO("Constructing DeBruijn graph");

	TRACE("Filling indices");
	size_t rl = 0;
	VERIFY_MSG(streams.size(), "No input streams specified");

	TRACE("... in parallel");
	typedef typename Index::InnerIndexT InnerIndex;
	typedef typename EdgeIndexHelper<InnerIndex>::CoverageFillingEdgeIndexBuilderT IndexBuilder;
	InnerIndex& debruijn = index.inner_index();
	//fixme hack
	rl = IndexBuilder().BuildIndexFromStream(debruijn, streams, (contigs_stream == 0) ? 0 : &(*contigs_stream));

	VERIFY(g.k() + 1== debruijn.k());
	// FIXME: output_dir here is damn ugly!

	TRACE("Filled indices");

	INFO("Condensing graph");
	DeBruijnGraphConstructor<Graph, InnerIndex> g_c(g, debruijn);
	TRACE("Constructor ok");
	g_c.ConstructGraph(100, 10000, 1.2); // TODO: move magic constants to config
	INFO("Graph condensed");

	return rl;
}

inline debruijn_config::construction CreateDefaultConstructionConfig() {
    debruijn_config::construction config;
    config.con_mode = construction_mode::con_extention;
    debruijn_config::construction::early_tip_clipper early_tc;
    early_tc.enable = false;
    config.early_tc = early_tc;
    config.keep_perfect_loops = true;
    config.read_buffer_size = 0;
    return config;
}

template<class ExtensionIndex>
void EarlyClipTips(size_t k, const debruijn_config::construction params, size_t rl, ExtensionIndex& ext) {
    if (params.early_tc.enable) {
        size_t length_bound = rl - k;
        if (params.early_tc.length_bound)
            length_bound = params.early_tc.length_bound.get();
        AlternativeEarlyTipClipper(ext, length_bound).ClipTips();
    }
}

template<class Graph, class Read, class Index>
size_t ConstructGraphUsingExtentionIndex(const debruijn_config::construction params,
		io::ReadStreamList<Read>& streams, Graph& g,
		Index& index, io::SingleStreamPtr contigs_stream = io::SingleStreamPtr(), size_t read_buffer_size = 0) {

    size_t k = g.k();
	INFO("Constructing DeBruijn graph for k=" << k);

	TRACE("Filling indices");
	VERIFY_MSG(streams.size(), "No input streams specified");

	TRACE("... in parallel");
	// FIXME: output_dir here is damn ugly!
	typedef DeBruijnExtensionIndex<> ExtensionIndex;
	typedef typename ExtensionIndexHelper<ExtensionIndex>::DeBruijnExtensionIndexBuilderT ExtensionIndexBuilder;
	ExtensionIndex ext((unsigned) k, index.inner_index().workdir());

	//fixme hack
	size_t rl = ExtensionIndexBuilder().BuildExtensionIndexFromStream(ext, streams, (contigs_stream == 0) ? 0 : &(*contigs_stream), params.read_buffer_size);

	EarlyClipTips(k, params, rl, ext);

	INFO("Condensing graph");
	index.Detach();
	DeBruijnGraphExtentionConstructor<Graph> g_c(g, ext);
	g_c.ConstructGraph(100, 10000, 1.2, params.keep_perfect_loops);//TODO move these parameters to config
	index.Attach();

    typedef typename Index::InnerIndexT InnerIndex;
    typedef typename EdgeIndexHelper<InnerIndex>::CoverageAndGraphPositionFillingIndexBuilderT IndexBuilder;
	INFO("Building index with coverage from graph")
	IndexBuilder().BuildIndexFromGraph(index.inner_index(), g, read_buffer_size);
	IndexBuilder().ParallelFillCoverage(index.inner_index(), streams);
	return rl;
}

template<class Graph, class Index, class Streams>
size_t ConstructGraph(const debruijn_config::construction &params,
                      Streams& streams, Graph& g,
		 Index& index, io::SingleStreamPtr contigs_stream = io::SingleStreamPtr()) {
	if(params.con_mode == construction_mode::con_extention) {
		return ConstructGraphUsingExtentionIndex(params, streams, g, index, contigs_stream, params.read_buffer_size);
//	} else if(params.con_mode == construction_mode::con_old){
//		return ConstructGraphUsingOldIndex(k, streams, g, index, contigs_stream);
	} else {
		INFO("Invalid construction mode")
		VERIFY(false);
		return 0;
	}
}

template<class Graph, class Index, class Streams>
size_t ConstructGraphWithCoverage(const debruijn_config::construction &params,
                                  Streams& streams, Graph& g,
                                  Index& index, FlankingCoverage<Graph>& flanking_cov,
                                  io::SingleStreamPtr contigs_stream = io::SingleStreamPtr()) {
	size_t rl = ConstructGraph(params, streams, g, index, contigs_stream);

	INFO("Filling coverage and flanking coverage from index");
	FillCoverageAndFlanking(index.inner_index(), g, flanking_cov);
	return rl;
}

//template<class Graph, class Reader, class Index>
//size_t ConstructGraphWithCoverageFromStream(size_t k,
//        Reader& stream, Graph& g,
//        Index& index, SingleReadStream* contigs_stream = 0) {
//    io::ReadStreamVector<io::IReader<typename Reader::read_type>> streams(stream);
//    return ConstructGraph(k, streams, g, index, contigs_stream);
//}

}

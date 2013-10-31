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

#ifndef GRAPH_CONSTRUCTION_HPP_
#define GRAPH_CONSTRUCTION_HPP_

#include "openmp_wrapper.h"

#include "io/multifile_reader.hpp"
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

namespace debruijn_graph {

typedef io::IReader<io::SingleRead> SingleReadStream;
typedef io::IReader<io::PairedRead> PairedReadStream;
typedef io::MultifileReader<io::SingleRead> CompositeSingleReadStream;
typedef io::ConvertingReaderWrapper UnitedStream;

template<class Graph, class Readers, class Index>
size_t ConstructGraphUsingOldIndex(size_t k,
		Readers& streams, Graph& g,
		Index& index, SingleReadStream* contigs_stream = 0) {
	INFO("Constructing DeBruijn graph");

	TRACE("Filling indices");
	size_t rl = 0;
	VERIFY_MSG(streams.size(), "No input streams specified");

	TRACE("... in parallel");
	typedef typename Index::InnerIndexT InnerIndex;
	typedef typename EdgeIndexHelper<InnerIndex>::CoverageFillingEdgeIndexBuilderT IndexBuilder;
	InnerIndex& debruijn = index.inner_index();
	rl = IndexBuilder().BuildIndexFromStream(debruijn, streams, contigs_stream);

	VERIFY(k + 1== debruijn.k());
	// FIXME: output_dir here is damn ugly!

	TRACE("Filled indices");

	INFO("Condensing graph");
	DeBruijnGraphConstructor<Graph, InnerIndex> g_c(g, debruijn, k);
	TRACE("Constructor ok");
	g_c.ConstructGraph(100, 10000, 1.2); // TODO: move magic constants to config
	TRACE("Graph condensed");

	return rl;
}

template<class ExtensionIndex>
void EarlyClipTips(size_t k, const debruijn_config::construction params, size_t rl, ExtensionIndex& ext) {
    if (params.early_tc.enable) {
        size_t length_bound = rl - k;
        if (params.early_tc.length_bound)
            length_bound = params.early_tc.length_bound.get();
        EarlyTipClipper(ext, length_bound).ClipTips();
    }
}

template<class Graph, class Read, class Index>
size_t ConstructGraphUsingExtentionIndex(size_t k, const debruijn_config::construction params,
		io::ReadStreamVector<io::IReader<Read> >& streams, Graph& g,
		Index& index, SingleReadStream* contigs_stream = 0) {

	INFO("Constructing DeBruijn graph");

	TRACE("Filling indices");
	VERIFY_MSG(streams.size(), "No input streams specified");

	TRACE("... in parallel");
	// FIXME: output_dir here is damn ugly!
	typedef DeBruijnExtensionIndex<> ExtensionIndex;
	typedef typename ExtensionIndexHelper<ExtensionIndex>::DeBruijnExtensionIndexBuilderT ExtensionIndexBuilder;
	ExtensionIndex ext((unsigned) k, index.inner_index().workdir());
	size_t rl = ExtensionIndexBuilder().BuildExtensionIndexFromStream(ext, streams, contigs_stream);

	EarlyClipTips(k, params, rl, ext);

	INFO("Condensing graph");
	index.Detach();
	DeBruijnGraphExtentionConstructor<Graph> g_c(g, ext, k);
	g_c.ConstructGraph(100, 10000, 1.2, params.keep_perfect_loops);//TODO move these parameters to config
	index.Attach();

    typedef typename Index::InnerIndexT InnerIndex;
    typedef typename EdgeIndexHelper<InnerIndex>::CoverageAndGraphPositionFillingIndexBuilderT IndexBuilder;
	INFO("Building index with coverage from graph")
	IndexBuilder().BuildIndexFromGraph(index.inner_index(), g);
	IndexBuilder().ParallelFillCoverage(index.inner_index(), streams);
	return rl;
}

template<class Graph, class Index, class Streams>
size_t ConstructGraph(size_t k, const debruijn_config::construction &params,
                      Streams& streams, Graph& g,
		 Index& index, SingleReadStream* contigs_stream = 0) {
	if(params.con_mode == construction_mode::con_extention) {
		return ConstructGraphUsingExtentionIndex(k, params, streams, g, index, contigs_stream);
//	} else if(params.con_mode == construction_mode::con_old){
//		return ConstructGraphUsingOldIndex(k, streams, g, index, contigs_stream);
	} else {
		INFO("Invalid construction mode")
		VERIFY(false);
		return 0;
	}
}

template<class Graph, class Index, class Streams>
size_t ConstructGraphWithCoverage(size_t k, const debruijn_config::construction &params,
                                  Streams& streams, Graph& g,
                                  Index& index, FlankingCoverage<Graph>& flanking_cov,
                                  SingleReadStream* contigs_stream = 0) {
	size_t rl = ConstructGraph(k, params, streams, g, index, contigs_stream);

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

#endif /* GRAPH_CONSTRUCTION_HPP_ */

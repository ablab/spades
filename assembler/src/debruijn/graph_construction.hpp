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
#include "indices/debruijn_edge_index.hpp"
#include "debruijn_graph.hpp"
#include "graphio.hpp"
#include "graph_pack.hpp"
#include "utils.hpp"
#include "perfcounter.hpp"
#include "early_simplification.hpp"

#include "read_converter.hpp"

namespace debruijn_graph {

typedef io::IReader<io::SingleRead> SingleReadStream;
typedef io::IReader<io::PairedRead> PairedReadStream;
typedef io::MultifileReader<io::SingleRead> CompositeSingleReadStream;
typedef io::ConvertingReaderWrapper UnitedStream;

//fixme     NewExtendedSequenceMapper<Graph> mapper(g, index, kmer_mapper, k + 1);
template<class PairedRead, class Graph, class Mapper>
void FillPairedIndexWithReadCountMetric(const Graph &g,
                                        const Mapper& mapper,
                                        PairedInfoIndexT<Graph>& paired_info_index,
                                        io::ReadStreamVector<io::IReader<PairedRead> >& streams) {

	INFO("Counting paired info with read count weight");
	LatePairedIndexFiller<Graph, Mapper,
			io::IReader<PairedRead>> pif(g, mapper, streams,
			PairedReadCountWeight);

	pif.FillIndex(paired_info_index);
	DEBUG("Paired info with read count weight counted");
}

//fixme	NewExtendedSequenceMapper<Graph> mapper(g, index, kmer_mapper, k + 1);
template<class PairedRead, class Graph, class Mapper>
void FillPairedIndexWithProductMetric(const Graph &g,
                                      const Mapper& mapper,
                                      PairedInfoIndexT<Graph>& paired_info_index,
                                      io::ReadStreamVector<io::IReader<PairedRead> >& streams) {

	INFO("Counting paired info with product weight");

	LatePairedIndexFiller<Graph, Mapper,
			io::IReader<PairedRead> > pif(g, mapper, streams,
			KmerCountProductWeight);
	pif.FillIndex(paired_info_index);
	DEBUG("Paired info with product weight counted");
}

template<class Graph, class Index>
void FillEtalonPairedIndex(PairedInfoIndexT<Graph>& etalon_paired_index,
		const Graph &g, const Index& index,
		const KmerMapper<Graph>& kmer_mapper, size_t is, size_t rs,
		size_t delta, const Sequence& genome, size_t k)
{
	VERIFY_MSG(genome.size() > 0,
			"The genome seems not to be loaded, program will exit");
	INFO((string) (FormattedString("Counting etalon paired info for genome of length=%i, k=%i, is=%i, rs=%i, delta=%i")
	        << genome.size() << k << is << rs << delta));

	EtalonPairedInfoCounter<Graph, Index> etalon_paired_info_counter(g, index, kmer_mapper, is, rs, delta, k);
	etalon_paired_info_counter.FillEtalonPairedInfo(genome, etalon_paired_index);

	DEBUG("Etalon paired info counted");
}

template<class Graph, class Index>
void FillEtalonPairedIndex(PairedInfoIndexT<Graph>& etalon_paired_index,
		const Graph &g, const Index& index,
		const KmerMapper<Graph>& kmer_mapper, const Sequence& genome,
		size_t k) {

	FillEtalonPairedIndex(etalon_paired_index, g, index, kmer_mapper,
			(size_t)math::round(*cfg::get().ds.IS), *cfg::get().ds.RL, size_t(*cfg::get().ds.is_var),
			genome, k);
	//////////////////DEBUG
	//	SimpleSequenceMapper<k + 1, Graph> simple_mapper(g, index);
	//	Path<EdgeId> path = simple_mapper.MapSequence(genome);
	//	SequenceBuilder sequence_builder;
	//	sequence_builder.append(Seq<k>(g.EdgeNucls(path[0])));
	//	for (auto it = path.begin(); it != path.end(); ++it) {
	//		sequence_builder.append(g.EdgeNucls(*it).Subseq(k));
	//	}
	//	Sequence new_genome = sequence_builder.BuildSequence();
	//	NewEtalonPairedInfoCounter<k, Graph> new_etalon_paired_info_counter(g, index,
	//			insert_size, read_length, insert_size * 0.1);
	//	PairedInfoIndexT<Graph> new_paired_info_index(g);
	//	new_etalon_paired_info_counter.FillEtalonPairedInfo(new_genome, new_paired_info_index);
	//	CheckInfoEquality(etalon_paired_index, new_paired_info_index);
	//////////////////DEBUG
//	INFO("Etalon paired info counted");
}

template<class Index>
void FillCoverageFromIndex(Index& index) {
	const auto& inner_index = index.inner_index();

    for (auto I = inner_index.value_cbegin(), E = inner_index.value_cend();
            I != E; ++I) {
        const auto& edge_info = *I;
        VERIFY(edge_info.offset != -1u);
    }

	for (auto I = inner_index.value_cbegin(), E = inner_index.value_cend();
			I != E; ++I) {
		const auto& edge_info = *I;
		VERIFY(edge_info.edge_id.get() != NULL);
		edge_info.edge_id->IncCoverage(edge_info.count);
	}

	DEBUG("Coverage counted");
}

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
	typedef typename NewEdgeIndexHelper<Index>::CoverageFillingEdgeIndexBuilderT IndexBuilder;
	InnerIndex& debruijn = index.inner_index();
	rl = IndexBuilder().BuildIndexFromStream(debruijn, streams, contigs_stream);

	VERIFY(k + 1== debruijn.K());
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
	ExtensionIndex ext(k, index.inner_index().workdir());
	size_t rl = ExtensionIndexBuilder().BuildExtensionIndexFromStream(ext, streams, contigs_stream);

	EarlyClipTips(k, params, rl, ext);

	INFO("Condensing graph");
	index.Detach();
	DeBruijnGraphExtentionConstructor<Graph> g_c(g, ext, k);
	g_c.ConstructGraph(100, 10000, 1.2, params.keep_perfect_loops);//TODO move these parameters to config
	index.Attach();

    typedef typename Index::InnerIndexT InnerIndex;
    typedef typename NewEdgeIndexHelper</*Inner*/Index>::CoverageAndGraphPositionFillingIndexBuilderT IndexBuilder;
	INFO("Building index with coverage from graph")
	IndexBuilder().BuildIndexWithCoverageFromGraph(g, index.inner_index(), streams, contigs_stream);
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
                                  Index& index, SingleReadStream* contigs_stream = 0) {
	size_t rl = ConstructGraph(k, params, streams, g, index, contigs_stream);

	INFO("Filling coverage from index")
	FillCoverageFromIndex(index);

	return rl;
}

template<class Graph, class Reader, class Index>
size_t ConstructGraphFromStream(size_t k, const debruijn_config::construction params,
        Reader& stream, Graph& g,
        Index& index, SingleReadStream* contigs_stream = 0) {
    io::ReadStreamVector<io::IReader<typename Reader::read_type>> streams(stream);
    return ConstructGraph(k, params, streams, g, index, contigs_stream);
}

template<class Graph, class Reader, class Index>
size_t ConstructGraphWithCoverageFromStream(size_t k,
        Reader& stream, Graph& g,
        Index& index, SingleReadStream* contigs_stream = 0) {
    io::ReadStreamVector<io::IReader<typename Reader::read_type>> streams(stream);
    return ConstructGraph(k, streams, g, index, contigs_stream);
}

}

#endif /* GRAPH_CONSTRUCTION_HPP_ */

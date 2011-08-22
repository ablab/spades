/*
 * graph_construction.hpp
 *
 *  Created on: Aug 12, 2011
 *      Author: sergey
 */

#ifndef GRAPH_CONSTRUCTION_HPP_
#define GRAPH_CONSTRUCTION_HPP_

#include "debruijn_graph_constructor.hpp"
#include "omni/edges_position_handler.hpp"
#include "new_debruijn.hpp"
#include "omni/paired_info.hpp"
#include "graphio.hpp"
#include "utils.hpp"

namespace debruijn_graph {
typedef io::IReader<io::SingleRead> SingleReadStream;
typedef io::IReader<io::PairedRead> PairedReadStream;
typedef io::ConvertingReaderWrapper UnitedStream;

template<size_t k, class ReadStream>
void FillPairedIndex(const Graph &g, const EdgeIndex<k + 1, Graph>& index
		, PairedInfoIndex<Graph>& paired_info_index,
		ReadStream& stream) {
	typedef SimpleSequenceMapper<k + 1, Graph> SequenceMapper;
	INFO("-----------------------------------------");
	stream.reset();
	INFO("Counting paired info");
	SequenceMapper mapper(g, index);
	PairedIndexFiller<k + 1, Graph, SequenceMapper, ReadStream> pif(g, mapper,
			stream);
	pif.FillIndex(paired_info_index);
	INFO("Paired info counted");
}

template<size_t k, class ReadStream>
void FillPairedIndexWithReadCountMetric(const Graph &g,
		const EdgeIndex<k + 1, Graph>& index
		, const KmerMapper<k + 1, Graph>& kmer_mapper
		, PairedInfoIndex<Graph>& paired_info_index , ReadStream& stream) {
	INFO("-----------------------------------------");
	stream.reset();
	INFO("Counting paired info with read count weight");
	ExtendedSequenceMapper<k + 1, Graph> mapper(g, index, kmer_mapper);
	ReadCountPairedIndexFiller<k + 1, Graph, ReadStream> pif(g, mapper, stream);
	pif.FillIndex(paired_info_index);
	INFO("Paired info with read count weight counted");
}

template<size_t k>
void FillEtalonPairedIndex(const Graph &g,
		PairedInfoIndex<Graph>& etalon_paired_index,
		const EdgeIndex<k + 1, Graph>& index,
		const Sequence& genome) {
	INFO("-----------------------------------------");
	INFO("Counting etalon paired info");

	EtalonPairedInfoCounter<k, Graph> etalon_paired_info_counter(g, index,
			cfg::get().ds.IS, cfg::get().ds.RL, cfg::get().ds.IS * 0.1);
	etalon_paired_info_counter.FillEtalonPairedInfo(genome,
			etalon_paired_index);
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
//	PairedInfoIndex<Graph> new_paired_info_index(g);
//	new_etalon_paired_info_counter.FillEtalonPairedInfo(new_genome, new_paired_info_index);
//	CheckInfoEquality(etalon_paired_index, new_paired_info_index);
	//////////////////DEBUG
	INFO("Etalon paired info counted");
}

template<size_t k>
void FillEtalonPairedIndex(const Graph &g,
		PairedInfoIndex<Graph>& etalon_paired_index,
		const EdgeIndex<k + 1, Graph>& index,
		size_t insert_size, size_t read_size,
		const Sequence& genome) {
	INFO("-----------------------------------------");
	INFO("Counting etalon paired info");

	EtalonPairedInfoCounter<k, Graph> etalon_paired_info_counter(g, index,
			insert_size, read_size, insert_size * 0.1);
	etalon_paired_info_counter.FillEtalonPairedInfo(genome,
			etalon_paired_index);

	INFO("Paired info counted");
}

template<size_t k, class ReadStream>
void FillCoverage(Graph& g, ReadStream& stream,
		EdgeIndex<k + 1, Graph>& index) {
	typedef SimpleSequenceMapper<k + 1, Graph> SequenceMapper;
	INFO("-----------------------------------------");
	stream.reset();
	INFO("Counting coverage");
	SequenceMapper read_threader(g, index);
	g.coverage_index().FillIndex<SequenceMapper>(stream,
			read_threader);
	INFO("Coverage counted");
}

template<size_t k, class ReadStream>
void ConstructGraph(Graph& g, EdgeIndex<k + 1, Graph>& index,
ReadStream& stream) {
	typedef SeqMap<k + 1, typename Graph::EdgeId> DeBruijn;
	INFO("-----------------------------------------");
	INFO("Constructing DeBruijn graph");
	DeBruijn& debruijn = index.inner_index();
	INFO("Filling DeBruijn graph");
	debruijn.Fill(stream);
	INFO("DeBruijn graph constructed");

	INFO("Condensing graph");
	DeBruijnGraphConstructor<k, Graph> g_c(debruijn);
	g_c.ConstructGraph(g, index);
	INFO("Graph condensed");
}

template<size_t k>
void ConstructGraphWithCoverage(Graph& g, EdgeIndex<k + 1, Graph>& index
		, SingleReadStream& stream) {
	ConstructGraph<k>(g, index, stream);
	FillCoverage<k>(g, stream, index);
}

template<size_t k>
void ConstructGraphWithPairedInfo(Graph& g, EdgeIndex<k + 1, Graph>& index,
		PairedInfoIndex<Graph>& paired_index, PairedReadStream& stream) {
		UnitedStream united_stream(stream);
		ConstructGraphWithCoverage<k>(g, index,
				united_stream);
		FillPairedIndex<k>(g, index, paired_index, stream);
}

template<size_t k>
void ConstructGraphWithEtalonPairedInfo(Graph& g, EdgeIndex<k + 1, Graph>& index,
		PairedInfoIndex<Graph>& paired_index,
		PairedReadStream& stream, const Sequence& genome) {
	UnitedStream united_stream(stream);
	ConstructGraphWithCoverage<k>(g, index,
			united_stream);
	FillEtalonPairedIndex<k>(g, paired_index, index,
			genome);
}

}

#endif /* GRAPH_CONSTRUCTION_HPP_ */

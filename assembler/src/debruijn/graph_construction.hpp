/*
 * graph_construction.hpp
 *
 *  Created on: Aug 12, 2011
 *      Author: sergey
 */

#ifndef GRAPH_CONSTRUCTION_HPP_
#define GRAPH_CONSTRUCTION_HPP_

#include "io/multifile_reader.hpp"
#include "debruijn_graph_constructor.hpp"
#include "omni/edges_position_handler.hpp"
#include "new_debruijn.hpp"
#include "omni/paired_info.hpp"
#include "graphio.hpp"
#include "graph_pack.hpp"
#include "utils.hpp"

namespace debruijn_graph {
typedef io::IReader<io::SingleRead> SingleReadStream;
typedef io::IReader<io::PairedRead> PairedReadStream;
typedef io::MultifileReader<io::SingleRead> CompositeSingleReadStream;
typedef io::ConvertingReaderWrapper UnitedStream;

template<size_t k, class ReadStream>

void FillPairedIndexWithReadCountMetric(const Graph &g,
		const EdgeIndex<k + 1, Graph>& index
		, const KmerMapper<k + 1, Graph>& kmer_mapper
		, PairedInfoIndex<Graph>& paired_info_index , ReadStream& stream) {
	INFO("-----------------------------------------");
	stream.reset();
	INFO("Counting paired info with read count weight");
	ExtendedSequenceMapper<k + 1, Graph> mapper(g, index, kmer_mapper);
	LatePairedIndexFiller<k + 1, Graph, ReadStream> pif(g, mapper, stream, PairedReadCountWeight);
	pif.FillIndex(paired_info_index);
	INFO("Paired info with read count weight counted");
}

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

template<size_t k>
void FillEtalonPairedIndex(const Graph &g,
		PairedInfoIndex<Graph>& etalon_paired_index,
		const EdgeIndex<k + 1, Graph>& index, const Sequence& genome) {
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
		size_t is, size_t rs,
		const Sequence& genome) {
	INFO("-----------------------------------------");
	INFO("Counting etalon paired info");

	EtalonPairedInfoCounter<k, Graph> etalon_paired_info_counter(g, index,
			is, rs, is * 0.1);
	etalon_paired_info_counter.FillEtalonPairedInfo(genome,
			etalon_paired_index);

	INFO("Etalon paired info counted");
}

template<size_t k>
void FillCoverage(Graph& g, SingleReadStream& stream,
		EdgeIndex<k + 1, Graph>& index) {
	typedef SimpleSequenceMapper<k + 1, Graph> SequenceMapper;
	INFO("-----------------------------------------");
	stream.reset();
	INFO("Counting coverage");
	SequenceMapper read_threader(g, index);
	g.coverage_index().FillIndex<SequenceMapper>(stream, read_threader);
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
void ConstructGraphWithCoverage(Graph& g, EdgeIndex<k + 1, Graph>& index,
SingleReadStream& stream, SingleReadStream* contigs_stream = 0) {
	vector<SingleReadStream*> streams;
	streams.push_back(&stream);
	if (contigs_stream) {
		INFO("Additional contigs stream added for construction");
		streams.push_back(contigs_stream);
	}
	CompositeSingleReadStream composite_stream(streams);
	ConstructGraph<k>(g, index, composite_stream);
	//It is not a bug!!! Don't use composite_stream here!!!
	FillCoverage<k>(g, stream, index);
}

template<size_t k>
void ConstructGraphWithPairedInfo(conj_graph_pack& gp,
		PairedInfoIndex<Graph>& paired_index, PairedReadStream& stream,
		SingleReadStream* contigs_stream = 0) {
	UnitedStream united_stream(stream);
	ConstructGraphWithCoverage<k>(gp.g, gp.index, united_stream,
			contigs_stream);
	if (cfg::get().etalon_info_mode)
		FillEtalonPairedIndex<k>(gp.g, paired_index, gp.index, gp.genome);
	else
		FillPairedIndex<k>(gp.g, gp.index, paired_index, stream);
}

}

#endif /* GRAPH_CONSTRUCTION_HPP_ */

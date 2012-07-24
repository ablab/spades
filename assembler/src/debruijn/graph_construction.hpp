//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
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
#include "debruijn_graph_constructor.hpp"
#include "omni/edges_position_handler.hpp"
#include "new_debruijn.hpp"
#include "omni/paired_info.hpp"
#include "graphio.hpp"
#include "graph_pack.hpp"
#include "utils.hpp"
#include "parallel_seq_map.hpp"
#include "perfcounter.hpp"
#include "omni/parallel_unordered_map.hpp"

#include "read_converter.hpp"

namespace debruijn_graph {

typedef io::IReader<io::SingleRead> SingleReadStream;
typedef io::IReader<io::PairedRead> PairedReadStream;
typedef io::MultifileReader<io::SingleRead> CompositeSingleReadStream;
typedef io::ConvertingReaderWrapper UnitedStream;

template<class PairedRead>
void FillPairedIndexWithReadCountMetric(const Graph &g,
		const IdTrackHandler<Graph>& int_ids, const EdgeIndex<Graph>& index,
		const KmerMapper<Graph>& kmer_mapper,
		PairedInfoIndex<Graph>& paired_info_index,
		io::ReadStreamVector< io::IReader<PairedRead> >& streams, size_t k) {

	INFO("Counting paired info with read count weight");
	NewExtendedSequenceMapper<Graph> mapper(g, index, kmer_mapper, k + 1);
	LatePairedIndexFiller<Graph, NewExtendedSequenceMapper<Graph>,
			io::IReader<PairedRead> > pif(g, mapper, streams,
			PairedReadCountWeight);

//	ExtendedSequenceMapper<k + 1, Graph> mapper(g, int_ids, index, kmer_mapper);
//	LatePairedIndexFiller<k + 1, Graph, ExtendedSequenceMapper<k + 1, Graph>, ReadStream> pif(g, mapper, stream, PairedReadCountWeight);

	pif.FillIndex(paired_info_index);
	DEBUG("Paired info with read count weight counted");
}

template<class PairedRead>
void FillPairedIndexWithProductMetric(const Graph &g,
		const EdgeIndex<Graph>& index, const KmerMapper<Graph>& kmer_mapper,
		PairedInfoIndex<Graph>& paired_info_index,
		io::ReadStreamVector< io::IReader<PairedRead> >& streams, size_t k) {

	INFO("Counting paired info with product weight");

	//	ExtendedSequenceMapper<k + 1, Graph> mapper(g, int_ids, index, kmer_mapper);
	//	LatePairedIndexFiller<k + 1, Graph, ExtendedSequenceMapper<k + 1, Graph>, ReadStream> pif(g, mapper, stream, PairedReadCountWeight);

	NewExtendedSequenceMapper<Graph> mapper(g, index, kmer_mapper, k + 1);
	LatePairedIndexFiller<Graph, NewExtendedSequenceMapper<Graph>,
			io::IReader<PairedRead> > pif(g, mapper, streams,
			KmerCountProductWeight);
	pif.FillIndex(paired_info_index);
	DEBUG("Paired info with product weight counted");
}

void FillPairedIndex(const Graph &g, const EdgeIndex<Graph>& index,
		PairedInfoIndex<Graph>& paired_info_index,
		io::IReader<io::PairedRead>& stream, size_t k) {
	typedef SimpleSequenceMapper<Graph> SequenceMapper;
	stream.reset();
	INFO("Counting paired info");
	SequenceMapper mapper(g, index, k + 1);
	PairedIndexFiller<Graph, SequenceMapper, io::IReader<io::PairedRead> > pif(
			g, mapper, stream);
	pif.FillIndex(paired_info_index);
	DEBUG("Paired info counted");
}

void FillEtalonPairedIndex(PairedInfoIndex<Graph>& etalon_paired_index,
		const Graph &g, const EdgeIndex<Graph>& index,
		const KmerMapper<Graph>& kmer_mapper, size_t is, size_t rs,
		size_t delta, const Sequence& genome, size_t k) {

	INFO(
			(string) (FormattedString(
					"Counting etalon paired info for genome of length=%i, k=%i, is=%i, rs=%i, delta=%i")
					<< genome.size() << k << is << rs << delta));

	EtalonPairedInfoCounter<Graph> etalon_paired_info_counter(g, index,
			kmer_mapper, is, rs, delta, k);
	etalon_paired_info_counter.FillEtalonPairedInfo(genome,
			etalon_paired_index);

	DEBUG("Etalon paired info counted");
}

void FillEtalonPairedIndex(PairedInfoIndex<Graph>& etalon_paired_index,
		const Graph &g, const EdgeIndex<Graph>& index,
		const KmerMapper<Graph>& kmer_mapper, const Sequence& genome,
		size_t k) {

	FillEtalonPairedIndex(etalon_paired_index, g, index, kmer_mapper,
			*cfg::get().ds.IS, *cfg::get().ds.RL, size_t(*cfg::get().ds.is_var),
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
	//	PairedInfoIndex<Graph> new_paired_info_index(g);
	//	new_etalon_paired_info_counter.FillEtalonPairedInfo(new_genome, new_paired_info_index);
	//	CheckInfoEquality(etalon_paired_index, new_paired_info_index);
	//////////////////DEBUG
//	INFO("Etalon paired info counted");
}

template<class Read>
void FillCoverage(io::ReadStreamVector< io::IReader<Read> >& streams, Graph& g,
		EdgeIndex<Graph>& index, size_t k) {

	typedef SimpleSequenceMapper<Graph> SequenceMapper;

	INFO("Counting coverage");
	SequenceMapper read_threader(g, index, k + 1);

	if (streams.size() > 1) {
		g.coverage_index().FillFastParallelIndex<SequenceMapper, Read>(streams,
				read_threader);
	} else if (streams.size() == 1) {
		g.coverage_index().FillIndex<SequenceMapper, Read>(streams.back(),
				read_threader);
	}

	DEBUG("Coverage counted");
}

template<class Graph, class Read>
size_t FillUsusalIndex(io::IReader<Read>& stream,
		SeqMap<typename Graph::EdgeId>& debruijn, size_t k) {
	INFO("Processing reads (takes a while)");

	size_t counter = 0;
	size_t rl = 0;
	stream.reset();

	Read r;
	while (!stream.eof()) {
		stream >> r;
		Sequence s = r.sequence();
		debruijn.CountSequence(s);
		rl = max(rl, s.size());
		VERBOSE_POWER(++counter, " reads processed");
	}

	INFO("DeBruijn graph constructed, " << counter << " reads used");

	return rl;
}

template<class Graph, class Read>
size_t FillIterativeParallelIndex(io::ReadStreamVector< io::IReader<Read> >& streams,
		SeqMap<typename Graph::EdgeId>& debruijn, size_t k) {
	size_t k_plus_1 = k + 1;

	size_t nthreads = streams.size();
	for (size_t i = 0; i < nthreads; ++i) {
		streams[i].reset();
	}

	vector<typename ParallelSeqVector::destination_container_t> temp_sets;
	for (size_t i = 0; i < nthreads; ++i) {
		temp_sets.push_back(runtime_k::GetSet(k_plus_1));
	}

	perf_counter pc;

	INFO("Processing reads (takes a while)");
	std::vector<size_t> rls(nthreads, 0);
	size_t counter = 0;

	{
		size_t cell_size = cfg::get().buffer_size
				/ (nthreads * nthreads
						* ((k_plus_1 / (4 * sizeof(seq_element_type)) + 1)
								* sizeof(seq_element_type)));

		ParallelSeqVector par_debruijn(k_plus_1, nthreads, cell_size);
		while (!streams.eof()) {

#pragma omp parallel num_threads(nthreads)
			{
#pragma omp for reduction(+ : counter)
				for (size_t i = 0; i < nthreads; ++i) {

					Read r;
					io::IReader<Read>& stream = streams[i];

					while (!par_debruijn.IsFull(i) && !stream.eof()) {
						stream >> r;
						if (r.size() > rls[i]) {
							rls[i] = r.size();
						}
						++counter;

						par_debruijn.CountSequence(r.sequence(), i);
					}
				}

#pragma omp barrier

				//Merge maps
#pragma omp for
				for (size_t i = 0; i < nthreads; ++i) {
					par_debruijn.Dump(temp_sets[i], i);
				}
			}

		}
	}

	size_t total_kmers = 0;
	for (size_t i = 0; i < nthreads; ++i) {
		total_kmers += temp_sets[i].size();
	}

	//Merging into final map
	INFO("Merging final maps");

	debruijn.nodes().rehash(total_kmers);
	for (size_t i = 0; i < nthreads; ++i) {
		debruijn.transfer(temp_sets[i]);
		temp_sets[i].clear();
	}

	INFO("Elapsed time: " << pc.time_ms());

	INFO("DeBruijn graph constructed, reads used: " << counter);

	size_t rl = 0;
	for (size_t i = 0; i < nthreads; ++i) {
		if (rl < rls[i]) {
			rl = rls[i];
		}
	}

	return rl;
}

template<class Graph, class Read>
size_t ConstructGraph(size_t k, io::ReadStreamVector< io::IReader<Read> >& streams,
		Graph& g, EdgeIndex<Graph>& index,
		SingleReadStream* contigs_stream = 0) {

	typedef SeqMap<typename Graph::EdgeId> DeBruijn;

	INFO("Constructing DeBruijn graph");

	DeBruijn& debruijn = index.inner_index();

	TRACE("Filling indices");
	size_t rl = 0;
	if (streams.size() > 1) {
		TRACE("... in parallel");
		rl = FillIterativeParallelIndex<Graph, Read>(streams, debruijn, k);
	} else if (streams.size() == 1) {
		rl = FillUsusalIndex<Graph, Read>(streams.back(), debruijn, k);
	} else {
		VERIFY_MSG(false, "No input streams specified");
	}
	TRACE("Filled indices");

	io::SingleRead r;
	if (contigs_stream) {
		INFO("Adding contigs from previous K");
		while (!contigs_stream->eof()) {
			*contigs_stream >> r;
			Sequence s = r.sequence();
			debruijn.CountSequence(s);
		}
		TRACE("Added contigs from previous K");
	}

	INFO("Condensing graph");
	DeBruijnGraphConstructor<Graph> g_c(debruijn, k);
	g_c.ConstructGraph(g, index);
	TRACE("Graph condensed");

	return rl;
}

template<class Read>
size_t ConstructGraphWithCoverage(size_t k,
        io::ReadStreamVector< io::IReader<Read> >& streams, Graph& g,
		EdgeIndex<Graph>& index, SingleReadStream* contigs_stream = 0) {
	size_t rl = ConstructGraph(k, streams, g, index, contigs_stream);
	FillCoverage<Read>(streams, g, index, k);
	return rl;
}

size_t ConstructGraphWithPairedInfo(size_t k, Graph& g, EdgeIndex<Graph>& index,
		PairedInfoIndex<Graph>& paired_index, PairedReadStream& paired_stream,
		SingleReadStream* single_stream = 0, SingleReadStream* contigs_stream =
				0) {
	DEBUG("Constructing DeBruijn graph with paired info");

	UnitedStream united_stream(paired_stream);

	vector<SingleReadStream*> streams;
	streams.push_back(&united_stream);
	if (single_stream) {
		streams.push_back(single_stream);
	}
	CompositeSingleReadStream reads_stream(streams);
	io::ReadStreamVector<SingleReadStream> strs(&reads_stream);

	size_t rl = ConstructGraphWithCoverage<io::SingleRead>(k, strs, g, index,
			contigs_stream);
	FillPairedIndex(g, index, paired_index, paired_stream, k);
	return rl;
}

}

#endif /* GRAPH_CONSTRUCTION_HPP_ */

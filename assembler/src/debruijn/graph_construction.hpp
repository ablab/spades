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

#include <omp.h>
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

template<size_t k>
void FillPairedIndexWithReadCountMetric(const Graph &g, const IdTrackHandler<Graph>& int_ids,
		const EdgeIndex<k + 1, Graph>& index
		, const KmerMapper<k + 1, Graph>& kmer_mapper
		, PairedInfoIndex<Graph>& paired_info_index , io::IReader<io::PairedRead>& stream) {
	stream.reset();
	INFO("Counting paired info with read count weight");
	NewExtendedSequenceMapper<k + 1, Graph> mapper(g, index, kmer_mapper);
	LatePairedIndexFiller<k + 1, Graph, NewExtendedSequenceMapper<k + 1, Graph>> pif(g, mapper, stream, PairedReadCountWeight);
//	ExtendedSequenceMapper<k + 1, Graph> mapper(g, int_ids, index, kmer_mapper);
//	LatePairedIndexFiller<k + 1, Graph, ExtendedSequenceMapper<k + 1, Graph>, ReadStream> pif(g, mapper, stream, PairedReadCountWeight);
	pif.FillIndex(paired_info_index);
	DEBUG("Paired info with read count weight counted");
}

template<size_t k>
void FillPairedIndexWithProductMetric(const Graph &g
		, const EdgeIndex<k + 1, Graph>& index
		, const KmerMapper<k + 1, Graph>& kmer_mapper
		, PairedInfoIndex<Graph>& paired_info_index , io::IReader<io::PairedRead>& stream) {
	stream.reset();
	INFO("Counting paired info with product weight");
	//	ExtendedSequenceMapper<k + 1, Graph> mapper(g, int_ids, index, kmer_mapper);
	//	LatePairedIndexFiller<k + 1, Graph, ExtendedSequenceMapper<k + 1, Graph>, ReadStream> pif(g, mapper, stream, PairedReadCountWeight);
	NewExtendedSequenceMapper<k + 1, Graph> mapper(g, index, kmer_mapper);
	LatePairedIndexFiller<k + 1, Graph, NewExtendedSequenceMapper<k + 1, Graph>> pif(g, mapper, stream, KmerCountProductWeight);
	pif.FillIndex(paired_info_index);
	DEBUG("Paired info with product weight counted");
}

template<size_t k>
void FillPairedIndex(const Graph &g, const EdgeIndex<k + 1, Graph>& index
		, PairedInfoIndex<Graph>& paired_info_index,
		io::IReader<io::PairedRead>& stream) {
	typedef SimpleSequenceMapper<k + 1, Graph> SequenceMapper;
	stream.reset();
	INFO("Counting paired info");
	SequenceMapper mapper(g, index);
	PairedIndexFiller<k + 1, Graph, SequenceMapper> pif(g, mapper,
			stream);
	pif.FillIndex(paired_info_index);
	DEBUG("Paired info counted");
}

template<size_t k>
void FillEtalonPairedIndex(PairedInfoIndex<Graph>& etalon_paired_index,
		const Graph &g,
		const EdgeIndex<k + 1, Graph>& index,
		const KmerMapper<k + 1, Graph>& kmer_mapper,
		size_t is, size_t rs,
        size_t delta,
		const Sequence& genome){
	INFO((string)(FormattedString("Counting etalon paired info for genome of length=%i, k=%i, is=%i, rs=%i, delta=%i")
			<< genome.size() << k << is << rs << delta));
	EtalonPairedInfoCounter<k, Graph> etalon_paired_info_counter(g, index, kmer_mapper,
			is, rs, delta);
	etalon_paired_info_counter.FillEtalonPairedInfo(genome,
			etalon_paired_index);

	DEBUG("Etalon paired info counted");
}

template<size_t k>
void FillEtalonPairedIndex(PairedInfoIndex<Graph>& etalon_paired_index,
		const Graph &g,
		const EdgeIndex<k + 1, Graph>& index,
		const KmerMapper<k+1, Graph>& kmer_mapper,
		const Sequence& genome) {
	FillEtalonPairedIndex<k>(etalon_paired_index, g, index, kmer_mapper, *cfg::get().ds.IS, *cfg::get().ds.RL, size_t(*cfg::get().ds.is_var), genome);
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

template<size_t k, class Read>
void FillCoverage(std::vector<io::IReader<Read>* >& streams, Graph& g, EdgeIndex<k + 1, Graph>& index) {

	typedef SimpleSequenceMapper<k + 1, Graph> SequenceMapper;

	INFO("Counting coverage");
	SequenceMapper read_threader(g, index);

	if (streams.size() > 1) {
        g.coverage_index().FillParallelIndex<SequenceMapper, Read>(streams, read_threader);
	}
	else if (streams.size() == 1) {
	    g.coverage_index().FillIndex<SequenceMapper, Read>(*streams.back(), read_threader);
	}

	DEBUG("Coverage counted");
}

template<size_t k, class Graph, class Read>
size_t FillUsusalIndex(io::IReader<Read>& stream, SeqMap<k + 1, typename Graph::EdgeId>& debruijn) {

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


template<size_t k, class Graph, class Container>
void InsertIntoDebruijn(SeqMap<k + 1, typename Graph::EdgeId>& debruijn, const Container& container) {
    for (auto iter = container.begin(); iter != container.end(); ++iter) {
        debruijn.addEdge(*iter);
    }
}


template<size_t k, class Graph, class Read>
size_t FillIterativeParallelIndex(std::vector<io::IReader<Read>* >& streams, SeqMap<k + 1, typename Graph::EdgeId>& debruijn) {

    typedef ParallelSeqVector<k + 1> ParallelDeBruijn;
    typedef Seq<k + 1> Kmer;

    size_t nthreads = streams.size();
    for (size_t i = 0; i < nthreads; ++i) {
        streams[i]->reset();
    }

    vector<typename ParallelDeBruijn::destination_container_t> temp_maps(nthreads);

    perf_counter pc;

    INFO("Processing reads (takes a while)");
    std::vector<size_t> rls(nthreads, 0);
    size_t counter = 0;

    {
        size_t cell_size = cfg::get().buffer_size / (nthreads * nthreads * sizeof(Kmer) * 2);

        ParallelDeBruijn par_debruijn(nthreads, cell_size);
        while (!ParllelStreamEOF(streams)) {

            #pragma omp parallel num_threads(nthreads)
            {
                #pragma omp for reduction(+ : counter)
                for (size_t i = 0; i < nthreads; ++i) {

                    Read r;
                    io::IReader<Read>& stream = *streams[i];

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
                    par_debruijn.Dump(temp_maps[i], i);
                }
            }

        }
    }

    size_t total_kmers = 0;
    for (size_t i = 0; i < nthreads; ++i) {
        total_kmers += temp_maps[i].size();
    }

    //Merging into final map
    INFO("Merging final maps");
    debruijn.nodes().rehash(total_kmers);
    for (size_t i = 0; i < nthreads; ++i) {
        InsertIntoDebruijn<k, Graph, typename ParallelDeBruijn::destination_container_t>(debruijn, temp_maps[i]);
        temp_maps[i].erase(temp_maps[i].begin(), temp_maps[i].end());
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


template<size_t k, class Graph, class Read>
void ConstructGraph(std::vector<io::IReader<Read>* >& streams, Graph& g, EdgeIndex<k + 1, Graph>& index,
		SingleReadStream* contigs_stream = 0) {

    typedef SeqMap<k + 1, typename Graph::EdgeId> DeBruijn;
    typedef Seq<k + 1> Kmer;

    INFO("Constructing DeBruijn graph");
	
    DeBruijn& debruijn = index.inner_index();

    size_t rl = 0;
    if (streams.size() > 0) {
        rl = FillIterativeParallelIndex<k, Graph, Read>(streams, debruijn);
    }
    else if (streams.size() == 1) {
        rl = FillUsusalIndex<k, Graph, Read>(*streams.back(), debruijn);
    }

    io::SingleRead r;
    if (contigs_stream) {
        INFO("Adding contigs from previous K");
        while (!contigs_stream->eof()) {
            *contigs_stream >> r;
            Sequence s = r.sequence();
            debruijn.CountSequence(s);
        }
    }

    if (!cfg::get().ds.RL.is_initialized()) {
        INFO("Figured out: read length = " << rl);
        cfg::get_writable().ds.RL = rl;
    } else if (*cfg::get().ds.RL != rl) {
        WARN("In datasets.info, wrong RL is specified: " << cfg::get().ds.RL << ", not " << rl);
    }

    INFO("Condensing graph");
    DeBruijnGraphConstructor<k, Graph> g_c(debruijn);
    g_c.ConstructGraph(g, index);
    DEBUG("Graph condensed");
}

template<size_t k, class Read>
void ConstructGraphWithCoverage(std::vector<io::IReader<Read>* >& streams, Graph& g, EdgeIndex<k + 1, Graph>& index,
		SingleReadStream* contigs_stream = 0) {
	ConstructGraph<k>(streams, g, index, contigs_stream);
	FillCoverage<k, Read>(streams, g, index);
}

template<size_t k>
void ConstructGraphWithPairedInfo(graph_pack<ConjugateDeBruijnGraph, k>& gp,
		PairedInfoIndex<Graph>& paired_index, PairedReadStream& paired_stream,
		SingleReadStream* single_stream = 0,
		SingleReadStream* contigs_stream = 0) {
	UnitedStream united_stream(paired_stream);

	vector<SingleReadStream*> streams;
	if(!cfg::get().etalon_graph_mode) {
		streams.push_back(&united_stream);
	}
	if (single_stream) {
		streams.push_back(single_stream);
	}
	CompositeSingleReadStream reads_stream(streams);
	vector<SingleReadStream*> strs;
	strs.push_back(&reads_stream);
	ConstructGraphWithCoverage<k, io::SingleRead>(strs, gp.g, gp.index, reads_stream, contigs_stream);

	if (cfg::get().etalon_info_mode || cfg::get().etalon_graph_mode)
		FillEtalonPairedIndex<k>(paired_index, gp.g, gp.index, gp.kmer_mapper, gp.genome);
	else
		FillPairedIndex<k>(gp.g, gp.index, paired_index, paired_stream);
}

}

#endif /* GRAPH_CONSTRUCTION_HPP_ */

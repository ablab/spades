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

template<size_t k>
void FillCoverage(Graph& g, SingleReadStream& stream,
		EdgeIndex<k + 1, Graph>& index) {
	typedef SimpleSequenceMapper<k + 1, Graph> SequenceMapper;
	stream.reset();
	INFO("Counting coverage");
	SequenceMapper read_threader(g, index);
	g.coverage_index().FillIndex<SequenceMapper>(stream, read_threader);
	DEBUG("Coverage counted");
}

template<size_t k, class Graph> 
size_t FillParallelIndex(SeqMap<k + 1, typename Graph::EdgeId>& debruijn, SingleReadStream& reads_stream, size_t nthreads, size_t buf_size) {
        
    INFO("openmp " << omp_get_max_threads());

    typedef ParallelSeqMap<k + 1, typename Graph::EdgeId> ParallelDeBruijn;
    typedef Seq<k + 1> Kmer;

    ParallelDeBruijn par_debruijn(nthreads);

    INFO("Processing reads (takes a while)");
    size_t counter = 0;
    size_t rl = 0;
    
    Sequence** vector_seq = new Sequence*[nthreads];
    for (size_t i = 0; i < nthreads; ++i)
        vector_seq[i] = new Sequence[buf_size];
    std::vector<size_t> sizes(nthreads, 0);
    for (size_t i = 0; i<nthreads; ++i) VERIFY(sizes[i] == 0);

    io::SingleRead r;
    {
        perf_counter pc;

        while (!reads_stream.eof()) {

            //for (size_t i = 0; i < nthreads; ++i) 
                //vector_seq[i].reserve(buf_size);
            INFO("Filling Buffer");
            for (size_t j = 0; j < buf_size; ++j) {
                for(size_t i = 0; i < nthreads && !reads_stream.eof(); ++i) {
                    reads_stream >> r;
                    const Sequence& seq = r.sequence();
                    vector_seq[i][j] = seq;
                    rl = max(rl, seq.size());
                    sizes[i]++;
                    VERBOSE_POWER(++counter, " reads processed");
                }
            }
            INFO("Finished Filling");

            INFO("Threading Buffer");

            #pragma omp parallel num_threads(nthreads)
            {
                #pragma omp for
                for (size_t i = 0; i < nthreads; ++i) {
                    for (size_t j = 0; j < sizes[i]; ++j) {
                        par_debruijn.CountSequence(vector_seq[i][j], i);
                    }
                }
            }

            INFO("Finished Threading");
            std::fill(sizes.begin(), sizes.end(), 0);
        }
        INFO("Started Merging");

        for (size_t i = 0; i < nthreads; ++i) 
            delete[] vector_seq[i];
        delete[] vector_seq;

        par_debruijn.MergeMaps(debruijn.nodes());

        INFO("Finished Merging, size of map " << debruijn.nodes().size());

        INFO("Elapsed time: " << pc.time_ms());
    }

    INFO("DeBruijn graph constructed, " << counter << " reads used");

    return rl;
}

template<size_t k, class Graph>
void ConstructGraph(Graph& g, EdgeIndex<k + 1, Graph>& index,
		SingleReadStream& reads_stream, SingleReadStream* contigs_stream = 0) {

    typedef SeqMap<k + 1, typename Graph::EdgeId> DeBruijn;
    typedef Seq<k + 1> Kmer;

    size_t nthreads = 4;
    size_t buf_size = 10000000;
    
    cout << "# threads" << endl;
    cin >>  nthreads;

    cout << "buf size" << endl;
    cin >> buf_size;

    INFO("Constructing DeBruijn graph");
	
    DeBruijn& debruijn = index.inner_index();

    size_t rl = FillParallelIndex<k, Graph>(debruijn, reads_stream, nthreads, buf_size);

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

template<size_t k>
void ConstructGraphWithCoverage(Graph& g, EdgeIndex<k + 1, Graph>& index,
		SingleReadStream& reads_stream, SingleReadStream* contigs_stream = 0) {
	ConstructGraph<k>(g, index, reads_stream, contigs_stream);
	FillCoverage<k>(g, reads_stream, index);
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
	ConstructGraphWithCoverage<k>(gp.g, gp.index, reads_stream, contigs_stream);

	if (cfg::get().etalon_info_mode || cfg::get().etalon_graph_mode)
		FillEtalonPairedIndex<k>(paired_index, gp.g, gp.index, gp.kmer_mapper, gp.genome);
	else
		FillPairedIndex<k>(gp.g, gp.index, paired_index, paired_stream);
}

}

#endif /* GRAPH_CONSTRUCTION_HPP_ */

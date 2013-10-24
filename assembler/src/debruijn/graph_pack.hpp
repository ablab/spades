//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * graph_pack.hpp
 *
 *  Created on: Aug 18, 2011
 *      Author: sergey
 */

#pragma once

#include "omni/edge_labels_handler.hpp"
#include "omni/id_track_handler.hpp"
#include "omni/edges_position_handler.hpp"
#include "de/paired_info.hpp"
#include "debruijn_graph.hpp"
#include "config_struct.hpp"
#include "edge_index.hpp"
#include "genomic_quality.hpp"
#include "sequence_mapper.hpp"
#include "genomic_info.hpp"
#include "long_read_storage.hpp"
#include "detail_coverage.hpp"

namespace debruijn_graph {

/*KmerFree*//*KmerStoring*/
template<class Graph, class SeqType, class KmerEdgeIndex = DeBruijnEdgeIndex<KmerStoringDeBruijnEdgeIndex<Graph, SeqType>>>
struct graph_pack: private boost::noncopyable {
    typedef Graph graph_t;
    typedef SeqType seq_t;
    typedef EdgeIndex<graph_t, seq_t, KmerEdgeIndex> index_t;
    typedef omnigraph::de::PairedInfoIndicesT<Graph> PairedInfoIndicesT;
    typedef LongReadContainer<Graph> LongReadContainerT;

    size_t k_value;

    graph_t g;
    index_t index;
    IdTrackHandler<graph_t> int_ids;
    EdgesPositionHandler<graph_t> edge_pos;
//	PairedInfoIndex<graph_t> etalon_paired_index;
    KmerMapper<graph_t, seq_t> kmer_mapper;
    FlankingCoverage<graph_t> flanking_cov;
    PairedInfoIndicesT paired_indices;
    PairedInfoIndicesT clustered_indices;
    PairedInfoIndicesT scaffolding_indices;
    LongReadContainerT single_long_reads;

    GenomicInfo ginfo;
    Sequence genome;
	EdgeQuality<Graph> edge_qual;

    graph_pack(size_t k, const std::string &workdir, size_t lib_count,
                        Sequence genome = Sequence(), size_t single_gap = 0,
                        bool careful_labeling = false, bool use_inner_ids = false,
                        size_t flanking_range = 50)
            : k_value(k), g(k), index(g, k + 1, workdir),
              int_ids(g, use_inner_ids), edge_pos(g, (int) single_gap, careful_labeling),
              kmer_mapper(g, k + 1),
              flanking_cov(g, flanking_range),
              paired_indices(g, lib_count),
              clustered_indices(g, lib_count),
              scaffolding_indices(g, lib_count),
              single_long_reads(g, lib_count),
              genome(genome),
              edge_qual(g)
    { }

    void FillQuality() {
        edge_qual.Fill(index, kmer_mapper, genome);
    }

    //todo remove with usages after checking
    void ClearQuality() {
        edge_qual.clear();
    }

//    void FillFlankingCoverage() {
//        flanking_cov.Fill(index.inner_index());
//    }
};

typedef graph_pack<ConjugateDeBruijnGraph, runtime_k::RtSeq,
                   DeBruijnEdgeIndex<KmerFreeDeBruijnEdgeIndex<ConjugateDeBruijnGraph, runtime_k::RtSeq>>> conj_graph_pack;
typedef conj_graph_pack::index_t Index;

typedef conj_graph_pack::PairedInfoIndicesT PairedIndicesT;
typedef conj_graph_pack::LongReadContainerT LongReadContainerT;
typedef omnigraph::de::PairedInfoIndexT<ConjugateDeBruijnGraph> PairedIndexT;


} // namespace debruijn_graph

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
#include "sequence_mapper.hpp"

namespace debruijn_graph {

/*KmerFree*//*KmerStoring*/
template<class Graph, class SeqType, class KmerEdgeIndex = DeBruijnEdgeIndex<KmerStoringDeBruijnEdgeIndex<Graph, SeqType>>>
struct graph_pack: private boost::noncopyable {
    typedef Graph graph_t;
    typedef SeqType seq_t;
    typedef EdgeIndex<graph_t, seq_t, KmerEdgeIndex> index_t;
    typedef omnigraph::de::PairedInfoIndicesT<Graph> PairedInfoIndicesT;

    size_t k_value;

    graph_t g;
    index_t index;
    IdTrackHandler<graph_t> int_ids;
    EdgesPositionHandler<graph_t> edge_pos;
//	PairedInfoIndex<graph_t> etalon_paired_index;
    KmerMapper<graph_t, seq_t> kmer_mapper;
    PairedInfoIndicesT paired_indices;
    PairedInfoIndicesT clustered_indices;

    Sequence genome;

    explicit graph_pack(unsigned k, const std::string &workdir,
                        Sequence genome = Sequence(), size_t single_gap = 0,
                        bool careful_labeling = false, bool use_inner_ids = false)
            : k_value(k), g(k), index(g, k + 1, workdir),
              int_ids(g, use_inner_ids), edge_pos(g, (int) single_gap, careful_labeling),
              kmer_mapper(g, k + 1),
              paired_indices(g, cfg::get().ds.reads.lib_count()), clustered_indices(g, cfg::get().ds.reads.lib_count()),
              genome(genome)
    { }
};

typedef graph_pack<ConjugateDeBruijnGraph, runtime_k::RtSeq,
                   DeBruijnEdgeIndex<KmerFreeDeBruijnEdgeIndex<ConjugateDeBruijnGraph, runtime_k::RtSeq>>> conj_graph_pack;
typedef conj_graph_pack::index_t Index;
typedef conj_graph_pack::PairedInfoIndicesT PairedIndicesT;
typedef omnigraph::de::PairedInfoIndexT<ConjugateDeBruijnGraph> PairedIndexT;

} // namespace debruijn_graph

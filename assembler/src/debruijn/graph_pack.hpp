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
#include "graphio.hpp"
#include "mismatch_masker.hpp"

namespace debruijn_graph {

typedef PairedInfoIndexT<ConjugateDeBruijnGraph> PairedIndexT;

template<class Graph, class SeqType>
struct graph_pack: private boost::noncopyable {
	typedef Graph graph_t;
	typedef SeqType seq_t;

	size_t k_value;

	graph_t g;
	EdgeIndex<graph_t, seq_t> index;
	IdTrackHandler<graph_t> int_ids;
	EdgesPositionHandler<graph_t> edge_pos;
//	PairedInfoIndex<graph_t> etalon_paired_index;
	KmerMapper<graph_t, seq_t> kmer_mapper;
	Sequence genome;
	MismatchMasker<graph_t> mismatch_masker;

	//todo review params
	explicit graph_pack(size_t k, const std::string &workdir,
			Sequence const& genome = Sequence(), size_t single_gap = 0,
			bool careful_labeling = false, bool use_inner_ids = false) :
			k_value(k), g(k), index(g, k + 1, workdir), int_ids(g,
					use_inner_ids), edge_pos(g, single_gap, careful_labeling), kmer_mapper(
					g, k + 1), genome(genome), mismatch_masker(g) {
	}
};

typedef graph_pack<ConjugateDeBruijnGraph, runtime_k::RtSeq> conj_graph_pack;
typedef graph_pack<NonconjugateDeBruijnGraph, runtime_k::RtSeq> nonconj_graph_pack;

inline void Convert(const conj_graph_pack& gp1,
		const PairedInfoIndexT<conj_graph_pack::graph_t>& clustered_index1,
		nonconj_graph_pack& gp2,
		PairedInfoIndexT<nonconj_graph_pack::graph_t>& clustered_index2) {
	string conv_folder = path::append_path(cfg::get().output_root,
			"temp_conversion");
	make_dir(conv_folder);
	string p = path::append_path(conv_folder, "conj_graph");
	PrintWithClusteredIndex(p, gp1, clustered_index1);
	ScanWithClusteredIndex(p, gp2, clustered_index2);
	remove_dir(conv_folder);
}

} // namespace debruijn_graph


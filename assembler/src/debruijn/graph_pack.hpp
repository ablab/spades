//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
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
#include "omni/paired_info.hpp"
#include "omni/id_track_handler.hpp"
#include "omni/edges_position_handler.hpp"
#include "new_debruijn.hpp"
#include "config_struct.hpp"
#include "graphio.hpp"

namespace debruijn_graph {

typedef PairedInfoIndex<ConjugateDeBruijnGraph> paired_info_index;

template<class Graph, size_t k>
struct graph_pack : private boost::noncopyable {
	typedef Graph graph_t;

	enum {
		k_value = k
	};

	graph_t g;
	EdgeIndex<k + 1, graph_t> index;
	IdTrackHandler<graph_t> int_ids;
	EdgesPositionHandler<graph_t> edge_pos;
	PairedInfoIndex<graph_t> etalon_paired_index;
	KmerMapper<k + 1, graph_t> kmer_mapper;

	Sequence const& genome;

	explicit graph_pack(Sequence const& genome) :
	g(k),
	index(g)
	, int_ids (g)
	, edge_pos(g, cfg::get().pos.max_single_gap), etalon_paired_index(g, 0), kmer_mapper(g),
	genome(genome) {
	}

	graph_pack() :
	g(k),
	index(g)
	, int_ids (g)
	, edge_pos(g, cfg::get().pos.max_single_gap), etalon_paired_index(g, 0), kmer_mapper(g),
	genome(Sequence()) {
	}
};

typedef graph_pack<ConjugateDeBruijnGraph, K> conj_graph_pack;
typedef graph_pack<NonconjugateDeBruijnGraph, K> nonconj_graph_pack;

inline void Convert(const conj_graph_pack& gp1, const PairedInfoIndex<conj_graph_pack::graph_t>& clustered_index1,
		nonconj_graph_pack& gp2, PairedInfoIndex<nonconj_graph_pack::graph_t>& clustered_index2) {
	fs::path conv_folder = fs::path(cfg::get().output_root) / "temp_conversion";
	make_dir(conv_folder.string());

	fs::path p = conv_folder / "conj_graph";

	PrintWithClusteredIndex(p.string(), gp1, clustered_index1);

	ScanWithClusteredIndex(p.string(), gp2, clustered_index2);

	remove_all(conv_folder);
}

} // namespace debruijn_graph


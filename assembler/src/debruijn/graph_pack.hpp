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

template<class Graph>
struct graph_pack : private boost::noncopyable {
	typedef Graph graph_t;

	size_t k_value;

	graph_t g;
	EdgeIndex<graph_t> index;
	IdTrackHandler<graph_t> int_ids;
	EdgesPositionHandler<graph_t> edge_pos;
//	PairedInfoIndex<graph_t> etalon_paired_index;
	KmerMapper<graph_t> kmer_mapper;

	Sequence genome;

	explicit graph_pack(size_t k, Sequence const& genome = Sequence(), size_t single_gap = 0, bool careful_labeling = false, bool use_inner_ids = false) :
    k_value(k),
	g(k),
	index(g, k + 1)
	, int_ids(g, use_inner_ids)
	, edge_pos(g, single_gap, careful_labeling), kmer_mapper(g, k + 1),
	genome(genome) {
	}
};

typedef graph_pack<ConjugateDeBruijnGraph> conj_graph_pack;
typedef graph_pack<NonconjugateDeBruijnGraph> nonconj_graph_pack;

inline void Convert(const conj_graph_pack& gp1, const PairedInfoIndex<conj_graph_pack::graph_t>& clustered_index1,
		nonconj_graph_pack& gp2, PairedInfoIndex<nonconj_graph_pack::graph_t>& clustered_index2) {

    string conv_folder = path::append_path(cfg::get().output_root, "temp_conversion");
    make_dir(conv_folder);

    string p = path::append_path(conv_folder, "conj_graph");

    PrintWithClusteredIndex(p, gp1, clustered_index1);
    ScanWithClusteredIndex (p, gp2, clustered_index2);

    remove_dir(conv_folder);
}

} // namespace debruijn_graph


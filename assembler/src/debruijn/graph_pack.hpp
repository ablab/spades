/*
 * graph_pack.hpp
 *
 *  Created on: Aug 18, 2011
 *      Author: sergey
 */

#pragma once

#include "omni/edge_labels_handler.hpp"
#include "omni/paired_info.hpp"
#include "omni/ID_track_handler.hpp"
#include "omni/edges_position_handler.hpp"
#include "new_debruijn.hpp"
#include "config_struct.hpp"
#include "graphio.hpp"

namespace debruijn_graph {

typedef PairedInfoIndex<ConjugateDeBruijnGraph> paired_info_index;

struct conj_graph_pack {
	typedef ConjugateDeBruijnGraph graph_t;

	graph_t g;
	EdgeIndex<K + 1, graph_t> index;
	IdTrackHandler<graph_t> int_ids;
	EdgesPositionHandler<graph_t> edge_pos;
	PairedInfoIndex<graph_t> etalon_paired_index;
	KmerMapper<K + 1, graph_t> kmer_mapper;

	Sequence const& genome;

	conj_graph_pack(Sequence const& genome) :
	g(K),
	index(g)
	, int_ids (g)
	, edge_pos(g), etalon_paired_index(g, 0), kmer_mapper(g),
	genome(genome) {
	}
};

struct nonconj_graph_pack {
	typedef NonconjugateDeBruijnGraph graph_t;

	graph_t g;
	IdTrackHandler<graph_t> int_ids;
	EdgeIndex<K + 1, graph_t> index;
	EdgesPositionHandler<graph_t> edge_pos;
	PairedInfoIndex<graph_t> clustered_index;

	nonconj_graph_pack()
	: g (K)
	, int_ids (g)
	, index (g)
	, edge_pos (g)
	, clustered_index (g)
	{
	}

	nonconj_graph_pack(conj_graph_pack const& gp, PairedInfoIndex<ConjugateDeBruijnGraph> const& prev_clustered_index)
	: g (K)
	, int_ids (g)
	, index (g)
	, edge_pos (g)
	, clustered_index (g)
	{
		fs::path conv_folder = fs::path(cfg::get().output_root) / "temp_conversion";
		mkdir(conv_folder.string().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH | S_IWOTH);

		fs::path p = conv_folder / "conj_graph";

		// todo: make printGraph const to its arguments
		printGraph<conj_graph_pack::graph_t>(
				gp.g,
				gp.int_ids,
				p.string(),
				gp.edge_pos,
				(paired_info_index*) 0,
				&gp.etalon_paired_index,
				&prev_clustered_index);

		scanNCGraph<graph_t>(g, int_ids, p.string(), 0, edge_pos, 0, &clustered_index);

		rmdir(conv_folder.string().c_str());
	}
};

} // namespace debruijn_graph


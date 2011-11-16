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
	, edge_pos(g, cfg::get().pos.max_single_gap), etalon_paired_index(g, 0), kmer_mapper(g),
	genome(genome) {
	}
};

struct nonconj_graph_pack {
	typedef NonconjugateDeBruijnGraph graph_t;

	graph_t g;
	EdgeIndex<K + 1, graph_t> index;
	IdTrackHandler<graph_t> int_ids;
	EdgesPositionHandler<graph_t> edge_pos;
	PairedInfoIndex<graph_t> etalon_paired_index;
	KmerMapper<K + 1, graph_t> kmer_mapper;
	//todo extract from here
	PairedInfoIndex<graph_t> clustered_index;

	nonconj_graph_pack()
	: g (K)
	, index (g)
	, int_ids (g)
	, edge_pos (g, cfg::get().pos.max_single_gap)
	, etalon_paired_index(g)
	, kmer_mapper(g)
	, clustered_index (g)
	{
	}

	//todo make separate tool and make unique graph_pack
	nonconj_graph_pack(conj_graph_pack const& gp, PairedInfoIndex<ConjugateDeBruijnGraph> const& prev_clustered_index)
	: g (K)
	, index (g)
	, int_ids (g)
	, edge_pos (g, cfg::get().pos.max_single_gap)
	, etalon_paired_index(g)
	, kmer_mapper(g)
	, clustered_index (g)
	{
		fs::path conv_folder = fs::path(cfg::get().output_root) / "temp_conversion";
		make_dir(conv_folder.string());
		
		fs::path p = conv_folder / "conj_graph";

		ConjugateDataPrinter<conj_graph_pack::graph_t> printer(gp.g, gp.int_ids);
		PrintGraphPack(p.string(), printer, gp);
		PrintClusteredIndex(p.string(), printer, prev_clustered_index);

		NonconjugateDataScanner<graph_t> scanner(g, int_ids);
		ScanGraphPack(p.string(), scanner, *this);
		ScanClusteredIndex(p.string(), scanner, clustered_index);

		remove_all(conv_folder);
	}
};

} // namespace debruijn_graph


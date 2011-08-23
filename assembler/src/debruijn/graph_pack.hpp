/*
 * graph_pack.hpp
 *
 *  Created on: Aug 18, 2011
 *      Author: sergey
 */

#ifndef GRAPH_PACK_HPP_
#define GRAPH_PACK_HPP_

#include "omni/edge_labels_handler.hpp"
#include "omni/paired_info.hpp"
#include "omni/ID_track_handler.hpp"
#include "omni/edges_position_handler.hpp"
#include "new_debruijn.hpp"

namespace debruijn_graph {

template<size_t k>
struct ConjugateGraphPack {
	Graph& g;
	EdgeIndex<k + 1, Graph>& index;
	IdTrackHandler<Graph>& int_ids;
	EdgesPositionHandler<Graph>& edge_pos;
	KmerMapper<k + 1, Graph>& kmer_mapper;
	PairedInfoIndex<Graph>& paired_index;
	PairedInfoIndex<Graph>& etalon_paired_index;
	const Sequence& genome;

	ConjugateGraphPack(Graph& g, EdgeIndex<k + 1, Graph>& index
			, IdTrackHandler<Graph>& int_ids
			, EdgesPositionHandler<Graph>& edge_pos
			, KmerMapper<k + 1, Graph>& kmer_mapper
			, PairedInfoIndex<Graph>& paired_index
			, PairedInfoIndex<Graph>& clustered_index
			, PairedInfoIndex<Graph>& etalon_paired_index,
			const Sequence& genome) :
			g(g), index(index), int_ids(int_ids), edge_pos(edge_pos), kmer_mapper(kmer_mapper), paired_index(
					paired_index), etalon_paired_index(
					etalon_paired_index), genome(genome) {
	}
};

}
#endif /* GRAPH_PACK_HPP_ */

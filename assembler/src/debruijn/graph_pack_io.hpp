#pragma once

#include "graph_pack.hpp"
#include "graphio.hpp"

namespace debruijn_graph {

void PrintConjugateGraphPack(const string& file_name, conj_graph_pack& gp,
		PairedInfoIndex<Graph> const* paired_index = 0,
		PairedInfoIndex<Graph> const* clustered_index = 0) {
	printGraph(gp.g, gp.int_ids, file_name, gp.edge_pos, paired_index,
			&gp.etalon_paired_index, clustered_index, &gp.kmer_mapper);
}

void ScanConjugateGraphPack(const string& file_name, conj_graph_pack& gp,
		PairedInfoIndex<Graph>* paired_index = 0, PairedInfoIndex<Graph>* clustered_index = 0) {
	scanConjugateGraph(&gp.g, &gp.int_ids, file_name,
			paired_index, &gp.edge_pos, &gp.etalon_paired_index, clustered_index, &gp.kmer_mapper);
}

}

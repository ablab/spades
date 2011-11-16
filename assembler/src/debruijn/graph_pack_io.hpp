#pragma once

#include "graph_pack.hpp"
#include "graphio.hpp"

namespace debruijn_graph {

//void PrintConjugateFull(const string& file_name, conj_graph_pack& gp,
//		PairedInfoIndex<Graph> const* paired_index = 0,
//		PairedInfoIndex<Graph> const* clustered_index = 0) {
//	ConjugateDataPrinter<conj_graph_pack::graph_t> printer(gp.g, gp.int_ids);
//	PrintGraphPack(file_name, printer, gp);
//	if (paired_index) {
//		PrintPairedIndex(file_name, printer, *paired_index);
//	}
//	if (clustered_index) {
//		PrintClusteredIndex(file_name, printer, *clustered_index);
//	}
//}

//void ScanConjugateGraphPack(const string& file_name, conj_graph_pack& gp,
//		PairedInfoIndex<Graph>* paired_index = 0, PairedInfoIndex<Graph>* clustered_index = 0) {
//	scanConjugateGraph(&gp.g, &gp.int_ids, file_name,
//			paired_index, &gp.edge_pos, &gp.etalon_paired_index, clustered_index, &gp.kmer_mapper);
//}

}

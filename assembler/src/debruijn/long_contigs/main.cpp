/*
 * main.cpp
 *
 *  Created on: Jul 11, 2011
 *      Author: andrey
 */

#include "lc_launch.hpp"

using namespace long_contigs;
using namespace debruijn_graph;

DECL_PROJECT_LOGGER("d");

int main() {
	cfg::create_instance(debruijn_graph::cfg_filename);
	lc_cfg::create_instance(lc_cfg_filename);

	checkFileExistenceFATAL(lc_cfg_filename);
	checkFileExistenceFATAL(debruijn_graph::cfg_filename);

	Graph g(K);
	EdgeIndex<K + 1, Graph> index(g);
	IdTrackHandler<Graph> intIds(g);
	PairedInfoIndex<Graph> pairedIndex(g, 0);
	PairedInfoIndices pairedInfos;
	KmerMapper<K+1, Graph> mapper(g);


	Sequence genome = long_contigs::load_genome();

	LoadFromFile(lc_cfg::get().ds.graph_file, g, intIds, mapper);

	make_dir(cfg::get().output_dir);

	if (cfg::get().etalon_info_mode) {
		AddEtalonInfo<K>(g, index, mapper, genome, pairedInfos);
	} else {
		pairedInfos.clear();
		AddRealInfo<K>(g, index, intIds, pairedInfos, mapper, lc_cfg::get().use_new_metrics);

		if (!cfg::get().etalon_info_mode && lc_cfg::get().write_real_paired_info) {
			SavePairedInfo(g, intIds, pairedInfos, cfg::get().output_dir + lc_cfg::get().paired_info_file_prefix, !lc_cfg::get().use_new_metrics);
		}

		if (lc_cfg::get().paired_info_only) {
			//pairedInfos.clear();
			//AddRealInfo<K>(g, index, intIds, pairedInfos, mapper, !lc_cfg::get().use_new_metrics);

			//if (lc_cfg::get().write_real_paired_info) {
			//	SavePairedInfo(g, pairedInfos, intIds, output_dir + lc_cfg::get().paired_info_file_prefix, lc_cfg::get().use_new_metrics);
			//}

			return 0;
		}
	}

	return 0;
}


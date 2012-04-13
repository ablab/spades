/*
 * main.cpp
 *
 *  Created on: Jul 11, 2011
 *      Author: andrey
 */

#include "lc_launch.hpp"

#include "../resolved_pair_info.hpp"

using namespace long_contigs;
using namespace debruijn_graph;


int main() {
//	cfg::create_instance(debruijn_graph::cfg_filename);
	lc_cfg::create_instance(lc_cfg_filename);

	checkFileExistenceFATAL(lc_cfg_filename);
	checkFileExistenceFATAL(debruijn_graph::cfg_filename);

	std::cerr << "!!!\n";

    Sequence genome = long_contigs::load_genome();
	conj_graph_pack gp(genome);

	std::cerr << "!!!\n";

	Graph& g = gp.g;
	EdgeIndex<K + 1, Graph>& index = gp.index;
	IdTrackHandler<Graph>& intIds = gp.int_ids;
	PairedInfoIndex<Graph> pairedIndex(g, 0);
	PairedInfoIndices pairedInfos;
	KmerMapper<K+1, Graph>& mapper = gp.kmer_mapper;

	std::cerr << "!!!\n";

	LoadFromFile(lc_cfg::get().ds.graph_file, g, intIds, mapper);

    std::string output_dir = "data/debruijn/" + lc_cfg::get().dataset_name + "/K" + ToString(K) + "/" + MakeLaunchTimeDirName() + "/alternative_rr/";
	make_dir(output_dir);

	if (lc_cfg::get().etalon_mode) {
		AddEtalonInfo<K>(g, index, mapper, genome, pairedInfos);
	} else {
		pairedInfos.clear();
		AddRealInfo<K>(g, index, intIds, pairedInfos, mapper, lc_cfg::get().use_new_metrics);

//		if (!cfg::get().etalon_info_mode && lc_cfg::get().write_real_paired_info) {
//			SavePairedInfo(g, intIds, pairedInfos, output_dir + lc_cfg::get().paired_info_file_prefix, !lc_cfg::get().use_new_metrics);
//		}

		if (lc_cfg::get().paired_info_for_resolved) {
		    Graph rg(K);
		    IdTrackHandler<Graph> r_intIds(rg);
		    KmerMapper<K+1, Graph> r_mapper(rg);
		    PairedInfoIndex<Graph> r_pairedIndex(rg, 0);
		    PairedInfoIndices r_pairedInfos;

            LoadFromFile(lc_cfg::get().resolved_graph, rg, r_intIds, r_mapper);

//            ResolvedGraphPairInfoCounter<Graph> resolved_graph_paired_info_counter(
//                    g, pairedInfos[0].pairedInfoIndex, rg, labels_after);

//            resolved_graph_paired_info_counter.FillResolvedGraphPairedInfo(
//                    r_pairedIndex);

            r_pairedInfos.push_back(PairedInfoIndexLibrary(rg, pairedInfos[0].readSize, pairedInfos[0].insertSize, pairedInfos[0].is_delta, pairedInfos[0].deDelta, pairedInfos[0].var,
                    &r_pairedIndex));

//            SavePairedInfo(rg, r_intIds, r_pairedInfos, output_dir + lc_cfg::get().paired_info_file_prefix + "resolved_", lc_cfg::get().use_new_metrics);

		    return 0;
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

	resolve_repeats_ml(gp, pairedInfos, output_dir, lc_cfg::get().params);

	return 0;
}


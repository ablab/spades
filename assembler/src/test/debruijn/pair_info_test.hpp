#pragma once

#include <boost/test/unit_test.hpp>
#include "test_utils.hpp"
#include "debruijn_stats.hpp"
#include "graphio.hpp"

namespace debruijn_graph {

BOOST_AUTO_TEST_SUITE(pair_info_tests)

BOOST_AUTO_TEST_CASE( EstimationFunctionalTest ) {
	conj_graph_pack gp;
	paired_info_index paired_index(gp.g);
	paired_info_index clustered_index(gp.g);

	size_t insert_size = 220;

//	put path here
	ScanAll("./data/debruijn/QUAKE_CROPPED_400K/K55/latest/saves/distance_estimation", gp, paired_index, clustered_index);
	paired_info_index etalon_paired_index(gp.g);
	FillAndCorrectEtalonPairedInfo(etalon_paired_index, gp, paired_index, 220, cfg::get().ds.RL, cfg::get().de.delta);
	INFO("Counting clustered info stats");
	EdgeQuality<Graph> edge_qual(gp.g, gp.index, gp.kmer_mapper, gp.genome);
	EstimationQualityStat<Graph> estimation_stat(gp.g, gp.int_ids, edge_qual, paired_index,
			clustered_index, etalon_paired_index, false);
	estimation_stat.Count();
	BOOST_CHECK_LE(estimation_stat.fpr(), 0.05);
	BOOST_CHECK_LE(estimation_stat.fnr(), 0.05);
}

BOOST_AUTO_TEST_SUITE_END()

}

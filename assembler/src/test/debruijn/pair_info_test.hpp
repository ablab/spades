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

//	put path here
	ScanAll("", gp, paired_index, clustered_index);
	paired_info_index etalon_paired_index(gp.g);
	FillAndCorrectEtalonPairedInfo(etalon_paired_index, gp, paired_index, output_folder);
	INFO("Counting clustered info stats");
	EdgeQuality<Graph> edge_qual(gp.g, gp.index, gp.kmer_mapper, gp.genome);
	EstimationQualityStat<Graph> estimation_stat(gp.g, gp.int_ids, edge_qual, paired_index,
			clustered_index, etalon_paired_index);
	estimation_stat.Count();
	estimation_stat.SaveStats();
}

BOOST_AUTO_TEST_SUITE_END()

}

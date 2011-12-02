#pragma once

#include <boost/test/unit_test.hpp>
#include "test_utils.hpp"
#include "debruijn_stats.hpp"
#include "graphio.hpp"

namespace debruijn_graph {

BOOST_AUTO_TEST_SUITE(pair_info_tests)

//todo get rid of magic constants
BOOST_AUTO_TEST_CASE( EstimationFunctionalTest ) {
	string genome_str;
	string genome_filename = "./data/input/E.Coli.K12.MG1655/MG1655-K12.fasta.gz";
	checkFileExistenceFATAL(genome_filename);
	io::Reader<io::SingleRead> genome_stream(genome_filename);
	io::SingleRead full_read;
	genome_stream >> full_read;
	genome_str = full_read/*.GetSequenceString().substr(0,
			400000)*/;

	Sequence genome(genome_str);
	conj_graph_pack gp(genome);
	paired_info_index paired_index(gp.g);
	paired_info_index clustered_index(gp.g);

//	put path here
	ScanAll("./data/debruijn/QUAKE_FULL/K55/latest/saves/distance_estimation", gp, paired_index, clustered_index);
	paired_info_index etalon_paired_index(gp.g);
	FillAndCorrectEtalonPairedInfo(etalon_paired_index, gp, paired_index, 220, 100, 10);
	INFO("Counting clustered info stats");
	EdgeQuality<Graph> edge_qual(gp.g, gp.index, gp.kmer_mapper, gp.genome);
	EstimationQualityStat<Graph> estimation_stat(gp.g, gp.int_ids, edge_qual, paired_index,
			clustered_index, etalon_paired_index);
	estimation_stat.Count();
	INFO("fpr " << estimation_stat.fpr());
	INFO("fnr " << estimation_stat.fnr());
	BOOST_CHECK_LE(estimation_stat.fpr(), 0.05);
	BOOST_CHECK_LE(estimation_stat.fnr(), 0.05);
}

BOOST_AUTO_TEST_SUITE_END()

}

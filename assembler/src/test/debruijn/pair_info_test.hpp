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
	string genome_filename = "./data/input/E.coli/MG1655-K12.fasta.gz";
	checkFileExistenceFATAL(genome_filename);
	io::Reader<io::SingleRead> genome_stream(genome_filename);
	io::SingleRead full_read;
	genome_stream >> full_read;
//	genome_str = full_read/*.GetSequenceString().substr(0,
//			400000)*/;

	Sequence genome(full_read.sequence());
	conj_graph_pack gp(genome);
	paired_info_index paired_index(gp.g);
	paired_info_index clustered_index(gp.g);

//	put path here
//	ScanAll("./data/debruijn/ECOLI_IS480_QUAKE/K55/latest/saves/distance_estimation", gp, paired_index, clustered_index);
//	ScanAll("./data/debruijn/ECOLI_IS480_QUAKE/K55/no_clust/saves/distance_estimation", gp, paired_index, clustered_index);
	ScanAll("./data/debruijn/ECOLI_IS480_QUAKE/K55/clust/saves/distance_estimation", gp, paired_index, clustered_index);
	paired_info_index etalon_paired_index(gp.g);
	FillAndCorrectEtalonPairedInfo(etalon_paired_index, gp, paired_index, /*is*/500, /*rl*/100, /*delta*/25);
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

void CheckSymmetry(const paired_info_index &paired_index) {
	for(auto it = paired_index.begin(); it != paired_index.end(); ++it) {
		auto info = *it;
		auto symmetric_info = paired_index.GetEdgePairInfo(info[0].second, info[0].first);
		BOOST_CHECK_EQUAL(info.size(), symmetric_info.size());
		for(size_t i = 0; i < info.size(); i++) {
			BOOST_CHECK_EQUAL(info[i].first, symmetric_info[symmetric_info.size() - 1 - i].second);
			BOOST_CHECK_EQUAL(info[i].second, symmetric_info[symmetric_info.size() - 1 - i].first);
			BOOST_CHECK_EQUAL(info[i].weight, symmetric_info[symmetric_info.size() - 1 - i].weight);
			BOOST_CHECK_EQUAL(info[i].d, -symmetric_info[symmetric_info.size() - 1 - i].d);
		}
	}
}

void CheckRCSymmetry(const ConjugateDeBruijnGraph &graph, const paired_info_index &paired_index) {
	for(auto it = paired_index.begin(); it != paired_index.end(); ++it) {
		auto info = *it;
		auto symmetric_info = paired_index.GetEdgePairInfo(graph.conjugate(info[0].second), graph.conjugate(info[0].first));
		BOOST_CHECK_EQUAL(info.size(), symmetric_info.size());
		for(size_t i = 0; i < info.size(); i++) {
			BOOST_CHECK_EQUAL(info[i].first, graph.conjugate(symmetric_info[i].second));
			BOOST_CHECK_EQUAL(info[i].second, graph.conjugate(symmetric_info[i].first));
			BOOST_CHECK_EQUAL(info[i].weight, symmetric_info[i].weight);
			BOOST_CHECK_EQUAL(info[i].d, -graph.length(info[i].first) + graph.length(info[i].second) + symmetric_info[i].d);
		}
	}
}

//BOOST_AUTO_TEST_CASE( CheckPairInfoSimmetry ) {
//	INFO("CheckPairInfoSimmetry started");
//	string genome_str;
//	string genome_filename = "./data/input/E.Coli.K12.MG1655/MG1655-K12.fasta.gz";
//	checkFileExistenceFATAL(genome_filename);
//	io::Reader<io::SingleRead> genome_stream(genome_filename);
//	io::SingleRead full_read;
//	genome_stream >> full_read;
//	Sequence genome(full_read.sequence());
//	conj_graph_pack gp(genome);
//	paired_info_index paired_index(gp.g);
//	paired_info_index clustered_index(gp.g);
//
//	ScanAll("./data/debruijn/ECOLI_IS220_QUAKE/K55/latest/saves/distance_estimation", gp, paired_index, clustered_index);
//	paired_info_index etalon_paired_index(gp.g);
//	CheckSymmetry(paired_index);
//	CheckRCSymmetry(gp.g, paired_index);
//	CheckSymmetry(clustered_index);
//	CheckRCSymmetry(gp.g, clustered_index);
//	INFO("CheckPairInfoSimmetry finished");
//}

BOOST_AUTO_TEST_SUITE_END()

}

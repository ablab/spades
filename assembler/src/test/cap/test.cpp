//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "compare_standard.hpp"
#include "logger/log_writers.hpp"
#include "graphio.hpp"
#include <boost/test/unit_test.hpp>

#include "comparison_utils.hpp"
#include "diff_masking.hpp"
#include "repeat_masking.hpp"
#include "assembly_compare.hpp"
#include "test_utils.hpp"

::boost::unit_test::test_suite*	init_unit_test_suite( int, char* [] )
{
	//logging::create_logger("", logging::L_DEBUG);
	//logging::__logger()->add_writer(make_shared<logging::console_writer>());
  logging::logger *log = logging::create_logger("", logging::L_INFO);
	log->add_writer(std::make_shared<logging::console_writer>());
	logging::attach_logger(log);

    using namespace ::boost::unit_test;
	char module_name [] = "cap_test";

    assign_op( framework::master_test_suite().p_name.value, basic_cstring<char>(module_name), 0 );

	return 0;
}

#include "lseq_test.hpp"

namespace cap {

inline void CheckDiffs(const string& actual_prefix, const string& etalon_prefix, bool exact_match = true) {
    string comparison_type_string;
    if (exact_match) {
        comparison_type_string = "exact match";
    } else {
        comparison_type_string = "graph isomorphism";
    }
    INFO("Checking differences for graphs: " + comparison_type_string);

    if (exact_match) {
        string suffixes[4] = {"grp", "clr", "sqn", "pos"};
        for (size_t i = 0; i < 4; ++i) {
            BOOST_CHECK_MESSAGE(CheckFileDiff(actual_prefix + "." + suffixes[i], etalon_prefix + "." + suffixes[i]),
                    "Check for suffix " + suffixes[i] + " failed");
        }
    } else {
        // TODO load both graphs and compare them as graphs.
        // In order to incapsulate reading of greaph and coloring.
        BOOST_CHECK_MESSAGE(CheckColoredGraphIsomorphism(actual_prefix, etalon_prefix), "GRAPHS DIFFER");
    }
}


template<size_t k, size_t K, class BuildSeq>
inline void LoadAndRunBPG(const string& filename, const string& output_dir, const string& etalon_root,
		const string& example_id = "", bool regenerate_etalon = false) {
	string add_saves_path = "";
	if (regenerate_etalon) {
		remove_dir(etalon_root);
		add_saves_path = etalon_root;
		make_dir(add_saves_path);
	}
	using boost::property_tree::ptree;
	ptree pt;
	read_xml(filename, pt);
//		size_t example_cnt = 0;
	BOOST_FOREACH (const ptree::value_type& example, pt.get_child("examples")) {
		const ptree& genomes_node = example.second;
		//todo change name
		string n = genomes_node.get<string>("<xmlattr>.n");
		if (example_id != "" && example_id != n) {
			INFO("Ignoring example " << n);
			continue;
		}

		vector<vector<io::SingleRead>> genomes;
		BOOST_FOREACH (const ptree::value_type& genome,	genomes_node) {
			if (genome.first == "genome") {
				const ptree& contigs_node = genome.second;
				vector<io::SingleRead> contigs;
				size_t contig_cnt = 0;
				BOOST_FOREACH (const ptree::value_type& contig,	contigs_node) {
					contigs.push_back(io::SingleRead("contig_" + ToString(contig_cnt++), contig.second.data()));
				}
				genomes.push_back(contigs);
			}
		}
		INFO("--------------------------------------------");
		INFO("Processing example " << n);
		VERIFY(genomes.size() == 2);
		io::VectorReader<io::SingleRead> stream_1(genomes[0]);
		io::VectorReader<io::SingleRead> stream_2(genomes[1]);
    RunBPComparison<k, K, BuildSeq>(
        stream_1,
        stream_2,
        "genome_0_",
        "genome_1_", /*refine*/
        true, /*untangle*/
        false,
        output_dir + "example_" + n /*ToString(++example_cnt)*/
          + "/", /*detailed_output*/true
          ,5, Sequence(), (add_saves_path != "") ? add_saves_path + n : "");
    if (!regenerate_etalon) {
      CheckDiffs(output_dir + "example_" + n + "/saves/colored_split_graph", etalon_root + n, false);
    }
  }
}

/*
BOOST_AUTO_TEST_CASE( ColoringStringsTest ) {
    ColorGenerator::instance().GenerateColors(20);
    for (size_t i = 0; i < 20; ++i) {
        INFO("Color #" << i << ": " << ColorGenerator::instance().GetIthColor(i));
    }
}
*/

BOOST_AUTO_TEST_CASE( SyntheticExamplesTestsRtSeq ) {
  utils::TmpFolderFixture _("tmp");
	make_dir("bp_graph_test");
	LoadAndRunBPG<15, 25, runtime_k::RtSeq>("./src/test/cap/tests/synthetic/tests.xml",
			"bp_graph_test/simulated_common/", "./src/test/cap/tests/synthetic/etalon/", "", false);
	remove_dir("bp_graph_test");
}

//BOOST_AUTO_TEST_CASE( SyntheticExamplesTestsLSeq ) {
//
//  utils::TmpFolderFixture _("tmp");
//	make_dir("bp_graph_test");
//
//	LoadAndRunBPG<15, 25, LSeq>("./src/test/cap/tests/synthetic/tests.xml",
//			"bp_graph_test/simulated_common/", "./src/test/cap/tests/synthetic/etalon/", "", false);
//	remove_dir("bp_graph_test");
//}

/*
BOOST_AUTO_TEST_CASE( SyntheticExamplesWithErrorsTests ) {
  utils::TmpFolderFixture _("tmp");
	make_dir("bp_graph_test");
	LoadAndRunBPG<15, 25, runtime_k::RtSeq>("./src/test/cap/tests/synthetic_with_err/tests2.xml",
			"bp_graph_test/simulated_common_err/", "./src/test/cap/tests/synthetic_with_err/etalon/", "1_err", false);
	LoadAndRunBPG<15, 25, LSeq>("./src/test/cap/tests/synthetic_with_err_lseq/tests2.xml",
			"bp_graph_test/simulated_common_err/", "./src/test/cap/tests/synthetic_with_err_lseq/etalon/", "1_err_lseq", false);
	remove_dir("bp_graph_test");
}
*/
// todo parts for more tests
//	Sequence genome = ReadGenome("data/input/E.coli/MG1655-K12.fasta.gz");
//	INFO("Running comparison against mutated genome");
//	RunBPComparison<17, 250>(genome, IntroduceMutations(genome, 0.01), "init", "mut"
//			, /*refine*/true, /*untangle*/false, "bp_graph_test/mutated_ref/", /*detailed*/false);

//	INFO("Running comparison against genome with reversals");

//	RunBPComparison<25, 250>(genome, IntroduceReversals(genome, 10, 1000, 2000), "init", "rev"
//			, /*refine*/false, /*untangle*/false, "bp_graph_test/reversaled_ref/", /*detailed*/false);

//	INFO("Running comparison against mutated genome with reversals");
//	RunBPComparison<25, 250>(genome, IntroduceMutations(IntroduceReversals(genome, 10, 1000, 2000), 0.01), "init", "mut_rev"
//			, /*refine*/true, /*untangle*/true, "bp_graph_test/reversaled_mut_ref/", /*detailed*/false);

//	typedef graph_pack<ConjugateDeBruijnGraph, 55> gp_t;
//	INFO("Running comparison against repeat graph contigs");
//	RunBPComparison<25, 250>(genome, RepeatGraphEdges<gp_t>(genome), "init", "repeat_g_cont"
//			, /*refine*/false, /*untangle*/false, "bp_graph_test/repeat_graph_edges_ref/", /*detailed*/false);
//
//	INFO("Running comparison against reversaled repeat graph contigs");
//	RunBPComparison<25, 250>(genome, RepeatGraphEdges<gp_t>(IntroduceReversals(genome, 10, 1000, 10000)), "init", "rev_repeat_g_cont"
//			, /*refine*/false, /*untangle*/false, "bp_graph_test/rev_repeat_graph_edges_ref/", /*detailed*/false);
}

//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "standard_base.hpp"
#include "logger/log_writers.hpp"
#include "graphio.hpp"
#include <boost/test/unit_test.hpp>

#include "comparison_utils.hpp"
#include "assembly_compare.hpp"
#include "diff_masking.hpp"
#include "repeat_masking.hpp"

#include "tests.hpp"

::boost::unit_test::test_suite*	init_unit_test_suite( int, char* [] )
{
	logging::create_logger("", logging::L_DEBUG);
	logging::__logger()->add_writer(make_shared<logging::console_writer>());

    using namespace ::boost::unit_test;
	char module_name [] = "debruijn_tools";

    assign_op( framework::master_test_suite().p_name.value, basic_cstring<char>(module_name), 0 );

	return 0;
}

namespace compare {

//BOOST_AUTO_TEST_CASE( RefVSAssemblyComparison ) {
//	static const size_t k = 55;
//	static const size_t K = 55;
//	Sequence ref = ReadGenome("/home/snurk/MRSA/USA300_FPR3757.fasta");
//	io::Reader contig_stream("/home/snurk/MRSA/MRSA_RCH_I56.fasta");
//	string folder = "mrsa_comp/RCH_I56/";
//	make_dir(folder);
//	RunBPComparison<k, K>(
//			ref,
//			contig_stream,
//			"ref",
//			"assembly",
//			true/*refine*/,
//			false/*untangle*/,
//			folder,
//			true/*detailed_output*/,
//			20);
//}

//BOOST_AUTO_TEST_CASE( TwoAssemblyComparison ) {
//	static const size_t k = 19;
//	static const size_t K = 55;
////	static const size_t K = 57;
////	static const size_t K = 53;
//
////	io::Reader stream_1("/home/snurk/gingi/2.fasta");
////	io::Reader stream_2("/home/snurk/gingi/3.fasta");
//
//	io::Reader stream_1("/home/anton/idba_compare/idba.fasta");
//	io::Reader stream_2("/home/anton/idba_compare/hammer21_dis_tuned_simpl_try_improve.fasta");
//	string ref = "/home/anton/idba_compare/MG1655-K12.fasta";
//	string folder = "/home/anton/idba_compare/hammer21_dis_tuned_simpl_vs_idba/";
////	string folder = "assembly_comp/gingi_new_3_vs_jeff/";
//	make_dir(folder);
//
//	RunBPComparison<k, K>(
//		stream_1,
//		stream_2,
//		"idba_",
//		"k21ts_",
//		true/*refine*/,
//		false/*untangle*/,
//		folder,
//		true/*detailed_output*/,
//		5/*delta*/,
//		ReadGenome(ref));
//}

//BOOST_AUTO_TEST_CASE( TwoAssemblyComparison ) {
//	static const size_t k = 19;
//	static const size_t K = 55;
////	static const size_t K = 57;
////	static const size_t K = 53;
//
////	io::Reader stream_1("/home/snurk/gingi/2.fasta");
////	io::Reader stream_2("/home/snurk/gingi/3.fasta");
//
////		io::Reader stream_1("/home/snurk/gingi/PGINGIVALIS_LANE2_BH.fasta");
////	io::Reader stream_1("/home/snurk/gingi/PGINGIVALIS_LANE3_BH.fasta");
//	io::Reader stream_2("/home/snurk/gingi/jeff.fasta");
//
////	io::Reader stream_2("/home/snurk/gingi/PGINGIVALIS_LANE2_BH.fasta");
//
//	io::Reader stream_1("/home/snurk/gingi/PGINGIVALIS_LANE3_BH.fasta");
////	io::Reader stream_2("/home/snurk/gingi/lane2_evsc.fasta");
//
////	string folder = "assembly_comp/gingi_new_3_vs_new_2/";
//	string folder = "assembly_comp/gingi_new_3_vs_jeff/";
//	make_dir(folder);
//
//	RunBPComparison<k, K>(
//		stream_1,
//		stream_2,
////		"2",
////		"jeff",
//		"3_new_",
////		"2_new_",
//		"jeff_",
//		true/*refine*/,
//		false/*untangle*/,
//		folder,
//		true/*detailed_output*/);
//}

//BOOST_AUTO_TEST_CASE( AssemblyRefComparison ) {
//	static const size_t k = 21;
//	static const size_t K = 201/*55*//*201*/;
////	static const size_t K = 57;
////	static const size_t K = 53;
//
////	io::Reader stream_1("/home/snurk/gingi/2.fasta");
////	io::Reader stream_2("/home/snurk/gingi/3.fasta");
//
//	io::Reader stream_1("/home/snurk/Dropbox/gingi/jeff.fasta");
//	io::Reader stream_2("/home/snurk/Dropbox/gingi/TDC60.fasta");
//	string ref = "/home/snurk/Dropbox/gingi/TDC60.fasta";
//
////	string folder = "assembly_comp/gingi_jeff_vs_tdc60_55/";
//	string folder = "assembly_comp/gingi_jeff_vs_tdc60_201/";
//	make_dir(folder);
//
//	RunBPComparison<k, K>(
//		stream_1,
//		stream_2,
//		"jeff_",
//		"tdc_",
//		true/*refine*/,
//		false/*untangle*/,
//		folder,
//		true/*detailed_output*/,
//		5/*delta*/,
//		ReadGenome(ref));
//}

//BOOST_AUTO_TEST_CASE( IDBA_vs_SPADES ) {
//	static const size_t k = 19;
//	static const size_t K = 55;
////	static const size_t K = 57;
////	static const size_t K = 53;
//
////	io::Reader stream_1("/home/snurk/gingi/2.fasta");
////	io::Reader stream_2("/home/snurk/gingi/3.fasta");
//
//	io::Reader stream_1("/home/snurk/idba_comp/idba-contig-100.fa");
//	io::Reader stream_2("/home/snurk/idba_comp/k21nodiscard.fasta");
//	string ref = "/home/snurk/idba_comp/MG1655-K12.fasta";
//
//	string folder = "/home/snurk/idba_comp/results/";
//	make_dir(folder);
//
//	RunBPComparison<k, K>(
//		stream_1,
//		stream_2,
//		"idba_",
//		"bh21_",
//		true/*refine*/,
//		false/*untangle*/,
//		folder,
//		true/*detailed_output*/,
//		5/*delta*/,
//		ReadGenome(ref));
//}

//BOOST_AUTO_TEST_CASE( TwoStrainComparisonWR ) {
//	make_dir("bp_graph_test");
//	INFO("Running comparison of two strains");
//	pair<Sequence, Sequence> genomes = CorrectGenomes<55>(CorrectGenomes<21>(ReadGenome("data/input/E.coli/MG1655-K12.fasta.gz")
//			, ReadGenome("data/input/E.coli/DH10B-K12.fasta")), 200);
//	INFO("Genomes ready");
//
//	CompareGenomes<701>(genomes.first, genomes.second, "bp_graph_test/two_strain_comp_wr/");
//	INFO("Finished");
//}

//inline void StrainComparisonWOR(const string& strain_1, const string& strain_2, const string& output_folder) {
//	make_dir("bp_graph_test");
//	INFO("Running comparison of two strains");
//	pair<Sequence, Sequence> genomes = CorrectGenomes<55>(TotallyClearGenomes<55>(CorrectGenomes<21>(ReadGenome(strain_1)
//			, ReadGenome(strain_2))), 30);
////	genomes = TotallyClearGenomes<701>(genomes);
//	VERIFY(CheckNoRepeats<301>(genomes.first));
//	VERIFY(CheckNoRepeats<301>(genomes.second));
//	INFO("Genomes ready");
//
//	CompareGenomes<701>(genomes.first, genomes.second, output_folder);
//}

//BOOST_AUTO_TEST_CASE( TwoStrainComparisonWOR ) {
//	StrainComparisonWOR("data/input/E.coli/MG1655-K12.fasta.gz"
//, "data/input/E.coli/DH10B-K12.fasta", "bp_graph_test/two_strain_comp_wo_repeats/");
//}

//BOOST_AUTO_TEST_CASE( CompareAllMRSA ) {
//	string mrsa_root = "/home/snurk/MRSA/more_strains/";
//	ifstream stream;
//	stream.open(mrsa_root + "list.txt");
//	string s1;
//	string s2;
//	VERIFY(!stream.eof());
//	stream >> s1;
//	while (!stream.eof()) {
//		stream >> s2;
//		StrainComparisonWOR(mrsa_root + s1 + ".fasta", mrsa_root + s2 + ".fasta"
//				, mrsa_root + "results/" +  s1 + "_vs_" + s2 + "/");
//	}
//	stream.close();
//}

//
//BOOST_AUTO_TEST_CASE( TwoStrainComparisonFirstWOR ) {
//	make_dir("bp_graph_test");
//	INFO("Running comparison of two strains");
//	pair<Sequence, Sequence> genomes = CorrectGenomes<21>(ReadGenome("data/input/E.coli/MG1655-K12.fasta.gz")
//			, ReadGenome("data/input/E.coli/DH10B-K12.fasta"));
//	genomes = CorrectGenomes<55>(genomes, 200);
//	genomes.first = ClearGenome<55>(genomes.first);
//	genomes = CorrectGenomes<55>(genomes, 200);
//	VERIFY(CheckNoRepeats<301>(genomes.first));
//	INFO("Genomes ready");
//
//	CompareGenomes<701>(genomes.first, genomes.second, "bp_graph_test/two_strain_comp_first_wo_repeats/");
//}

//BOOST_AUTO_TEST_CASE( StrainVSRepeatGraphComparison ) {
//	static const size_t repeat_clearing_k = 55;
//	static const size_t repeat_graph_k = 101;
//	static const size_t refining_k1 = 25;
//	static const size_t refining_k2 = 151;
//	static const size_t bp_k = 301;
//
//	make_dir("bp_graph_test");
//	INFO("Running comparison of strain vs repeat graph for other strain");
//	Sequence genome1 = ReadGenome("data/input/E.coli/MG1655-K12.fasta.gz");
//	Sequence genome2 = ReadGenome("data/input/E.coli/DH10B-K12.fasta");
//
//	typedef graph_pack<ConjugateDeBruijnGraph, repeat_graph_k> repeat_gp_t;
//	vector<Sequence> repeat_graph_edges = RepeatGraphEdges<repeat_gp_t>(genome2);
//
//	pair<Sequence, vector<Sequence>> refined_data1 = RefineData<refining_k2>(RefineData<refining_k1>(
//			make_pair(genome1, repeat_graph_edges)));
//
//	pair<Sequence, vector<Sequence>> cleared = Clear<repeat_clearing_k>(refined_data1.first, refined_data1.second);
//
//	pair<Sequence, vector<Sequence>> refined_data2 = RefineData<bp_k>(RefineData<refining_k2>(RefineData<refining_k1>(cleared)));
//
//	io::VectorReader<io::SingleRead> final_stream1(
//			io::SingleRead("first", refined_data2.first.str()));
//	io::VectorReader<io::SingleRead> final_stream2(MakeReads(refined_data2.second));
//
//	typedef graph_pack<ConjugateDeBruijnGraph, bp_k> comparing_gp_t;
//	INFO("Running assembly comparer");
//	AssemblyComparer<comparing_gp_t> comparer(final_stream1, final_stream2, "strain1", "strain2", /*untangle*/false);
//	comparer.CompareAssemblies("bp_graph_test/strain_vs_repeat_graph_comp/", /*detailed_output*/true);
//	INFO("Finished");
//}
//
//BOOST_AUTO_TEST_CASE( StrainVSRepeatGraphComparison2 ) {
//	static const size_t repeat_clearing_k = 55;
//	static const size_t repeat_graph_k = 201;
//	static const size_t refining_k1 = 25;
//	static const size_t refining_k2 = 151;
//	static const size_t bp_k = 201;
//
//	make_dir("bp_graph_test");
//	INFO("Running comparison of strain vs repeat graph for other strain");
//	Sequence genome1 = ReadGenome("data/input/E.coli/MG1655-K12.fasta.gz");
//	Sequence genome2 = ReadGenome("data/input/E.coli/DH10B-K12.fasta");
//
//	typedef graph_pack<ConjugateDeBruijnGraph, repeat_graph_k> repeat_gp_t;
//	vector<Sequence> repeat_graph_edges = RepeatGraphEdges<repeat_gp_t>(genome2);
//
//	pair<Sequence, vector<Sequence>> refined_data1 = RefineData<refining_k2>(RefineData<refining_k1>(
//			make_pair(genome1, repeat_graph_edges)));
//
////	pair<Sequence, vector<Sequence>> cleared = Clear<repeat_clearing_k>(refined_data1.first, refined_data1.second);
//	Sequence cleared = ClearGenome<repeat_clearing_k>(refined_data1.first);
//
//	pair<Sequence, vector<Sequence>> refined_data2 = RefineData<bp_k>(RefineData<refining_k2>(RefineData<refining_k1>(make_pair(cleared, refined_data1.second))));
//
//	io::VectorReader<io::SingleRead> final_stream1(
//			io::SingleRead("first", refined_data2.first.str()));
//	io::VectorReader<io::SingleRead> final_stream2(MakeReads(refined_data2.second));
//
//	typedef graph_pack<ConjugateDeBruijnGraph, bp_k> comparing_gp_t;
//	INFO("Running assembly comparer");
//	AssemblyComparer<comparing_gp_t> comparer(final_stream1, final_stream2, "strain1", "strain2", /*untangle*/false);
//	comparer.CompareAssemblies("bp_graph_test/strain_vs_repeat_graph_comp2/", /*detailed_output*/true);
//	INFO("Finished");
//}

//BOOST_AUTO_TEST_CASE( ThreadingContigsOverGraph ) {
//	typedef graph_pack<ConjugateDeBruijnGraph, 55> gp_t;
//	io::EasyReader base_contigs("/home/anton/gitrep/algorithmic-biology/assembler/data/tmp/andrew_nurk.fasta");
//	io::EasyReader other_contigs("/home/anton/gitrep/algorithmic-biology/assembler/data/tmp/velvet-sc.fasta");
//	string base_saves = "/home/anton/gitrep/algorithmic-biology/assembler/data/debruijn/LBOUILLONII_QUAKE/saves/simplified_graph";
//	string output_dir = "bul_comparison/ac";
//	make_dir(output_dir);
//	ThreadAssemblies<gp_t>(base_saves, base_contigs, "spades", other_contigs, "velvet", output_dir);
//}

//BOOST_AUTO_TEST_CASE( GenerateGraphFragment ) {
//	std::string input_path = "./data/debruijn/HMP_LANE_3_0/K55/latest/saves/simplified_graph";
//	std::string output_path = "./src/test/debruijn/graph_fragments/topology_ec/iter_unique_path";
//	size_t split_threshold = 230;
//	int int_edge_id = 5573881;
//	graph_pack<ConjugateDeBruijnGraph, 55> gp;
//	ScanGraphPack(input_path, gp);
//	//prints only basic graph structure
//	PrintGraphComponentContainingEdge(output_path, gp.g,
//			split_threshold, gp.int_ids, int_edge_id);
//
//	//long way to write to dot file
//	Graph g(55);
//	IdTrackHandler<Graph> int_ids(g);
//	ScanBasicGraph(output_path, g, int_ids);
//	total_labeler_graph_struct graph_struct(g, &int_ids, (const EdgesPositionHandler<Graph>*)0);
//	total_labeler tot_lab(&graph_struct);
//	WriteToDotFile(g,
//			tot_lab, output_path + ".dot",
//			"mygraph", Path<EdgeId>(), Path<EdgeId>());
//}

}

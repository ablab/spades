#define BOOST_TEST_MODULE debruijn_tools

#include <boost/test/unit_test.hpp>
#include "assembly_compare.hpp"
#include "comparison_utils.hpp"

DECL_PROJECT_LOGGER("dtls")

namespace debruijn_graph {

//BOOST_AUTO_TEST_CASE( BreakPointGraph ) {
//	static const size_t k = 19;
//	static const size_t K = 201;
////    	io::EasyReader stream1("/home/snurk/assembly_compare/geba_0001_vsc.fasta.gz");
////    	io::EasyReader stream2("/home/snurk/assembly_compare/geba_0001_spades.fasta.gz");
//	//todo split N's
////	io::EasyReader stream1("/home/sergey/assembly_compare/geba_0002_allpaths.fasta.gz");
////	io::EasyReader stream2("/home/sergey/assembly_compare/geba_0002_spades.fasta.gz");
//
////	io::EasyReader stream1("/home/anton/gitrep/algorithmic-biology/assembler/data/PGINGIVALIS_LANE2_BH_split.fasta.gz");
////	io::EasyReader stream2("/home/anton/gitrep/algorithmic-biology/assembler/data/input/P.gingivalis/TDC60.fasta");
//// 	comparer.CompareAssemblies(stream1, stream2, "spades_", "ref_");
//
//	io::Reader<io::SingleRead> stream_1("/home/snurk/gingi/gingi_lane2_it.fasta");
////	io::Reader<io::SingleRead> stream_2("/home/snurk/gingi/lane2_evsc.fasta");
//	io::Reader<io::SingleRead> stream_2("/home/snurk/gingi/MDA2_clc_new.fasta");
//
//	RunBPComparison<k, K>(
//		stream_1,
//		stream_2,
//		"spades",
//		"evsc",
//		true/*refine*/,
//		true/*untangle*/,
//		"assembly_compare/",
//		true/*detailed_output*/);
//}

template<size_t k, size_t K>
inline void LoadAndRunBPG(const string& filename, const string& output_dir,
		const string& example_id = "") {
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
		RunBPComparison<k, K>(
			stream_1,
			stream_2,
			"genome_0_",
			"genome_1_", /*refine*/
			true, /*untangle*/
			true,
			output_dir + "example_" + n /*ToString(++example_cnt)*/
					+ "/", /*detailed_output*/false);
	}
}

BOOST_AUTO_TEST_CASE( TwoStrainComparisonWR ) {
	make_dir("bp_graph_test");
	INFO("Running comparison of two strains");
	pair<Sequence, Sequence> genomes = CorrectGenomes<21>(ReadGenome("data/input/E.coli/MG1655-K12.fasta.gz")
			, ReadGenome("data/input/E.coli/DH10B-K12.fasta"));
	INFO("Genomes ready");

	io::VectorReader<io::SingleRead> stream1(
			io::SingleRead("first", genomes.first.str()));
	io::VectorReader<io::SingleRead> stream2(
			io::SingleRead("second", genomes.second.str()));
	typedef graph_pack</*Nonc*/ConjugateDeBruijnGraph, 301> comparing_gp_t;
	INFO("Running assembly comparer");
	AssemblyComparer<comparing_gp_t> comparer(stream1, stream2, "strain1", "strain2", /*untangle*/false);
	comparer.CompareAssemblies("bp_graph_test/two_strain_comp/", /*detailed_output*/true);
	INFO("Finished");
}

BOOST_AUTO_TEST_CASE( TwoStrainComparisonWOR ) {
	make_dir("bp_graph_test");
	INFO("Running comparison of two strains");
	pair<Sequence, Sequence> genomes = CorrectGenomes<55>(TotallyClearGenomes<55>(CorrectGenomes<21>(ReadGenome("data/input/E.coli/MG1655-K12.fasta.gz")
			, ReadGenome("data/input/E.coli/DH10B-K12.fasta"))), 30);
	VERIFY(CheckNoRepeats<301>(genomes.first));
	VERIFY(CheckNoRepeats<301>(genomes.second));
	INFO("Genomes ready");

	io::VectorReader<io::SingleRead> stream1(
			io::SingleRead("first", genomes.first.str()));
	io::VectorReader<io::SingleRead> stream2(
			io::SingleRead("second", genomes.second.str()));
	typedef graph_pack</*Nonc*/ConjugateDeBruijnGraph, 301> comparing_gp_t;
	INFO("Running assembly comparer");
	AssemblyComparer<comparing_gp_t> comparer(stream1, stream2, "strain1", "strain2", /*untangle*/false);
	comparer.CompareAssemblies("bp_graph_test/two_strain_comp_wo_repeats/", /*detailed_output*/true);
	INFO("Finished");
}

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

BOOST_AUTO_TEST_CASE( BreakPointGraphTests ) {
//	make_dir("bp_graph_test");
//	INFO("Running simulated examples");
//	LoadAndRunBPG<7, 25>("/home/snurk/assembly_compare/tests2.xml",
//			"bp_graph_test/simulated_common/");
//
//	INFO("Running simulated examples with introduced errors");
//	LoadAndRunBPG<7, 25>("/home/snurk/assembly_compare/tests2.xml",
//			"bp_graph_test/simulated_common_err/", "1_err");
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


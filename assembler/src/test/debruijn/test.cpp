#define BOOST_TEST_MODULE debruijn_test

#include "graphio.hpp"
#include <iostream>
#include "logging.hpp"
#include "test_utils.hpp"
#include "assembly_compare.hpp"
#include "io/splitting_wrapper.hpp"


//headers with tests
#include "debruijn_graph_test.hpp"
#include "simplification_test.hpp"
//#include "pair_info_test.hpp"

DECL_PROJECT_LOGGER("dt")

namespace debruijn_graph {

	template<class gp_t>
	inline void ThreadAssemblies(const string& base_saves, ContigStream& base_assembly, const string& base_prefix
			, ContigStream& assembly_to_thread, const string& to_thread_prefix, const string& output_dir) {
		typedef typename gp_t::graph_t Graph;
		gp_t gp;
//		ConstructGraph<gp_t::k_value, Graph>(gp.g, gp.index, base_assembly);
		ScanGraphPack(base_saves, gp);
		base_assembly.reset();
		FillPos(gp, base_assembly, base_prefix);
		FillPos(gp, assembly_to_thread, to_thread_prefix);

		EdgePosGraphLabeler<Graph> pos_labeler(gp.g, gp.edge_pos);
		StrGraphLabeler<Graph> str_labeler(gp.g);
		CompositeLabeler<Graph> labeler(pos_labeler, str_labeler);

		NewExtendedSequenceMapper<gp_t::k_value + 1, Graph> mapper(gp.g, gp.index, gp.kmer_mapper);

		assembly_to_thread.reset();
		io::SingleRead read;
		while (!assembly_to_thread.eof()) {
			assembly_to_thread >> read;
			WriteComponentsAlongPath(gp.g, labeler, output_dir + read.name() + ".dot"
					, read.name(), /*split_edge_length*/400, mapper.MapSequence(read.sequence())
					, Path<typename Graph::EdgeId>(), Path<typename Graph::EdgeId>(), true);
		}
	}

BOOST_AUTO_TEST_CASE( BreakPointGraph ) {
	AssemblyComparer<graph_pack<NonconjugateDeBruijnGraph, 55>> comparer;
//    	io::EasyReader stream1("/home/snurk/assembly_compare/geba_0001_vsc.fasta.gz");
//    	io::EasyReader stream2("/home/snurk/assembly_compare/geba_0001_spades.fasta.gz");
	//todo split N's
//	io::EasyReader stream1("/home/sergey/assembly_compare/geba_0002_allpaths.fasta.gz");
//	io::EasyReader stream2("/home/sergey/assembly_compare/geba_0002_spades.fasta.gz");

//	io::EasyReader stream1("/home/anton/gitrep/algorithmic-biology/assembler/data/PGINGIVALIS_LANE2_BH_split.fasta.gz");
//	io::EasyReader stream2("/home/anton/gitrep/algorithmic-biology/assembler/data/input/P.gingivalis/TDC60.fasta");
// 	comparer.CompareAssemblies(stream1, stream2, "spades_", "ref_");

//	io::Reader<io::SingleRead> raw_reader_1("/home/sergey/assembly_compare/PGINGIVALIS_LANE2_BH_split.fasta.gz");
//	io::Reader<io::SingleRead> raw_reader_2("/home/sergey/assembly_compare/TDC60.fasta");
//	io::SplittingWrapper filtered_reader_1(raw_reader_1);
//	io::SplittingWrapper filtered_reader_2(raw_reader_2);
//	io::RCReaderWrapper<io::SingleRead> rc_reader_1(filtered_reader_1);
//	io::RCReaderWrapper<io::SingleRead> rc_reader_2(filtered_reader_2);
//	comparer.CompareAssemblies(rc_reader_1, rc_reader_2, "spades_", "ref_");
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


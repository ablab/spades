#define BOOST_TEST_MODULE debruijn_tools

#include <boost/test/unit_test.hpp>
#include "graphio.hpp"
#include <iostream>
#include "logging.hpp"
#include "assembly_compare.hpp"
#include "io/splitting_wrapper.hpp"
#include "io/vector_reader.hpp"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

DECL_PROJECT_LOGGER("dtls")

namespace debruijn_graph {

	template<size_t k, size_t K>
	inline void RunBPComparison(ContigStream& raw_stream1, ContigStream& raw_stream2
			, const string& name1, const string& name2, bool untangle, const string& output_folder
			, bool detailed_output = true) {
		io::SplittingWrapper stream1(raw_stream1);
		io::SplittingWrapper stream2(raw_stream2);

//		typedef graph_pack<ConjugateDeBruijnGraph, k> refining_gp_t;
//		refining_gp_t refining_gp;
//		ConstructGPForRefinement(refining_gp, stream1, stream2);
//

		//todo turn refinement on!!!
//		ContigRefiner<refining_gp_t> refined_stream1(stream1, refining_gp);
//		ContigRefiner<refining_gp_t> refined_stream2(stream2, refining_gp);

		typedef graph_pack</*Nonc*/ConjugateDeBruijnGraph, K> comparing_gp_t;
		AssemblyComparer<comparing_gp_t> comparer(/*refined_*/stream1, /*refined_*/stream2, name1, name2, untangle);
		comparer.CompareAssemblies(output_folder, detailed_output);
	}

//BOOST_AUTO_TEST_CASE( BreakPointGraph ) {
//	static const size_t k = 55;
//	static const size_t K = 101;
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
//	io::Reader<io::SingleRead> raw_reader_1("/home/anton/gitrep/algorithmic-biology/assembler/data/PGINGIVALIS_LANE2_BH_split.fasta.gz");
//	io::Reader<io::SingleRead> raw_reader_2("/home/anton/gitrep/algorithmic-biology/assembler/data/input/P.gingivalis/TDC60.fasta");
//
//	RunBPComparison<k, K>(raw_reader_1, raw_reader_2, "spades_", "ref_", true, "assembly_comparison/");
//}

	template <size_t k, size_t K>
	inline void LoadAndRunBPG(const string& filename, const string& output_dir, const string& example_id = "") {
		using boost::property_tree::ptree;
		ptree pt;
		read_xml(filename, pt);
//		size_t example_cnt = 0;
		BOOST_FOREACH (const ptree::value_type& example,
				pt.get_child("examples")) {
			const ptree& genomes_node = example.second;
			//todo change name
			string n = genomes_node.get<string>("<xmlattr>.n");
			if (example_id != "" && example_id != n) {
				INFO("Ignoring example " << n);
				continue;
			}

			vector<vector<io::SingleRead>> genomes;
			BOOST_FOREACH (const ptree::value_type& genome,
							genomes_node) {
				if (genome.first == "genome") {
					const ptree& contigs_node = genome.second;
					vector<io::SingleRead> contigs;
					size_t contig_cnt = 0;
					BOOST_FOREACH (const ptree::value_type& contig,
									contigs_node) {
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
			RunBPComparison<k, K>(stream_1, stream_2
					, "genome_0_", "genome_1_", /*untangle*/true
					, output_dir + "example_" + n /*ToString(++example_cnt)*/ + "/", /*detailed_output*/true);
		}
	}

BOOST_AUTO_TEST_CASE( BreakPointGraphTests ) {
	//todo now disabled: refinement, untangling!!!
	LoadAndRunBPG<25, 25>("/home/snurk/assembly_compare/tests.xml", "bp_graph_test/"/*, "6"*/);
}

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
		make_dir(output_dir + read.name());
		WriteComponentsAlongPath(gp.g, labeler, output_dir + read.name() + "/.dot"
				, /*split_edge_length*/400, mapper.MapSequence(read.sequence())
				, Path<typename Graph::EdgeId>(), Path<typename Graph::EdgeId>(), true);
	}
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

//Prints only basic graph structure!!!
inline void PrintGraphComponentContainingEdge(const string& file_name, const Graph& g,
		size_t split_edge_length, const IdTrackHandler<Graph>& int_ids,
		int int_edge_id) {
	LongEdgesInclusiveSplitter<Graph> inner_splitter(g, split_edge_length);

	VERIFY_MSG(int_ids.ReturnEdgeId(int_edge_id) != NULL, "Couldn't find edge with id = " << int_edge_id);

	AnyEdgeContainFilter<Graph> filter(g, int_ids.ReturnEdgeId(int_edge_id));
	FilteringSplitterWrapper<Graph> splitter(inner_splitter, filter);
	vector<vector<VertexId>> components;
	while (!splitter.Finished()) {
		components.push_back(splitter.NextComponent());
	}
	VERIFY(components.size() == 1);
	ConjugateDataPrinter<Graph> printer(g, components.front().begin(),
			components.front().end(), int_ids);
	PrintBasicGraph<Graph>(file_name, printer);
}

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


#define BOOST_TEST_MODULE debruijn_tools

#include <boost/test/unit_test.hpp>
#include "graphio.hpp"
#include "xmath.h"
#include <iostream>
#include "logging.hpp"
#include "assembly_compare.hpp"
#include "io/splitting_wrapper.hpp"
#include "io/vector_reader.hpp"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>

DECL_PROJECT_LOGGER("dtls")

namespace debruijn_graph {

inline double uniform_01() {
	static boost::mt19937 rng(43);
	static boost::uniform_01<boost::mt19937> zeroone(rng);
	return zeroone();
}

inline bool event_happened(double rate) {
	return ls(uniform_01(), rate);
}

inline int rand_int(size_t min, size_t max) {
	static boost::mt19937 rng(43);
	boost::uniform_int<> un_int(min, max);
	boost::variate_generator<boost::mt19937&, boost::uniform_int<> >
	         die(rng, un_int);
	return die();
}

inline char switch_nucl(char n) {
	VERIFY(is_nucl(n));
	return nucl((dignucl(n) + rand_int(1,3)) % 4);
}

inline Sequence IntroduceReversal(const Sequence& s, size_t min_len, size_t max_len) {
	VERIFY(s.size() > min_len);
	//inclusive
	size_t start = rand_int(0, s.size() - min_len);
	size_t len = rand_int(min_len, std::min(max_len, s.size() - start));
	//exclusive
	size_t end = start + len;
	INFO("Reversing fragment of length " << len << " from " << start << " to " << end);
	return s.Subseq(0, start) + !s.Subseq(start, end) + s.Subseq(end);
}

inline Sequence IntroduceReversals(const Sequence& s, size_t rev_count, size_t min_len, size_t max_len) {
	Sequence res = s;
	for (size_t i = 0; i < rev_count; ++i) {
		res = IntroduceReversal(res, min_len, max_len);
	}
	return res;
}

template<class gp_t>
inline void ConstructRepeatGraph(gp_t& gp) {
	io::VectorReader<io::SingleRead> stream(
			io::SingleRead("genome", gp.genome.str()));
	io::RCReaderWrapper<io::SingleRead> rc_stream(stream);
	ConstructGraph<gp_t::k_value, typename gp_t::graph_t>(gp.g, gp.index, rc_stream);
}

template<class Graph>
inline vector<Sequence> EdgesSequences(const Graph& g) {
	vector<Sequence> res;
	set<EdgeId> edges;
	for (auto it = g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
		if (edges.find(*it) == edges.end()) {
			res.push_back(g.EdgeNucls(*it));
			edges.insert(g.conjugate(*it));
		}
	}
	return res;
}

template<class gp_t>
inline vector<Sequence> RepeatGraphEdges(const Sequence& genome) {
	typedef typename gp_t::graph_t Graph;
	typedef typename Graph::EdgeId EdgeId;

	gp_t gp(genome);
	ConstructRepeatGraph(gp);
	return EdgesSequences(gp.g);
}

inline Sequence ReadGenome(const string& filename) {
	checkFileExistenceFATAL(filename);
	io::Reader<io::SingleRead> genome_stream(filename);
	io::SingleRead genome;
	genome_stream >> genome;
	return genome.sequence();
}

inline Sequence IntroduceMutations(const Sequence& s, double rate) {
	VERIFY(ge(rate, 0.) && ls(rate, 1.0));
	string as_str = s.str();
	for (size_t i = 0; i < s.size(); ++i) {
		if (event_happened(rate)) {
			as_str[i] = switch_nucl(as_str[i]);
		}
	}
	return Sequence(as_str);
}

template<size_t k, size_t K>
inline void RunBPComparison(ContigStream& raw_stream1,
		ContigStream& raw_stream2, const string& name1, const string& name2,
		bool refine, bool untangle, const string& output_folder,
		bool detailed_output = true) {
	io::SplittingWrapper stream1(raw_stream1);
	io::SplittingWrapper stream2(raw_stream2);

	typedef graph_pack</*Nonc*/ConjugateDeBruijnGraph, K> comparing_gp_t;
	if (refine) {
		typedef graph_pack<ConjugateDeBruijnGraph, k> refining_gp_t;
		refining_gp_t refining_gp;
		ConstructGPForRefinement(refining_gp, stream1, stream2);

		ContigRefiner<refining_gp_t> refined_stream1(stream1, refining_gp);
		ContigRefiner<refining_gp_t> refined_stream2(stream2, refining_gp);

		AssemblyComparer<comparing_gp_t> comparer(refined_stream1,
				refined_stream2, name1, name2, untangle);
		comparer.CompareAssemblies(output_folder, detailed_output);
	} else {
		AssemblyComparer<comparing_gp_t> comparer(stream1,
				stream2, name1, name2, untangle);
		comparer.CompareAssemblies(output_folder, detailed_output);
	}
}

template<size_t k, size_t K>
inline void RunBPComparison(const Sequence& ref,
		ContigStream& stream, const string& name1, const string& name2,
		bool refine, bool untangle, const string& output_folder,
		bool detailed_output = true) {
	io::VectorReader<io::SingleRead> ref_stream(
			io::SingleRead(name1, ref.str()));
	RunBPComparison<k, K>(ref_stream, stream, name1, name2,
			refine, untangle, output_folder,
			detailed_output);
}

template<size_t k, size_t K>
inline void RunBPComparison(const Sequence& s1,
		const Sequence& s2, const string& name1, const string& name2,
		bool refine, bool untangle, const string& output_folder,
		bool detailed_output = true) {
	io::VectorReader<io::SingleRead> stream(
			io::SingleRead(name2, s2.str()));
	RunBPComparison<k, K>(s1, stream, name1, name2,
			refine, untangle, output_folder,
			detailed_output);
}

const vector<io::SingleRead> MakeReads(const vector<Sequence>& ss) {
	vector<io::SingleRead> ans;
	for (size_t i = 0; i < ss.size(); ++i) {
		ans.push_back(io::SingleRead("read_" + ToString(i), ss[i].str()));
	}
	return ans;
}

template<size_t k, size_t K>
inline void RunBPComparison(const Sequence& ref,
		const vector<Sequence>& contigs, const string& name1, const string& name2,
		bool refine, bool untangle, const string& output_folder,
		bool detailed_output = true) {
	io::VectorReader<io::SingleRead> stream(MakeReads(contigs));
	RunBPComparison<k, K>(ref, stream, name1, name2,
			refine, untangle, output_folder,
			detailed_output);
}

BOOST_AUTO_TEST_CASE( BreakPointGraph ) {
	static const size_t k = 19;
	static const size_t K = 201;
//    	io::EasyReader stream1("/home/snurk/assembly_compare/geba_0001_vsc.fasta.gz");
//    	io::EasyReader stream2("/home/snurk/assembly_compare/geba_0001_spades.fasta.gz");
	//todo split N's
//	io::EasyReader stream1("/home/sergey/assembly_compare/geba_0002_allpaths.fasta.gz");
//	io::EasyReader stream2("/home/sergey/assembly_compare/geba_0002_spades.fasta.gz");

//	io::EasyReader stream1("/home/anton/gitrep/algorithmic-biology/assembler/data/PGINGIVALIS_LANE2_BH_split.fasta.gz");
//	io::EasyReader stream2("/home/anton/gitrep/algorithmic-biology/assembler/data/input/P.gingivalis/TDC60.fasta");
// 	comparer.CompareAssemblies(stream1, stream2, "spades_", "ref_");

	io::Reader<io::SingleRead> stream_1("/home/snurk/gingi/gingi_lane2_it.fasta");
//	io::Reader<io::SingleRead> stream_2("/home/snurk/gingi/lane2_evsc.fasta");
	io::Reader<io::SingleRead> stream_2("/home/snurk/gingi/MDA2_clc_new.fasta");

	RunBPComparison<k, K>(
		stream_1,
		stream_2,
		"spades",
		"evsc",
		true/*refine*/,
		true/*untangle*/,
		"assembly_compare/",
		true/*detailed_output*/);
}

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

//BOOST_AUTO_TEST_CASE( BreakPointGraphTests ) {
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

//
//	INFO("Running comparison against other strain");
//	RunBPComparison<21, 301>(genome, ReadGenome("data/input/E.coli/DH10B-K12.fasta"), "init", "other_strain"
//			, /*refine*/true, /*untangle*/true, "bp_graph_test/other_strain_comp/", /*detailed*/false);
//}

template<class gp_t>
inline void ThreadAssemblies(const string& base_saves,
		ContigStream& base_assembly, const string& base_prefix,
		ContigStream& assembly_to_thread, const string& to_thread_prefix,
		const string& output_dir) {
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

	NewExtendedSequenceMapper<gp_t::k_value + 1, Graph> mapper(gp.g, gp.index,
			gp.kmer_mapper);

	assembly_to_thread.reset();
	io::SingleRead read;
	while (!assembly_to_thread.eof()) {
		assembly_to_thread >> read;
		make_dir(output_dir + read.name());
		WriteComponentsAlongPath(gp.g, labeler,
				output_dir + read.name() + "/.dot", /*split_edge_length*/400,
				mapper.MapSequence(read.sequence()),
				Path<typename Graph::EdgeId>(), Path<typename Graph::EdgeId>(),
				true);
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
//todo rewrite with normal splitter usage instead of filtering
inline void PrintGraphComponentContainingEdge(const string& file_name,
		const Graph& g, size_t split_edge_length,
		const IdTrackHandler<Graph>& int_ids, int int_edge_id) {
	LongEdgesInclusiveSplitter<Graph> inner_splitter(g, split_edge_length);

	VERIFY_MSG(int_ids.ReturnEdgeId(int_edge_id) != NULL,
			"Couldn't find edge with id = " << int_edge_id);

	AnyEdgeContainFilter<Graph> filter(g, int_ids.ReturnEdgeId(int_edge_id));
	FilteringSplitterWrapper<Graph> splitter(inner_splitter, filter);
	vector<vector<VertexId>> components;
	while (!splitter.Finished()) {
		components.push_back(splitter.NextComponent());
	}VERIFY(components.size() == 1);
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


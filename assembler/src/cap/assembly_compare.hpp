//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "standard.hpp"
#include "utils.hpp"
#include "graph_pack.hpp"
#include "simple_tools.hpp"
#include "omni/omni_utils.hpp"
#include "debruijn_stats.hpp"
#include "io/splitting_wrapper.hpp"
#include "io/multifile_reader.hpp"
#include <boost/algorithm/string/predicate.hpp>

#include "coloring.hpp"
#include "colored_graph_construction.hpp"
#include "untangling.hpp"
#include "assembly_problem_detection.hpp"
#include "stats.hpp"
#include "visualization.hpp"

namespace cap {

//class RCSplittingStream: public io::DelegatingReaderWrapper<io::SingleRead> {
//private:
//	io::SplittingWrapper filtered_reader_;
//	io::RCReaderWrapper<io::SingleRead> stream_;
//public:
//	RCSplittingStream(ContigStream &base_stream) :
//			filtered_reader_(base_stream), stream_(filtered_reader_) {
//		Init(stream_);
//	}
//};

//todo finish later
//template<class gp_t1, class gp_t2>
//void ConvertToBPGraphPack(const gp_t1& gp
//		, const ColorHandler<typename gp_t1::graph_t>& coloring
//		, gp_t2& bp_gp) {
//	string tmp_dir = "/home/snurk/tmp/";
//	string filename = tmp_dir + "tmp";
//	make_dir(tmp_dir);
//	PrintGraphPack(filename, gp);
//	typename ScannerTraits<typename gp_t2::graph_t>::Scanner scanner(bp_gp.g,
//				bp_gp.int_ids);
//	ScanBasicGraph(filename, scanner);
//	scanner.loadPositions(filename, bp_gp.edge_pos);
//	//
//}

template<class Graph>
class GapsRemover {
	typedef typename Graph::VertexId VertexId;
	Graph& g_;
	const ColorHandler<Graph>& coloring_;
	const edge_type gap_color_;
	const size_t length_bound_;
public:
	GapsRemover(Graph& g, const ColorHandler<Graph>& coloring,
			edge_type gap_color, size_t length_bound) :
			g_(g), coloring_(coloring), gap_color_(gap_color), length_bound_(
					length_bound) {

	}

	void RemoveGaps() {
		for (auto it = g_.SmartEdgeBegin(); !it.IsEnd(); ++it) {
			if (coloring_.Color(*it) == gap_color_
					&& g_.length(*it) <= length_bound_
					&& g_.CanCompressVertex(g_.EdgeStart(*it))
					&& g_.CanCompressVertex(g_.EdgeEnd(*it))) {
				VertexId start = g_.EdgeStart(*it);
				VertexId end = g_.EdgeEnd(*it);
				if (!g_.RelatedVertices(start, end)) {
					g_.CompressVertex(end);
				}
				g_.CompressVertex(start);
			}
		}
	}
};

template<class Graph>
void SimplifyGraph(Graph& g, size_t br_delta) {
	debruijn_config::simplification::bulge_remover br_config;
	br_config.max_bulge_length_coefficient = 2;
	br_config.max_coverage = 1000.;
	br_config.max_relative_coverage = 1.2;
	br_config.max_delta = br_delta;
	br_config.max_relative_delta = 0.1;
	INFO("Removing bulges");
	RemoveBulges(g, br_config);

//		debruijn_config::simplification::tip_clipper tc;
//		tc.max_coverage = 1000;
//		tc.max_relative_coverage = 1000;
//		tc.max_tip_length_coefficient = 6;
//		ClipTips(gp.g, tc, 10 * gp.g.k());
}

template<class gp_t>
void ConstructColoredGraph(gp_t& gp,
		ColorHandler<typename gp_t::graph_t>& coloring,
		vector<ContigStream*>& streams, int br_delta = -1) {
	typedef typename gp_t::graph_t Graph;
	const size_t k = gp_t::k_value;
	typedef NewExtendedSequenceMapper<k + 1, Graph> Mapper;

	INFO("Constructing de Bruijn graph for k=" << k);

	//dirty hack because parallel construction uses cfg::get!!!
	io::MultifileReader<Contig> stream(streams);
	ConstructGraph<k, Graph>(gp.g, gp.index, stream);

	//TODO do we still need it?
	if (br_delta > 0)
		SimplifyGraph(gp.g, br_delta);

	ColoredGraphConstructor<Graph, Mapper> colored_graph_constructor(gp.g,
			coloring, *MapperInstance < gp_t > (gp));
	colored_graph_constructor.ConstructGraph(streams);

	INFO("Filling contig positions");
	for (auto it = streams.begin(); it != streams.end(); ++it) {
		ContigStream& stream = **it;
		stream.reset();
		FillPos(gp, stream);
	}
}

template<class gp_t>
class AssemblyComparer {
private:
	typedef typename gp_t::graph_t Graph;
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef NewExtendedSequenceMapper<gp_t::k_value + 1, Graph> Mapper;

//	io::IReader<io::SingleRead> &stream1_;
//	io::IReader<io::SingleRead> &stream2_;
	gp_t gp_;
	ColorHandler<Graph> coloring_;
	io::RCReaderWrapper<io::SingleRead> rc_stream1_;
	io::RCReaderWrapper<io::SingleRead> rc_stream2_;
	string name1_;
	PrefixAddingReaderWrapper stream1_;
	string name2_;
	PrefixAddingReaderWrapper stream2_;
	bool untangle_;

	template<class gp_t2>
	void UniversalSaveGP(
			const gp_t2& gp/*, const GraphColorer<typename gp_t2::graph_t> coloring*/,
			const string& filename) {
		typename PrinterTraits<Graph>::Printer printer(gp.g, gp.int_ids);
		INFO("Saving graph to " << filename);
		printer.saveGraph(filename);
		printer.saveEdgeSequences(filename);
		printer.savePositions(filename, gp.edge_pos);
//		SaveColoring(gp.g
//				, gp.int_ids
//				, coloring
//				, filename);

//		LengthIdGraphLabeler<Graph> labeler(gp.g);
//		WriteSimple(gp.g, labeler, filename + ".dot");
	}

	void DeleteVioletEdges(Graph& g, const ColorHandler<Graph>& coloring) {
		for (auto it = g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
			if (coloring.Color(*it) == edge_type::violet) {
				g.DeleteEdge(*it);
			}
		}
		Cleaner<Graph>(g).Clean();
	}

	void SaveOldGraph(const string& path) {
		INFO("Saving graph to " << path);
		PrintGraphPack(path, gp_);
//		LengthIdGraphLabeler<Graph> labeler(gp_.g);
//		WriteToDotFile(gp_.g, labeler, path + ".dot");
	}

	template<class gp_t2>
	void ProduceResults(gp_t2& gp, const ColorHandler<Graph>& coloring,
			const string& output_folder, bool detailed_output) {
		INFO("Removing unnecessary edges");
		DeleteVioletEdges(gp.g, coloring);

//		if (detailed_output) {
//			PrintColoredGraph(gp.g, coloring, gp.edge_pos,
//					output_folder + "initial_pics/purple_removed.dot");
//			UniversalSaveGP(gp, output_folder + "saves/purple_removed");
//		}

//		ReliableSplitter<Graph> splitter(gp.g, /*max_size*/100, /*edge_length_bound*/5000);
//		BreakPointsFilter<Graph> filter(gp.g, coloring, 3);
		INFO("Counting stats, outputting pictures");
		BPGraphStatCounter<Graph> counter(gp.g, coloring, output_folder);
		LengthIdGraphLabeler<Graph> labeler(gp.g);
		counter.CountStats(labeler, detailed_output);
	}

	void PrepareDirs(const string& output_folder, bool detailed_output) {
		rm_dir(output_folder);
		make_dir(output_folder);
		if (detailed_output) {
			make_dir(output_folder + "initial_pics/");
			make_dir(output_folder + "saves/");
			make_dir(output_folder + "purple_edges_pics/");
		}
	}

public:

	AssemblyComparer(io::IReader<io::SingleRead> &stream1,
			io::IReader<io::SingleRead> &stream2, const string& name1,
			const string& name2, bool untangle = false,
			const Sequence& reference = Sequence()) :
			gp_(reference, 200, true), coloring_(gp_.g), rc_stream1_(stream1), rc_stream2_(
					stream2), name1_(name1), stream1_(rc_stream1_, name1), name2_(
					name2), stream2_(rc_stream2_, name2), untangle_(untangle) {
	}

	void CompareAssemblies(const string& output_folder, bool detailed_output =
			true, bool one_many_resolve = false, int br_delta = -1,
			const string& add_saves_path = "") {
//		VERIFY(gp_.genome.size() > 0);
		//todo ???
		stream1_.reset();
		stream2_.reset();

		PrepareDirs(output_folder, detailed_output);

//		INFO("Constructing graph");
//		INFO("K = " << gp_t::k_value);
//		ConstructGraph<gp_t::k_value, Graph>(gp_.g, gp_.index, stream1_,
//				stream2_);
//
//		//TODO do we still need it?
//		if (br_delta > 0)
//			SimplifyGraph(gp_.g, (size_t) br_delta);
//
//		ColorHandler<Graph> coloring(gp_.g);
//		ColoredGraphConstructor<Graph, Mapper> colored_graph_constructor(gp_.g,
//				coloring, *MapperInstance < gp_t > (gp_));
//		colored_graph_constructor.ConstructGraph(vector<ContigStream*> {
//				&stream1_, &stream2_ });
//
//		INFO("Filling contig positions");
//		stream1_.reset();
//		FillPos(gp_, stream1_);
//		stream2_.reset();
//		FillPos(gp_, stream2_);

		vector<ContigStream*> streams = { &stream1_, &stream2_ };
		ConstructColoredGraph(gp_, coloring_, streams, br_delta);

		if (gp_.genome.size() > 0) {
			INFO("Filling ref pos " << gp_.genome.size());
//			FillPos(gp_, gp_.genome, "ref_0");
//			FillPos(gp_, !gp_.genome, "ref_1");

//			SimpleInDelAnalyzer<Graph> del_analyzer(
//					gp_.g,
//					coloring_,
//					gp_.edge_pos,
//					(*MapperInstance < gp_t > (gp_)).MapSequence(gp_.genome).simple_path().sequence(),
//					edge_type::red);
//			del_analyzer.Analyze();

//			AlternatingPathsCounter<Graph> alt_count(gp_.g, coloring);
//			alt_count.CountPaths();

			ContigBlockStats<Graph, Mapper> block_stats(gp_.g, gp_.edge_pos,
					*MapperInstance(gp_), gp_.genome, stream1_);
			block_stats.Count();

			MissingGenesAnalyser<Graph, Mapper> missed_genes(gp_.g, coloring_,
					gp_.edge_pos, gp_.genome, *MapperInstance(gp_),
					vector<pair<bool, pair<size_t, size_t>>> {
					make_pair(true, make_pair(260354, 260644)),
					make_pair(true, make_pair(300641, 300904)),
					make_pair(true, make_pair(300904, 301920)),
					make_pair(true, make_pair(301917, 302348)),
					make_pair(true, make_pair(260354, 260644)),
					make_pair(true, make_pair(300641, 300904)),
					make_pair(true, make_pair(300904, 301920)),
					make_pair(true, make_pair(301917, 302348)),
					make_pair(true, make_pair(302449, 304752)),
					make_pair(true, make_pair(263821, 264594)),
					make_pair(true, make_pair(265025, 265726)),
					make_pair(true, make_pair(265740, 266951))
				}
				, output_folder + "missed_genes/");

			missed_genes.Analyze();
			}
//		2339834
//		INFO("Removing gaps");
//		GapsRemover<Graph> gaps_remover(gp_.g, coloring, edge_type::blue, 700);
//		gaps_remover.RemoveGaps();
//		INFO("Gaps removed");

		if (boost::starts_with(name1_, "idba")) {
			IDBADiffAnalyzer<gp_t> diff_analyzer(gp_, coloring_, name1_, name2_,
					output_folder + "/idba_analysis/");
			diff_analyzer.Analyze(stream1_, stream2_);
		}

		if (one_many_resolve) {
			VERIFY(!untangle_);
			RestrictedOneManyResolver<Graph> resolver(gp_.g, coloring_,
					edge_type::violet);
			resolver.Resolve();
		}

		if (detailed_output) {
			if (gp_.genome.size() > 0) {
				PrintColoredGraphAlongRef(gp_, coloring_, gp_.edge_pos,
						gp_.genome,
						output_folder + "initial_pics/colored_split_graph.dot");
			} else {
				PrintColoredGraph(gp_.g, coloring_, gp_.edge_pos,
						output_folder + "initial_pics/colored_split_graph.dot");
			}
//			if (add_saves_path != "") {
//				UniversalSaveGP(gp_, //coloring,
//						add_saves_path);
//				SaveColoring(gp_.g, gp_.int_ids, coloring, add_saves_path);
//				PrintColoredGraph(gp_.g, coloring, gp_.edge_pos,
//						add_saves_path + ".dot");
//			}
			UniversalSaveGP(gp_, //coloring,
					output_folder + "saves/colored_split_graph");
			SaveColoring(gp_.g, gp_.int_ids, coloring_,
					output_folder + "saves/colored_split_graph");
			PrintColoredGraph(gp_.g, coloring_, gp_.edge_pos,
					output_folder + "saves/colored_split_graph.dot");
		}

		if (untangle_) {
			VERIFY(false);
//			INFO("Untangling graph");
//			bp_graph_pack<typename gp_t::graph_t> untangled_gp(gp_t::k_value);
//			UntangledGraphConstructor<gp_t> untangler(gp_, coloring,
//					untangled_gp, stream1_, stream2_);
//			//todo ???
//			//		SimplifyGraph(untangled_gp.g);
//
//			if (detailed_output) {
//				PrintColoredGraph(untangled_gp.g, untangled_gp.coloring,
//						untangled_gp.edge_pos,
//						output_folder + "initial_pics/untangled_graph.dot");
//				UniversalSaveGP(untangled_gp, //untangled_gp.coloring,
//						output_folder + "saves/untangled_graph");
//			}
//
//			ProduceResults(untangled_gp, untangled_gp.coloring, output_folder,
//					detailed_output);
		} else {
//			INFO("Analyzing gaps");
//			GapComparativeAnalyzer<Graph> gap_analyzer(gp_.g, coloring,
//					gp_.edge_pos);
//			gap_analyzer.ReportPotentialGapsCloses(
//					output_folder + "gap_closing_edges/");

			//trivial breakpoints
			string bp_folder = output_folder + "breakpoints/";
			make_dir(bp_folder);
			TrivialBreakpointFinder<Graph> bp_finder(gp_.g, coloring_,
					gp_.edge_pos);
			bp_finder.FindBreakPoints(bp_folder);

			//possible rearrangements
			string rearr_folder = output_folder + "rearrangements/";
			make_dir(rearr_folder);
			SimpleRearrangementDetector<gp_t> rearr_det(gp_, coloring_, "tdc_",
					rearr_folder);
			rearr_det.Detect();

			ProduceResults(gp_, coloring_, output_folder, detailed_output);
		}
	}

private:
	DECL_LOGGER("AssemblyComparer")
	;
};

template<size_t k, size_t K>
void RunBPComparison(ContigStream& raw_stream1, ContigStream& raw_stream2,
		const string& name1, const string& name2, bool refine, bool untangle,
		const string& output_folder, bool detailed_output = true, size_t delta =
				5, Sequence reference = Sequence(),
		const string& add_saves_path = "") {
	io::SplittingWrapper stream1(raw_stream1);
	io::SplittingWrapper stream2(raw_stream2);

	typedef debruijn_graph::graph_pack<
	/*Nonc*/debruijn_graph::ConjugateDeBruijnGraph, K> comparing_gp_t;

	if (refine) {
		typedef graph_pack<ConjugateDeBruijnGraph, k> refining_gp_t;
		refining_gp_t refining_gp;
		io::VectorReader<io::SingleRead> genome_stream(
				io::SingleRead("genome", reference.str()));
		vector<ContigStream*> comp_stream =
				{ &stream1, &stream2, &genome_stream };

		ConstructGPForRefinement(refining_gp, comp_stream, delta);

		io::ModifyingWrapper<io::SingleRead> refined_stream1(
				stream1,
				GraphReadCorrectorInstance(refining_gp.g,
						*MapperInstance(refining_gp)));
		io::ModifyingWrapper<io::SingleRead> refined_stream2(
				stream2,
				GraphReadCorrectorInstance(refining_gp.g,
						*MapperInstance(refining_gp)));
		io::ModifyingWrapper<io::SingleRead> reference_stream(
				genome_stream,
				GraphReadCorrectorInstance(refining_gp.g,
						*MapperInstance(refining_gp)));

		reference_stream.reset();
		AssemblyComparer<comparing_gp_t> comparer(refined_stream1,
				refined_stream2, name1, name2, untangle,
				ReadSequence(reference_stream));
		comparer.CompareAssemblies(output_folder, detailed_output, /*one_many_resolve*/
		false, 10, add_saves_path);
	} else {
		AssemblyComparer<comparing_gp_t> comparer(stream1, stream2, name1,
				name2, untangle, reference);
		comparer.CompareAssemblies(output_folder, detailed_output, /*one_many_resolve*/
		false, 10, add_saves_path);
	}
}

template<size_t k, size_t K>
void RunBPComparison(const Sequence& ref, ContigStream& stream,
		const string& name1, const string& name2, bool refine, bool untangle,
		const string& output_folder, bool detailed_output = true, size_t delta =
				5) {
	io::VectorReader<io::SingleRead> ref_stream(
			io::SingleRead(name1, ref.str()));
	RunBPComparison<k, K>(ref_stream, stream, name1, name2, refine, untangle,
			output_folder, detailed_output, delta);
}

template<size_t k, size_t K>
void RunBPComparison(const Sequence& s1, const Sequence& s2,
		const string& name1, const string& name2, bool refine, bool untangle,
		const string& output_folder, bool detailed_output = true) {
	io::VectorReader<io::SingleRead> stream(io::SingleRead(name2, s2.str()));
	RunBPComparison<k, K>(s1, stream, name1, name2, refine, untangle,
			output_folder, detailed_output);
}

template<size_t k, size_t K>
void RunBPComparison(const Sequence& ref, const vector<Sequence>& contigs,
		const string& name1, const string& name2, bool refine, bool untangle,
		const string& output_folder, bool detailed_output = true) {
	io::VectorReader<io::SingleRead> stream(MakeReads(contigs));
	RunBPComparison<k, K>(ref, stream, name1, name2, refine, untangle,
			output_folder, detailed_output);
}

template<size_t k>
void CompareGenomes(const Sequence& genome_1, const Sequence& genome_2,
		const string& output_dir) {
	INFO("Genome comparison started");
	io::VectorReader<io::SingleRead> stream1(
			io::SingleRead("", genome_1.str()));
	io::VectorReader<io::SingleRead> stream2(
			io::SingleRead("", genome_2.str()));
	typedef graph_pack</*Nonc*/ConjugateDeBruijnGraph, k> comparing_gp_t;
	INFO("Running assembly comparer");
	AssemblyComparer<comparing_gp_t> comparer(stream1, stream2, "strain1",
			"strain2", /*untangle*/false);
	comparer.CompareAssemblies(output_dir, /*detailed_output*/true, /*on_many_resolve*/
	false);
	INFO("Finished");
}

template<class gp_t>
void ThreadAssemblies(const string& base_saves, ContigStream& base_assembly,
		const string& base_prefix, ContigStream& assembly_to_thread,
		const string& to_thread_prefix, const string& output_dir) {
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

}

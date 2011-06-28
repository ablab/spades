/*
 * launch.hpp
 *
 *  Created on: May 6, 2011
 *      Author: sergey
 */

#ifndef LAUNCH_HPP_
#define LAUNCH_HPP_

#include "visualization_utils.hpp"
#include "ireadstream.hpp"

//#include "debruijn_graph.hpp"
#include "paired_info.hpp"
#include "debruijn_graph_constructor.hpp"
#include "tip_clipper.hpp"
#include "bulge_remover.hpp"
#include "erroneous_connection_remover.hpp"
#include "coverage_handler.hpp"
#include "repeat_resolver.hpp"
#include "omni_tools.hpp"
#include "seq_map.hpp"
#include "ID_track_handler.hpp"
#include "read/osequencestream.hpp"
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "new_debruijn.hpp"
#include "config.hpp"

namespace debruijn_graph {

DECL_LOGGER("debruijn_graph")

typedef NewConjugateDeBruijnGraph Graph;
typedef Graph::EdgeId EdgeId;
typedef Graph::VertexId VertexId;

using namespace omnigraph;

template<size_t k>
void CountStats(Graph& g, const EdgeIndex<k + 1, Graph>& index,
		const string& genome) {
	INFO("Counting stats");
	StatCounter<Graph, k> stat(g, index, genome);
	stat.Count();
	INFO("Stats counted");
}

void WriteToDotFile(Graph &g, const string& file_name, string graph_name,
		Path<EdgeId> path1/* = Path<EdgeId> ()*/, Path<EdgeId> path2/* = Path<EdgeId> ()*/) {
	INFO("Writing graph '" << graph_name << "' to file " << file_name);
	gvis::WritePaired(file_name, graph_name, g,	path1, path2);
	INFO("Graph '" << graph_name << "' written to file " << file_name);
}

void DetailedWriteToDot(Graph &g, const string& file_name, string graph_name,
		Path<EdgeId> path1/* = Path<EdgeId> ()*/, Path<EdgeId> path2/* = Path<EdgeId> ()*/) {
	INFO("Writing graph '" << graph_name << "' to file " << file_name);
	gvis::WriteToFile(file_name, graph_name, g,
			path1, path2);
	INFO("Graph '" << graph_name << "' written to file " << file_name);
}

template<size_t k>
Path<typename Graph::EdgeId> FindGenomePath(const string &genome,
		const Graph& g, const EdgeIndex<k + 1, Graph>& index) {
	SimpleSequenceMapper<k, Graph> srt(g, index);
	return srt.MapSequence(Sequence(genome));
}

template<size_t k>
void ProduceInfo(Graph& g, const EdgeIndex<k + 1, Graph>& index,
		const string& genome, const string& file_name, const string& graph_name) {
	CountStats<k> (g, index, genome);
	Path<typename Graph::EdgeId> path1 = FindGenomePath<k> (genome, g, index);
	Path<typename Graph::EdgeId> path2 = FindGenomePath<k> (ReverseComplement(genome), g, index);
	WriteToDotFile(g, file_name, graph_name, path1, path2);
}

template<size_t k>
void ProduceDetailedInfo(Graph& g, const EdgeIndex<k + 1, Graph>& index,
		const string& genome, const string& folder, const string& file_name, const string& graph_name) {
	CountStats<k> (g, index, genome);
	Path<typename Graph::EdgeId> path1 = FindGenomePath<k> (genome, g, index);
	Path<typename Graph::EdgeId> path2 = FindGenomePath<k> (ReverseComplement(genome), g, index);

	mkdir(folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH | S_IWOTH);
	DetailedWriteToDot(g, folder + file_name, graph_name, path1, path2);
}

void ClipTips(Graph &g) {
	INFO("\n-----------------------------------------");
	INFO("Clipping tips");
	TipComparator<Graph> comparator(g);
	size_t max_tip_length = CONFIG.read<size_t> ("tc_max_tip_length");
	size_t max_coverage = CONFIG.read<size_t> ("tc_max_coverage");
	double max_relative_coverage = CONFIG.read<double> (
			"tc_max_relative_coverage");
	TipClipper<Graph, TipComparator<Graph>> tc(g, comparator, max_tip_length,
			max_coverage, max_relative_coverage);
	tc.ClipTips();
	INFO("Clipping tips finished");
}

void RemoveBulges(Graph &g) {
	INFO("\n-----------------------------------------");
	INFO("Removing bulges");
	double max_coverage = CONFIG.read<double> ("br_max_coverage");
	double max_relative_coverage = CONFIG.read<double> (
			"br_max_relative_coverage");
	double max_delta = CONFIG.read<double> ("br_max_delta");
	double max_relative_delta = CONFIG.read<double> ("br_max_relative_delta");
	size_t max_length_div_K = CONFIG.read<int> ("br_max_length_div_K");
	BulgeRemover<Graph> bulge_remover(max_length_div_K * g.k(), max_coverage,
			max_relative_coverage, max_delta, max_relative_delta);
	bulge_remover.RemoveBulges(g);
	INFO("Bulges removed");
}

void RemoveLowCoverageEdges(Graph &g) {
	INFO("\n-----------------------------------------");
	INFO("Removing low coverage edges");
	double max_coverage = CONFIG.read<double> ("ec_max_coverage");
	int max_length_div_K = CONFIG.read<int> ("ec_max_length_div_K");
	LowCoverageEdgeRemover<Graph> erroneous_edge_remover(
			max_length_div_K * g.k(), max_coverage);
	erroneous_edge_remover.RemoveEdges(g);
	INFO("Low coverage edges removed");
}

void ResolveRepeats(Graph &g, IdTrackHandler<Graph> &old_IDs, PairedInfoIndex<Graph> &info, Graph &new_graph, IdTrackHandler<Graph> &new_IDs, const string& output_folder) {
	INFO("\n-----------------------------------------");
	INFO("Resolving primitive repeats");
	RepeatResolver<Graph> repeat_resolver(g, old_IDs, 0, info, new_graph, new_IDs);
	mkdir((output_folder).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH| S_IWOTH);
	repeat_resolver.ResolveRepeats(output_folder);
	INFO("Primitive repeats resolved");
}

template<size_t k, class ReadStream>
void FillPairedIndex(Graph &g, PairedInfoIndex<Graph>& paired_info_index,
		ReadStream& stream, EdgeIndex<k + 1, Graph>& index) {
	INFO("\n-----------------------------------------");
	stream.reset();
	INFO("Counting paired info");
	PairedIndexFiller<Graph, k, ReadStream> pif(g, index, stream);
	pif.FillIndex(paired_info_index);
	INFO("Paired info counted");
}

template<size_t k, class ReadStream>
void FillCoverage(Graph& g/*CoverageHandler<Graph> coverage_handler*/, ReadStream& stream,
		EdgeIndex<k + 1, Graph>& index) {
	typedef SimpleSequenceMapper<k, Graph> ReadThreader;
	INFO("\n-----------------------------------------");
	stream.reset();
	INFO("Counting coverage");
	ReadThreader read_threader(g, index);
	//todo temporary solution!
	g.FillCoverage<ReadStream, ReadThreader>(stream, read_threader);
//	coverage_handler.FillCoverage<k, ReadStream> (stream, index);
	INFO("Coverage counted");
}

template<size_t k, class ReadStream>
void ConstructGraph(Graph& g, EdgeIndex<k + 1, Graph>& index,
		ReadStream& stream) {
	typedef SeqMap<k + 1, typename Graph::EdgeId> DeBruijn;
	INFO("\n-----------------------------------------");
	INFO("Constructing DeBruijn graph");
	DeBruijn& debruijn = index.inner_index();
	INFO("Filling DeBruijn graph");
	debruijn.Fill(stream);
	INFO("DeBruijn graph constructed");

	INFO("Condensing graph");
	DeBruijnGraphConstructor<k, Graph> g_c(debruijn);
	g_c.ConstructGraph(g, index);
	INFO("Graph condensed");
}

template<size_t k, class ReadStream>
void ConstructGraphWithCoverage(Graph& g, EdgeIndex<k + 1, Graph>& index,
		/*CoverageHandler<Graph>& coverage_handler,*/ ReadStream& stream) {
	ConstructGraph<k, ReadStream> (g, index, stream);
	FillCoverage<k, ReadStream> (g/*coverage_handler*/, stream, index);
}

template<size_t k, class PairedReadStream>
void ConstructGraphWithPairedInfo(Graph& g, EdgeIndex<k + 1, Graph>& index,
		/*CoverageHandler<Graph>& coverage_handler,*/
		PairedInfoIndex<Graph>& paired_index, PairedReadStream& stream) {
	typedef SimpleReaderWrapper<PairedReadStream> UnitedStream;
	UnitedStream united_stream(stream);
	ConstructGraphWithCoverage<k, UnitedStream> (g, index/*, coverage_handler*/,
			united_stream);
	FillPairedIndex<k, PairedReadStream> (g, paired_index, stream, index);
}

template<size_t k>
void SimplifyGraph(Graph& g, EdgeIndex<k + 1, Graph>& index, size_t iteration_count, const string& genome, const string& output_folder) {
	INFO("\n-----------------------------------------");
	INFO("Graph simplification started");

	ProduceDetailedInfo<k> (g, index, genome,
			output_folder + "before_simplification/", "graph.dot",
			"non_simplified_graph");
	for (size_t i = 0; i < iteration_count; i++) {
		INFO("\n-----------------------------------------");
		INFO("Iteration " << i);

		ClipTips(g);
		ProduceDetailedInfo<k> (g, index, genome,
				output_folder + "tips_clipped_" + ToString(i) + "/", "graph.dot",
				"no_tip_graph");

		RemoveBulges(g);
		ProduceDetailedInfo<k> (g, index, genome,
				output_folder + "bulges_removed_" + ToString(i) + "/", "graph.dot",
				"no_bulge_graph");

		RemoveLowCoverageEdges(g);
		ProduceDetailedInfo<k> (
				g,
				index,
				genome,
				output_folder + "erroneous_edges_removed_" + ToString(i) + "/",
						"graph.dot", "no_erroneous_edges_graph");
	}
	INFO("Graph simplification finished");
}

void OutputContigs(Graph& g, const string& contigs_output_filename) {
	INFO("\n-----------------------------------------");
	INFO("Outputting contigs to " << contigs_output_filename);

	osequencestream oss(contigs_output_filename);
	//TipComparator<Graph> compare(g); // wtf, don't we have usual less for edges?
	for (auto it = g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
		oss << g.EdgeNucls(*it);
	}
	INFO("Contigs written");
}

template<size_t k, class ReadStream>
void DeBruijnGraphWithPairedInfoTool(ReadStream& stream, const string& genome,
		const string& output_folder, const string& work_tmp_dir) {
	INFO("Edge graph construction tool started");

	Graph g(k);
	EdgeIndex<k + 1, Graph> index(g);
	PairedInfoIndex<Graph> paired_index(g);
	IdTrackHandler<Graph> IntIds(g);

	ConstructGraphWithPairedInfo<k, ReadStream> (g, index, paired_index, stream);

	ProduceInfo<k> (g, index, genome, output_folder + "edge_graph.dot",
			"edge_graph");

	SimplifyGraph<k>(g, index, 3, genome, output_folder);

	ProduceInfo<k> (g, index, genome, output_folder + "simplified_graph.dot",
			"simplified_graph");

//	paired_index.OutputData(output_folder + "edges_dist.txt");

//	SimpleOfflineClusterer<Graph> clusterer(paired_index);
//	PairedInfoIndex<Graph> clustered_paired_index(g);
//	clusterer.cluster(clustered_paired_index);

	INFO("before ResolveRepeats");
	RealIdGraphLabeler<Graph> IdTrackLabelerBefore(g, IntIds);
	gvis::WriteSimple( output_folder + "repeats_resolved_simple_before.dot", "no_repeat_graph", g, IdTrackLabelerBefore);

	Graph new_graph(k);
	IdTrackHandler<Graph> NewIntIds(new_graph, IntIds.MaxVertexId(), IntIds.MaxEdgeId());
	ResolveRepeats(g, IntIds, paired_index, new_graph, NewIntIds, output_folder+"resolve/");
	INFO("before graph writing");
	RealIdGraphLabeler<Graph> IdTrackLabelerAfter(new_graph, NewIntIds);
	gvis::WriteSimple( output_folder + "repeats_resolved_simple_after.dot", "no_repeat_graph", new_graph, IdTrackLabelerAfter);
		INFO("repeat resolved grpah written");
	ProduceInfo<k> (new_graph, index, genome, output_folder + "repeats_resolved.dot",
			"no_repeat_graph");

	ProduceInfo<k> (new_graph, index, genome, work_tmp_dir + "repeats_resolved.dot",
			"no_repeat_graph");
	OutputContigs(g, output_folder + "contigs.fasta");
	INFO("Tool finished");

}

template<size_t k, class ReadStream>
void DeBruijnGraphTool(ReadStream& stream, const string& genome,
		const string& output_folder) {
	INFO("Edge graph construction tool started");

	Graph g(k);
	EdgeIndex<k + 1, Graph> index(g);
	IdTrackHandler<Graph> IntIds(g);

	typedef SimpleReaderWrapper<ReadStream> UnitedStream;
	UnitedStream united_stream(stream);
	ConstructGraphWithCoverage<k, UnitedStream> (g, index, united_stream);

	ProduceInfo<k> (g, index, genome, output_folder + "edge_graph.dot",
			"edge_graph");

	SimplifyGraph<k>(g, index, 3, genome, output_folder);

	ProduceInfo<k> (g, index, genome, output_folder + "simplified_graph.dot",
			"simplified_graph");

	OutputContigs(g, output_folder + "contigs.fasta");
	INFO("Tool finished");
}

}

#endif /* LAUNCH_HPP_ */

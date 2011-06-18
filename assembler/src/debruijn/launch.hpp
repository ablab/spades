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

#include "debruijn_graph.hpp"
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


namespace debruijn_graph {

DECL_LOGGER("debruijn_graph")

typedef DeBruijnGraph Graph;
typedef Graph::EdgeId EdgeId;
typedef Graph::VertexId VertexId;

using namespace omnigraph;

template<size_t k>
void CountStats(Graph& g, const EdgeIndex<k + 1, Graph>& index,
		const string& genome) {
	INFO("Counting stats");
	StatCounter<Graph, k> stat(g, index, genome);
	stat.CountStatistics();
	INFO("Stats counted");
}

void WriteToDotFile(Graph &g, const string& file_name, string graph_name,
		Path<EdgeId> path = Path<EdgeId> ()) {
	INFO("Writing graph '" << graph_name << "' to file " << file_name);
	gvis::WriteToFile(/*DE_BRUIJN_DATA_FOLDER + */file_name, graph_name, g,
			path);
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
	Path<typename Graph::EdgeId> path = FindGenomePath<k> (genome, g, index);
	WriteToDotFile(g, file_name, graph_name, path);
}

void ClipTips(Graph &g) {
	INFO("Clipping tips");
	TipComparator<Graph> comparator(g);
	size_t max_tip_length = CONFIG.read<size_t> ("tc_max_tip_length");
	size_t max_coverage = CONFIG.read<size_t> ("tc_max_coverage");
	double max_relative_coverage = CONFIG.read<double> (
			"tc_max_relative_coverage");
	TipClipper<Graph, TipComparator<Graph>> tc(g, comparator, max_tip_length,
			max_coverage, max_relative_coverage);
	tc.ClipTips();
}

void RemoveBulges(Graph &g) {
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
	INFO("Removing low coverage edges");
	double max_coverage = CONFIG.read<double> ("ec_max_coverage");
	int max_length_div_K = CONFIG.read<int> ("ec_max_length_div_K");
	LowCoverageEdgeRemover<Graph> erroneous_edge_remover(
			max_length_div_K * g.k(), max_coverage);
	erroneous_edge_remover.RemoveEdges(g);
	INFO("Low coverage edges removed");
}

void ResolveRepeats(Graph &g, PairedInfoIndex<Graph> &info, Graph &new_graph) {
	INFO("Resolving primitive repeats");
	RepeatResolver<Graph> repeat_resolver(g, 0, info, new_graph);
	repeat_resolver.ResolveRepeats();
	INFO("Primitive repeats resolved");
}

template<size_t k, class ReadStream>
void FillPairedIndex(PairedInfoIndex<Graph>& paired_info_index,
		ReadStream& stream, EdgeIndex<k + 1, Graph>& index) {
	stream.reset();
	INFO("Counting paired info");
	paired_info_index.FillIndex<k, ReadStream> (index, stream);
	INFO("Paired info counted");
}

template<size_t k, class ReadStream>
void FillCoverage(CoverageHandler<Graph> coverage_handler, ReadStream& stream,
		EdgeIndex<k + 1, Graph>& index) {
	stream.reset();
	INFO("Counting coverage");
	coverage_handler.FillCoverage<k, ReadStream> (stream, index);
	INFO("Coverage counted");
}

template<size_t k, class ReadStream>
void ConstructGraph(Graph& g, EdgeIndex<k + 1, Graph>& index,
		ReadStream& stream) {
	typedef SeqMap<k + 1, typename Graph::EdgeId> DeBruijn;

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
		CoverageHandler<Graph>& coverage_handler, ReadStream& stream) {
	ConstructGraph<k, ReadStream> (g, index, stream);
	FillCoverage<k, ReadStream> (coverage_handler, stream, index);
}

template<size_t k, class PairedReadStream>
void ConstructGraphWithPairedInfo(Graph& g, EdgeIndex<k + 1, Graph>& index,
		CoverageHandler<Graph>& coverage_handler,
		PairedInfoIndex<Graph>& paired_index, PairedReadStream& stream) {
	typedef SimpleReaderWrapper<PairedReadStream> UnitedStream;
	UnitedStream united_stream(stream);
	ConstructGraphWithCoverage<k, UnitedStream> (g, index, coverage_handler,
			united_stream);
	FillPairedIndex<k, PairedReadStream> (paired_index, stream, index);
}

template<size_t k, class ReadStream>
void DeBruijnGraphTool(ReadStream& stream, const string& genome,
		const string& output_folder) {
	INFO("Edge graph construction tool started");

	Graph g(k);
	EdgeIndex<k + 1, Graph> index(g);
	CoverageHandler<Graph> coverage_handler(g);
	PairedInfoIndex<Graph> paired_index(g);
	IdTrackHandler<Graph> IntIds(g);

	ConstructGraphWithPairedInfo<k, ReadStream> (g, index, coverage_handler,
			paired_index, stream);
	//	{
	//		typedef SimpleReaderWrapper<ReadStream> UnitedStream;
	//		UnitedStream united_stream(stream);
	//		ConstructGraphWithCoverage<k, UnitedStream> (g, index, coverage_handler, united_stream);
	//	}

	ProduceInfo<k> (g, index, genome, output_folder + "edge_graph.dot",
			"edge_graph");
	//	paired_index.OutputData(output_folder + "edges_dist.txt");

	for (size_t i = 0; i < 3; i++) {
		ClipTips(g);
		ProduceInfo<k> (g, index, genome,
				output_folder + "tips_clipped_" + ToString(i) + ".dot",
				"no_tip_graph");

		RemoveBulges(g);
		ProduceInfo<k> (g, index, genome,
				output_folder + "bulges_removed_" + ToString(i) + ".dot",
				"no_bulge_graph");

		RemoveLowCoverageEdges(g);
		ProduceInfo<k> (
				g,
				index,
				genome,
				output_folder + "erroneous_edges_removed_" + ToString(i)
						+ ".dot", "no_erroneous_edges_graph");
	}

//	SimpleOfflineClusterer<Graph> clusterer(paired_index);
//	PairedInfoIndex<Graph> clustered_paired_index(g);
//	clusterer.cluster(clustered_paired_index);
	INFO("before ResolveRepeats");
	Graph new_graph(k);
	ResolveRepeats(g, paired_index, new_graph);
	INFO("before graph writing");
	gvis::WriteSimple( output_folder + "repeats_resolved_siiimple.dot", "no_repeat_graph", new_graph);
	INFO("repeat resolved grpah written");
	ProduceInfo<k> (new_graph, index, genome, output_folder + "repeats_resolved.dot",
			"no_repeat_graph");
	INFO("Tool finished");

}

}

#endif /* LAUNCH_HPP_ */

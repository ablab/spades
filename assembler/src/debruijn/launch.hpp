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

#include "paired_info.hpp"
#include "edge_graph_constructor.hpp"
#include "tip_clipper.hpp"
#include "bulge_remover.hpp"
#include "coverage_handler.hpp"

namespace edge_graph {

using de_bruijn::EdgeIndex;
using de_bruijn::DeBruijnPlus;
using de_bruijn::CoverageHandler;
typedef de_bruijn::Path<EdgeId> Path;
typedef de_bruijn::PairedInfoIndex<EdgeGraph> PairedIndex;

void CountStats(const EdgeGraph& g) {
	INFO("Counting stats");
	de_bruijn::DFS<EdgeGraph> dfs(g);
	de_bruijn::SimpleStatCounter<EdgeGraph> stat_c;
	dfs.Traverse(&stat_c);
	INFO("Vertex count=" << stat_c.v_count() << "; Edge count="
			<< stat_c.e_count());
	INFO("Stats counted");
}

void WriteToDotFile(const EdgeGraph &g, const string& file_name,
		string graph_name,
		de_bruijn::Path<EdgeId> path = de_bruijn::Path<EdgeId>()) {
	INFO("Writing graph '" << graph_name << "' to file " << file_name);
	WriteToFile(/*DE_BRUIJN_DATA_FOLDER + */file_name, graph_name, g, path);
	INFO("Graph '" << graph_name << "' written to file " << file_name);
}

template<size_t k>
Path FindGenomePath(const string &genome, const EdgeGraph& g,
		const EdgeIndex<k + 1, EdgeGraph>& index) {
	de_bruijn::SimpleSequenceMapper<k, EdgeGraph> srt(g, index);
	return srt.MapSequence(Sequence(genome));
}

template<size_t k>
void ProduceInfo(const EdgeGraph& g, const EdgeIndex<k + 1, EdgeGraph>& index,
		const string& genome, const string& file_name, const string& graph_name) {
	CountStats(g);
	Path path = FindGenomePath<k> (genome, g, index);
	WriteToDotFile(g, file_name, graph_name, path);
}

void ClipTips(EdgeGraph &g) {
	INFO("Clipping tips");
	TipComparator<EdgeGraph> comparator(g);
	TipClipper<EdgeGraph, TipComparator<EdgeGraph>> tc(comparator);
	tc.ClipTips(g);
}

void RemoveBulges(EdgeGraph &g) {
	INFO("Removing bulges");
	de_bruijn::BulgeRemover<EdgeGraph> bulge_remover(5 * g.k());
	bulge_remover.RemoveBulges(g);
	INFO("Bulges removed");
}

template<size_t k, class ReadStream>
void FillPairedIndex(PairedIndex& paired_info_index, ReadStream& stream,
		EdgeIndex<k + 1, EdgeGraph>& index) {
	stream.reset();
	INFO("Counting paired info");
	paired_info_index.FillIndex<k, ReadStream> (index, stream);
	INFO("Paired info counted");
}

template<size_t k, class ReadStream>
void FillCoverage(de_bruijn::CoverageHandler<EdgeGraph> coverage_handler, ReadStream& stream,
		EdgeIndex<k + 1, EdgeGraph>& index) {
	stream.reset();
	INFO("Counting coverage");
	coverage_handler.FillCoverage<k, ReadStream> (stream, index);
	INFO("Coverage counted");
}

template <size_t k, class ReadStream>
void ConstructGraph(EdgeGraph& g, EdgeIndex<k + 1, EdgeGraph>& index, CoverageHandler<EdgeGraph>& coverage_handler, ReadStream& stream) {
	typedef de_bruijn::DeBruijnPlus<k + 1, EdgeId> DeBruijn;

	INFO("Constructing DeBruijn graph");
	DeBruijn debruijn(stream);
	INFO("DeBruijn graph constructed");

	INFO("Condensing graph");
	EdgeGraphConstructor<k> g_c(debruijn);
	g_c.ConstructGraph(g, index);
	INFO("Graph condensed");

	FillCoverage<k, ReadStream> (coverage_handler, stream, index);
}

template <size_t k, class PairedReadStream>
void ConstructGraphWithPairedInfo(EdgeGraph& g, EdgeIndex<k + 1, EdgeGraph>& index
		, CoverageHandler<EdgeGraph>& coverage_handler, PairedIndex& paired_index
		, PairedReadStream& stream) {
	typedef SimpleReaderWrapper<PairedReadStream> UnitedStream;
	UnitedStream united_stream(stream);
	ConstructGraph<k, UnitedStream>(g, index, coverage_handler, united_stream);
	FillPairedIndex<k, PairedReadStream> (paired_index, stream, index);
}

template<size_t k, class ReadStream>
void EdgeGraphTool(ReadStream& stream, const string& genome, const string& output_folder) {
	typedef de_bruijn::DeBruijnPlus<k + 1, EdgeId> DeBruijn;
	INFO("Edge graph construction tool started");

	EdgeGraph g(k);
	EdgeIndex<k + 1, EdgeGraph> index(g);
	de_bruijn::CoverageHandler<EdgeGraph> coverage_handler(g);
	PairedIndex paired_index(g);

	ConstructGraphWithPairedInfo<k, ReadStream>(g, index, coverage_handler, paired_index, stream);

	ProduceInfo<k> (g, index, genome, output_folder + "edge_graph.dot", "edge_graph");
	paired_index.OutputData(output_folder + "edges_dist.txt");

	ClipTips(g);
	ProduceInfo<k> (g, index, genome, output_folder + "tips_clipped.dot", "no_tip_graph");

	RemoveBulges(g);
	ProduceInfo<k> (g, index, genome, output_folder + "bulges_removed.dot", "no_bulge_graph");

	INFO("Tool finished")
}

}

#endif /* LAUNCH_HPP_ */

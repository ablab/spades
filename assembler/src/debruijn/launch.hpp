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
//typedef de_bruijn::DeBruijnPlus<K + 1, EdgeId> DeBruijn;
//typedef de_bruijn::EdgeIndex<K + 1, EdgeGraph> Index;

using de_bruijn::EdgeIndex;
using de_bruijn::DeBruijnPlus;
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

/*template<class ReadStream>
 void ConstructUncondensedGraph(DeBruijn& debruijn, ReadStream& stream) {
 INFO("Constructing DeBruijn graph");
 debruijn.ConstructGraphFromStream(stream);
 INFO("DeBruijn graph constructed");
 }*/

template<size_t k>
void CondenseGraph(DeBruijnPlus<k + 1, EdgeId>& debruijn, EdgeGraph& g,
		EdgeIndex<k + 1, EdgeGraph>& index) {
	INFO("Condensing graph");
	EdgeGraphConstructor<k> g_c(debruijn);
	g_c.ConstructGraph(g, index);
	INFO("Graph condensed");
}

void ClipTips(EdgeGraph &g) {
	INFO("Clipping tips");
	TipComparator<EdgeGraph> comparator(g);
	TipClipper<EdgeGraph, TipComparator<EdgeGraph>> tc(comparator);
	tc.ClipTips(g);
	INFO("Tips clipped");
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

template<size_t k, class ReadStream>
void EdgeGraphTool(ReadStream& stream, size_t insert_size, const string& genome) {
	typedef de_bruijn::DeBruijnPlus<k + 1, EdgeId> DeBruijn;
	INFO("Edge graph construction tool started");

	INFO("Constructing DeBruijn graph");
	typedef SimpleReaderWrapper<ReadStream> UnitedStream;
	UnitedStream unitedStream(stream);
	DeBruijn debruijn(unitedStream);
	INFO("DeBruijn graph constructed");

	EdgeGraph g(k);
	EdgeIndex<k + 1, EdgeGraph> index(g, debruijn);
	CondenseGraph<k> (debruijn, g, index);

	de_bruijn::CoverageHandler<EdgeGraph> coverage_handler(g);
	FillCoverage<k, UnitedStream> (coverage_handler, unitedStream, index);
	ProduceInfo<k> (g, index, genome, "edge_graph.dot", "edge_graph");

	PairedIndex paired_index(g, insert_size);
	FillPairedIndex<k, ReadStream> (paired_index, stream, index);

	ClipTips(g);
	ProduceInfo<k> (g, index, genome, "tips_clipped.dot", "no_tip_graph");

	RemoveBulges(g);
	ProduceInfo<k> (g, index, genome, "bulges_removed.dot", "no_bulge_graph");

	INFO("Tool finished")
}

}

#endif /* LAUNCH_HPP_ */

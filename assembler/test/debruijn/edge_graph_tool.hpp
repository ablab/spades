/*
 * edge_graph_tool.hpp
 *
 *  Created on: Apr 10, 2011
 *      Author: sergey
 */

#ifndef EDGE_GRAPH_TOOL_HPP_
#define EDGE_GRAPH_TOOL_HPP_
#include "tip_clipper.hpp"
#include "bulge_remover.hpp"
#include "coverage_counter.hpp"
#include "visualization_utils.hpp"
#include "paired_info.hpp"

namespace edge_graph {

typedef de_bruijn::DeBruijn<K> DeBruijn;
typedef SimpleIndex<K + 1, EdgeId> Index;
typedef de_bruijn::Path<EdgeId> Path;

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
	WriteToFile(DE_BRUIJN_DATA_FOLDER + file_name, graph_name, g, path);
	INFO("Graphraph '" << graph_name << "' written to file " << file_name);
}

Path FindGenomePath(const string &genome, const EdgeGraph& g,
		const Index& index) {
	de_bruijn::SimpleReadThreader<K, EdgeGraph> srt(g, index);
	return srt.ThreadRead(Sequence(genome));
}

void ProduceInfo(const EdgeGraph& g, const Index& index, const string& genome,
		const string& file_name, const string& graph_name) {
	CountStats(g);
	Path path = FindGenomePath(genome, g, index);
	WriteToDotFile(g, file_name, graph_name, path);
}

template<class ReadStream>
void ConstructUncondensedGraph(DeBruijn& debruijn, ReadStream& stream) {
	INFO("Constructing DeBruijn graph");
	debruijn.ConstructGraphFromStream(stream);
	INFO("DeBruijn graph constructed");
}

template<class ReadStream>
void CondenseGraph(const DeBruijn& debruijn, EdgeGraph& g, Index& index,
		ReadStream& stream, const string& genome) {
	INFO("Condensing graph");
	CondenseConstructor<K> g_c(debruijn);
	g_c.ConstructGraph(g, index);
	INFO("Graph condensed");

	INFO("Counting coverage");
	CoverageCounter<K, EdgeGraph> cc(g, index);
	cc.CountCoverage(stream);
	INFO("Coverage counted");

	ProduceInfo(g, index, genome, "edge_graph.dot", "edge_graph");
}

void ClipTips(EdgeGraph &g, Index &index, const string& genome,
		string dotFileName) {
	INFO("Clipping tips");
	TipComparator<EdgeGraph> comparator(g);
	TipClipper<EdgeGraph, TipComparator<EdgeGraph> > tc(comparator);
	tc.ClipTips(g);
	INFO("Tips clipped");

	ProduceInfo(g, index, genome, dotFileName, "no_tip_graph");
}

void RemoveBulges(EdgeGraph &g, Index &index, const string& genome,
		string dotFileName) {
	INFO("Removing bulges");
	de_bruijn::BulgeRemover<EdgeGraph> bulge_remover;
	bulge_remover.RemoveBulges(g);
	INFO("Bulges removed");

	ProduceInfo(g, index, genome, dotFileName, "no_bulge_graph");
}

void EdgeGraphTool(StrobeReader<2, Read, ireadstream>& reader, const string& genome) {
	INFO("Edge graph construction tool started");

	DeBruijn debruijn;

	SimpleReaderWrapper<2, Read, ireadstream> stream(reader);
	ConstructUncondensedGraph<SimpleReaderWrapper<2, Read, ireadstream> > (debruijn, stream);

	//	debruijn.show(genome);

	EdgeGraph g(K);
	SimpleIndex < K + 1, EdgeId > index;
	EdgeHashRenewer<K + 1, EdgeGraph> index_handler(g, index);
	g.AddActionHandler(&index_handler);

	stream.reset();
	CondenseGraph<SimpleReaderWrapper<2, Read, ireadstream> > (debruijn, g, index, stream, genome);

	reader.reset();
	de_bruijn::PairedInfoIndex<EdgeGraph> pairedInfoIndex(g, index, reader);

	ClipTips(g, index, genome, "tips_clipped.dot");

	RemoveBulges(g, index, genome, "bulges_removed.dot");

	g.RemoveActionHandler(&index_handler);

	INFO("Tool finished")
}

}

#endif /* EDGE_GRAPH_TOOL_HPP_ */

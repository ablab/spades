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

namespace edge_graph {

typedef de_bruijn::DeBruijn<K> DeBruijn;
typedef SimpleIndex<K + 1, EdgeId> Index;

void CountStats(const EdgeGraph& g) {
	INFO("Counting stats");
	de_bruijn::DFS<EdgeGraph> dfs(g);
	de_bruijn::SimpleStatCounter<EdgeGraph> stat_c;
	dfs.Traverse(&stat_c);
	INFO("Vertex count=" << stat_c.v_count() << "; Edge count="
			<< stat_c.e_count());
	INFO("Stats counted");
}

void WriteToDotFile(EdgeGraph &g, const string& file_name,
		const string& graph_name,
		de_bruijn::Path<EdgeId> path = de_bruijn::Path<EdgeId>()) {
	INFO("Writing to file");
	WriteToFile(file_name, graph_name, g, path);
}

template<class ReadStream>
void ConstructUncondensedGraph(DeBruijn& debruijn, ReadStream& stream) {
	INFO("Constructing DeBruijn graph");
	debruijn.ConstructGraphFromStream(stream);
	INFO("DeBruijn graph constructed");
}

const de_bruijn::Path<EdgeId> findGenomePath(const string &genome,
		const EdgeGraph& g, Index &index) {
	de_bruijn::SimpleReadThreader<K, EdgeGraph> srt(g, index);
	return srt.ThreadRead(Sequence(genome));
}

template<class ReadStream>
void CondenseGraph(DeBruijn& debruijn, EdgeGraph& g, Index& index,
		ReadStream& stream, string genome = "") {
	INFO("Condensing graph");
	CondenseConstructor<K> g_c(debruijn);
	g_c.ConstructGraph(g, index);
	INFO("Graph condensed");

	INFO("Counting coverage");
	CoverageCounter<K, EdgeGraph> cc(g, index);
	cc.CountCoverage(stream);
	INFO("Coverage counted");

	CountStats(g);
	de_bruijn::Path<EdgeId> path = findGenomePath(genome, g, index);
	WriteToDotFile(g, "edge_graph.dot", "edge_graph", path);
}

void ClipTips(EdgeGraph &g, Index &index, string genome = "") {
	INFO("Clipping tips");
	TipComparator<EdgeGraph> comparator(g);
	TipClipper<EdgeGraph, TipComparator<EdgeGraph> > tc(comparator);
//	cout << "oppa" << endl;
	tc.ClipTips(g);
//	cout << "oppa" << endl;
	INFO("Tips clipped");

	CountStats(g);
	de_bruijn::Path<EdgeId> path = findGenomePath(genome, g, index);
	WriteToDotFile(g, "tips_clipped.dot", "no_tip_graph", path);
}

void RemoveBulges(EdgeGraph &g, Index &index, string genome = "") {
	INFO("Removing bulges");
	de_bruijn::BulgeRemover<EdgeGraph> bulge_remover;
	bulge_remover.RemoveBulges(g);
	INFO("Bulges removed");

	CountStats(g);
	de_bruijn::Path<EdgeId> path = findGenomePath(genome, g, index);
	WriteToDotFile(g, "bulges_removed.dot", "no_bulge_graph", path);
}

template<class ReadStream>
void EdgeGraphTool(ReadStream& stream, string genome = "") {
	INFO("Edge graph construction tool started");

	DeBruijn debruijn;

	ConstructUncondensedGraph<ReadStream> (debruijn, stream);

	EdgeGraph *g = new EdgeGraph(K);
	SimpleIndex<K + 1, EdgeId> *index = new SimpleIndex<K + 1, EdgeId>();
	EdgeHashRenewer<K + 1, EdgeGraph> *indexHandler = new EdgeHashRenewer<K + 1, EdgeGraph>(*g, *index);
	g->AddActionHandler(indexHandler);

	stream.reset();
	CondenseGraph<ReadStream> (debruijn, *g, *index, stream, genome);

	ClipTips(*g, *index, genome);

	RemoveBulges(*g, *index, genome);

	g->RemoveActionHandler(indexHandler);
	delete indexHandler;
	delete g;
	delete index;
	INFO("Tool finished")
}

}

#endif /* EDGE_GRAPH_TOOL_HPP_ */

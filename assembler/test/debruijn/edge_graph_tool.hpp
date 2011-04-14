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

namespace edge_graph {

typedef de_bruijn::DeBruijn<K> DeBruijn;
typedef GraphConstructor<K>::Index Index;

void CountStats(const EdgeGraph& g) {
	INFO("Counting stats");
	de_bruijn::DFS<EdgeGraph> dfs(g);
	de_bruijn::SimpleStatCounter<EdgeGraph> stat_c;
	dfs.Traverse(&stat_c);
	INFO("Vertex count=" << stat_c.v_count() << "; Edge count="
			<< stat_c.e_count());
	INFO("Stats counted");
}

void WriteToDotFile(EdgeGraph* g, const string& file_name, const string& graph_name) {
	INFO("Writing to file");
	WriteToFile(file_name, graph_name, *g);
}

template <class ReadStream>
void ConstructUncondensedGraph(DeBruijn& debruijn, ReadStream& stream) {
	INFO("Constructing DeBruijn graph");
	debruijn.ConstructGraphFromStream(stream);
	INFO("DeBruijn graph constructed");
}

template <class ReadStream>
void CondenseGraph(DeBruijn& debruijn, EdgeGraph*& g, Index*& index, ReadStream& stream) {
	INFO("Condensing graph");
	CondenseConstructor<K> g_c(debruijn);
	g_c.ConstructGraph(g, index);
	INFO("Graph condensed");

	INFO("Counting coverage");
	CoverageCounter<K, EdgeGraph> cc(*g, *index);
	cc.CountCoverage(stream);
	INFO("Coverage counted");

	CountStats(*g);
	WriteToDotFile(g, "edge_graph.dot", "edge_graph");
}

void ClipTips(EdgeGraph* g) {
	INFO("Clipping tips");
	TipComparator<EdgeGraph> comparator(*g);
	TipClipper<EdgeGraph, TipComparator<EdgeGraph> > tc(comparator);
	tc.ClipTips(*g);
	INFO("Tips clipped");

	CountStats(*g);
	WriteToDotFile(g, "tips_clipped.dot", "no_tip_graph");
}

void RemoveBulges(EdgeGraph* g) {
	INFO("Removing bulges");
	de_bruijn::BulgeRemover<EdgeGraph> bulge_remover;
	bulge_remover.RemoveBulges(*g);
	INFO("Bulges removed");

	CountStats(*g);
	WriteToDotFile(g, "bulges_removed.dot", "no_bulge_graph");
}

template <class ReadStream>
void EdgeGraphTool(ReadStream& stream) {
	INFO("Edge graph construction tool started");

	DeBruijn debruijn;

	ConstructUncondensedGraph<ReadStream>(debruijn, stream);

	EdgeGraph *g;
	GraphConstructor<K>::Index *index;

	stream.reset();
	CondenseGraph<ReadStream>(debruijn, g, index, stream);

	ClipTips(g);

	RemoveBulges(g);

	delete g;
	delete index;
	INFO("Tool finished")
}

}


#endif /* EDGE_GRAPH_TOOL_HPP_ */

/*
 * edge_graph_tool.hpp
 *
 *  Created on: Apr 10, 2011
 *      Author: sergey
 */

#ifndef EDGE_GRAPH_TOOL_HPP_
#define EDGE_GRAPH_TOOL_HPP_

#include "read_generator.hpp"
#include "condensed_graph_constructor.hpp"
#include "ireadstream.hpp"
#include "test_utils.hpp"
#include "tip_clipper.hpp"
#include <algorithm>

//#define SUBSTR_LENGTH 1000
//#define COVERAGE 1
//#define R 35
#define K 27

using namespace std;

namespace edge_graph {

void CountStats(const EdgeGraph& g) {
	de_bruijn::DFS<EdgeGraph> dfs(g);
	de_bruijn::SimpleStatCounter<EdgeGraph> stat_c;
	dfs.Traverse(&stat_c);
	INFO("Vertex count=" << stat_c.v_count() << "; Edge count="
			<< stat_c.e_count());
}

void CondenseTool(DeBruijn<K>& debruijn) {
	INFO("Condensing graph");
	CondenseConstructor<K> g_c(debruijn);

	EdgeGraph *g;
	GraphConstructor<K>::Index *index;
	g_c.ConstructGraph(g, index);
	INFO("Graph condensed");

	INFO("Writing to file");
	WriteToFile("simulated_mistakes.dot", "simulated_mistakes_graph", *g);

	INFO("Counting stats");
	CountStats(*g);

	INFO("Clipping tips");
	TipComparator comparator(*g);
	TipClipper<TipComparator> tc(comparator);
	tc.ClipTips(*g);

	INFO("Writing to file");
	WriteToFile("tips_clipped.dot", "no_tips_graph", *g);
	delete g;
	delete index;
}

void ConstructGraphOnReads(vector<Read> reads) {
	INFO("Tool started");
	INFO("Constructing DeBruijn graph");
	DeBruijn<K> debruijn;
	debruijn.ConstructGraph(reads);
	INFO("DeBruijn graph constructed");
	CondenseTool(debruijn);
}

template <class ReadStream>
void ConstructGraphFromStream(ReadStream& stream) {
	INFO("Tool started");
	INFO("Constructing DeBruijn graph");
	DeBruijn<K> debruijn;
	debruijn.ConstructGraphFromStream(stream);
	INFO("DeBruijn graph constructed");
	CondenseTool(debruijn);
}


}


#endif /* EDGE_GRAPH_TOOL_HPP_ */

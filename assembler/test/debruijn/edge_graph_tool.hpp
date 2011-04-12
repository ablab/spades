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

void WriteToFile(const string& file_name, const string& graph_name,
		const EdgeGraph& g) {
	fstream filestr;
	filestr.open(file_name.c_str(), fstream::out);
	gvis::PairedGraphPrinter<VertexId> gp(
			"simulated_data_graph", filestr);
	ComplementGraphVisualizer gv(gp);
	gv.Visualize(g);
	filestr.close();

}

void SimulatedMistakesTool() {
	INFO("Tool started");
	vector<Read> reads = GenerateReadsWithMistakes();
	INFO("Constructing DeBruijn graph");
	DeBruijn<K> debruijn;
	debruijn.ConstructGraph(reads);
	INFO("DeBruijn graph constructed");

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
	delete g;
	delete index;
}

}


#endif /* EDGE_GRAPH_TOOL_HPP_ */

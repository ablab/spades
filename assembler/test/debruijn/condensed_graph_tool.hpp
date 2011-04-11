/*
 * condensed_graph_tool.hpp
 *
 *  Created on: Mar 25, 2011
 *      Author: sergey
 */

#ifndef CONDENSED_GRAPH_TOOL_HPP_
#define CONDENSED_GRAPH_TOOL_HPP_

#include "read_generator.hpp"
#include "condensed_graph_constructor.hpp"
#include "ireadstream.hpp"
#include <algorithm>

#define filename "./data/input/MG1655-K12.fasta.gz"
#define SUBSTR_LENGTH 1000
#define COVERAGE 1
#define R 35
#define K 15

using namespace std;

namespace condensed_graph {

void CountStats(const CondensedGraph& g) {
	DFS dfs(g);
	SimpleStatCounter stat_c;
	dfs.Traverse(stat_c);
	INFO("Vertex count=" << stat_c.v_count() << "; Edge count="
			<< stat_c.e_count());
}

void WriteToFile(const string& file_name, const string& graph_name,
		const CondensedGraph& g) {
	fstream filestr;
	filestr.open(file_name.c_str(), fstream::out);
	gvis::PairedGraphPrinter<const condensed_graph::Vertex*> gp(
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
	condensed_graph::CondenseConstructor<K> g_c(debruijn);

	condensed_graph::CondensedGraph *g;
	condensed_graph::SimpleIndex<K> *index;
	g_c.ConstructGraph(g, index);
	INFO("Graph condensed");

	INFO("Writing to file");
	condensed_graph::WriteToFile("simulated_mistakes.dot", "simulated_mistakes_graph", *g);

	INFO("Counting stats");
	condensed_graph::CountStats(*g);
//	delete g;
//	delete index;
}

}

#endif /* CONDENSED_GRAPH_TOOL_HPP_ */

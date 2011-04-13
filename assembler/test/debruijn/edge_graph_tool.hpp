/*
 * edge_graph_tool.hpp
 *
 *  Created on: Apr 10, 2011
 *      Author: sergey
 */

#ifndef EDGE_GRAPH_TOOL_HPP_
#define EDGE_GRAPH_TOOL_HPP_
<<<<<<< HEAD:assembler/test/debruijn/edge_graph_tool.hpp
#include "tip_clipper.hpp"
=======

>>>>>>> 7fd1100e920f1d86fca7dcd79637eaa633e8973c:assembler/test/debruijn/edge_graph_tool.hpp
namespace edge_graph {

void CountStats(const EdgeGraph& g) {
	de_bruijn::DFS<EdgeGraph> dfs(g);
	de_bruijn::SimpleStatCounter<EdgeGraph> stat_c;
	dfs.Traverse(&stat_c);
	INFO("Vertex count=" << stat_c.v_count() << "; Edge count="
			<< stat_c.e_count());
}

<<<<<<< HEAD:assembler/test/debruijn/edge_graph_tool.hpp
//void WriteToFile(const string& file_name, const string& graph_name,
//		const EdgeGraph& g) {
//	fstream filestr;
//	filestr.open(file_name.c_str(), fstream::out);
//	gvis::PairedGraphPrinter<VertexId> gp(
//			"simulated_data_graph", filestr);
//	ComplementGraphVisualizer gv(gp);
//	gv.Visualize(g);
//	filestr.close();
//}
=======
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
>>>>>>> 7fd1100e920f1d86fca7dcd79637eaa633e8973c:assembler/test/debruijn/edge_graph_tool.hpp

template <class ReadStream>
void ConstructionTool(ReadStream& stream) {
	INFO("Tool started");
	INFO("Constructing DeBruijn graph");
	DeBruijn<K> debruijn;
	debruijn.ConstructGraphFromStream(stream);
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

<<<<<<< HEAD:assembler/test/debruijn/edge_graph_tool.hpp
	INFO("Clipping tips");
	TipComparator comparator(*g);
	TipClipper<TipComparator> tc(comparator);
	cout << "create" << endl;
	tc.ClipTips(*g);

	INFO("Writing to file");
	WriteToFile("tips_clipped.dot", "no_tips_graph", *g);
	delete g;
	delete index;
}

=======
>>>>>>> 7fd1100e920f1d86fca7dcd79637eaa633e8973c:assembler/test/debruijn/edge_graph_tool.hpp
}


#endif /* EDGE_GRAPH_TOOL_HPP_ */

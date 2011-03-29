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
#include "ifaststream.hpp"

#define filename "some_name"
#define R 15
#define K 9

using namespace std;

namespace condensed_graph {


	void tool() {
		cerr << "Hello, I am assembler!" << endl;

		// read all 'read's

		cerr << "Reading " << filename << endl;
		ifaststream stream(filename);
		string name;
		string nucl_str;
		stream >> name >> nucl_str;

		stream.close();

		ReadGenerator<R> gen(nucl_str.substr(0, 4000), 15);

		vector<strobe_read<R> > reads;
		while (!gen.eof()) {
			strobe_read<R> read;
			gen >> read;
			reads.push_back(read);
		}

		DeBruijn<K> debruijn;
		debruijn.ConstructGraph(reads);
		condensed_graph::CondenseConstructor<K> g_c(debruijn);

		condensed_graph::CondensedGraph *g;
		condensed_graph::SimpleIndex<K> *index;
		g_c.ConstructGraph(g, index);
		fstream filestr;
		filestr.open("graph.dot", fstream::out);
		gvis::PairedGraphPrinter<const condensed_graph::Vertex*> gp(
				"simulated_data_graph", filestr);
		condensed_graph::ComplementGraphVisualizer gv(gp);
		gv.Visualize(*g);
		filestr.close();

		condensed_graph::DFS dfs(g);
		condensed_graph::SimpleStatCounter stat_c;
		dfs.Traverse(stat_c);
		cerr << "Vertex count=" << stat_c.v_count() << "; Edge count="
				<< stat_c.e_count() << endl;
	}
}

#endif /* CONDENSED_GRAPH_TOOL_HPP_ */

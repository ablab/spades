/*
 * assembler.cpp
 *
 *  Created on: 02.03.2011
 *      Author: vyahhi
 */

#include "ireadstream.hpp"
#include "condensedGraph.hpp"
#include "condensedGraphConstructor.hpp"
#include "debruijn.hpp"
#include <cassert>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <vector>
#include <ctime>
#include <string>
#include <stdio.h>

#define K 5//25
//todo make separate class to construct graph and remove R from here!!!
// read size:
#define R 9//100
#define N 9//100//11//100

using namespace std;

// input files:
//#define filename1 "./data/MG1655-K12_emul1.fasta.gz"
//#define filename2 "./data/MG1655-K12_emul2.fasta.gz"
#define filename1 "./test/data/s_6_1.fastq.gz"
#define filename2 "./test/data/s_6_2.fastq.gz"

int main(int argc, char *argv[]) {
	cerr << "Hello, I am assembler!" << endl;
	time_t now = time(NULL);

	// read all 'read's

	cerr << "Reading " << filename1 << " and " << filename2 << "..." << endl;
	ireadstream<R, 2, int> irs(filename1, filename2);
	vector<mate_read<R, int>::type> *v = irs.readAll(30000); // read not all `reads` (for faster debug)
	irs.close();
	cerr << "Total reads (mate, without Ns): " << v->size() << endl;
	cerr << "Current time: " << (time(NULL) - now) << " sec." << endl;

	// construct graph

	DeBruijn<K> debruijn;
	debruijn.ConstructGraph(*v);
	condensed_graph::CondenseConstructor<K> g_c(debruijn);

	condensed_graph::Graph *g;
	condensed_graph::SimpleHashTable<K> *index;
	g_c.ConstructGraph(g, index);
	fstream filestr;
	filestr.open("graph.dot", fstream::out);
	gvis::GraphPrinter<const condensed_graph::Vertex*> gp(
			"simulated data graph", filestr);
	condensed_graph::SimpleGraphVisualizer gv(gp);
	gv.Visualize(*g);
	filestr.close();

	condensed_graph::DFS dfs(g);
	condensed_graph::SimpleStatCounter stat_c;
	dfs.Traverse(stat_c);
	cerr << "Vertex count=" << stat_c.v_count() << "; Edge count="
			<< stat_c.e_count() << endl;

	/*
	 * Simple de Bruijn graph construction:
	 */
	/*cerr << "Constructing de Bruijn graph..." << endl;
	 DeBruijn<K> *graph = new DeBruijn<K>();
	 for (size_t i = 0; i < v->size(); ++i) {
	 for (size_t r = 0; r < 2; ++r) {
	 Seq<R> read = v->operator[](i)[r];
	 Seq<K> head = Seq<K>(read);
	 for (size_t j = K; j < R; ++j) {
	 Seq<K> tail = head << read[j];
	 graph->addEdge(head, tail);
	 head = tail;
	 }
	 }
	 }
	 cerr << "Total nodes: " << graph->size() << endl;*/
	cerr << "Current time: " << (time(NULL) - now) << " sec." << endl;

	// simplify graph

	// TODO
	cerr << "Current time: " << (time(NULL) - now) << " sec." << endl;

	// output graph

	cerr << "Total Time: " << (time(NULL) - now) << " sec." << endl;
	return 0;
}

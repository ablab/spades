/*
 * assembler.cpp
 *
 *  Created on: 02.03.2011
 *      Author: vyahhi
 */

#include "ireadstream.hpp"
#include "condensedGraph.hpp"
//#include "debruijn.hpp"
#include <cassert>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <vector>
#include <ctime>
#include <string>
#include <stdio.h>

using namespace std;

// read size:
#define R 100

// k-mer size:
//#define K 11

// input files:
#define filename1 "./data/MG1655-K12_emul1.fasta.gz"
#define filename2 "./data/MG1655-K12_emul2.fasta.gz"
//#define filename1 "./test/data/s_6_1.fastq.gz"
//#define filename2 "./test/data/s_6_2.fastq.gz"

int main(int argc, char *argv[]) {
	cerr << "Hello, I am assembler!" << endl;
	time_t now = time(NULL);

	// read all 'read's

	cerr << "Reading " << filename1 << " and " << filename2 << "..." << endl;
	ireadstream<R,2,int> irs(filename1, filename2);
	vector<mate_read<R,int>::type> *v = irs.readAll(30000); // read not all `reads` (for faster debug)
	irs.close();
	cerr << "Total reads (mate, without Ns): " << v->size() << endl;
	cerr << "Current time: " << (time(NULL) - now) << " sec." << endl;

	// construct graph

	time_t now2 = time(NULL);
	condensed_graph::Graph g;
	for (size_t i = 0; i < v->size(); ++i) {
		if (i % 10000 == 0) {
			cerr << "mate reads: " << i << ", time: " << (time(NULL) - now2) << endl;
		}
		g.ThreadRead((*v)[i][0]);
		g.ThreadRead((*v)[i][1]);
	}

	fstream filestr;
	filestr.open("graph.dot", fstream::out);
	gvis::OnlineGraphPrinter<condensed_graph::Vertex*> gp("simulated data graph", filestr);
	condensed_graph::SimpleGraphVisualizer gv(gp);
	gv.Visualize(g);
	filestr.close();

	condensed_graph::DFS dfs(g);
	condensed_graph::SimpleStatCounter h;
	dfs.Traverse(h);
	cerr<<"Vertex count="<<h.v_count()<<"; Edge count="<<h.e_count() << endl;

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

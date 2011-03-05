/*
 * assembler.cpp
 *
 *  Created on: 02.03.2011
 *      Author: vyahhi
 */

#include "ireadstream.hpp"
#include "condensedGraph.hpp"
#include <cassert>
#include <iostream>
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
//#define filename1 "./data/MG1655-K12_emul1.fasta.gz"
//#define filename2 "./data/MG1655-K12_emul2.fasta.gz"
#define filename1 "./test/data/s_6_1.fastq.gz"
#define filename2 "./test/data/s_6_2.fastq.gz"

int main(int argc, char *argv[]) {
	cerr << "Hello, I am assembler!" << endl;
	time_t now = time(NULL);

	// read all 'read's

	ireadstream<R,2,int> irs(filename1, filename2);
	vector<mate_read<R,int>::type> *v = irs.readAll();
	irs.close();
	cerr << "Total reads (mate, without Ns): " << v->size() << endl;
	cerr << "Current time: " << (time(NULL) - now) << " sec." << endl;

	// construct graph

	/*DeBruijn<K> *graph = new DeBruijn<K>();
	for (int i = 0; i < v->size; ++i) {
		for (int i = 0; )
	}*/
	condensed_graph::Graph g;
	for (size_t i = 0; i < v->size(); ++i) {
		g.ThreadRead(v->operator [](i)[0]);
		g.ThreadRead(v->operator [](i)[1]);
//		cout << v->operator [](i)[0].str() << endl;
//		cout << v->operator [](i)[1].str() << endl;
	}

	cerr << "Current time: " << (time(NULL) - now) << " sec." << endl;


	// simplify graph

	// TODO
	cerr << "Current time: " << (time(NULL) - now) << " sec." << endl;

	// output graph

	cerr << "Total Time: " << (time(NULL) - now) << " sec." << endl;
	return 0;
}

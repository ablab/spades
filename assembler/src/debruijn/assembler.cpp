#include <cassert>
#include <iostream>
#include <list>
#include <cstdio>
#include <vector>
#include <ctime>
#include <string>
#include "../parser.hpp"
#include "debruijn.hpp"

using namespace std;

pair<string,string> filenames = make_pair("./data/s_6_1.fastq.gz", "./data/s_6_2.fastq.gz");

#define MPSIZE 100
#define K 11

int main(int argc, char *argv[]) {

	std::cerr << "Hello, I am assembler!" << std::endl;

	time_t now = time(NULL);

	Seq<3> s = Seq<3>("ACG");
	s.shift_right(0);
	//return 0;

	// simple de Bruijn graph
	DeBruijn<K> graph;
	// start parsing...
	FASTQParser<MPSIZE>* fqp = new FASTQParser<MPSIZE>();
	fqp->open(filenames.first, filenames.second);
	int cnt = 0;
	while (!fqp->eof()) {
		MatePair<MPSIZE> mp = fqp->read(); // is it copy? :)
		if (mp.id != -1) { // don't have 'N' in reads
			Seq<K> head = mp.seq1.head<K>();
			Seq<K> tail;
			for (size_t i = K; i < MPSIZE; ++i) {
				cerr << head.str() << endl;
				tail = head.shift_right(mp.seq1[i]);
				graph.addEdge(head, tail);
				head = tail;
			}
		}
		cnt++;
	}
	cout << "Total reads: " << cnt << endl;
	cout << "Total nodes: " << graph._nodes.size() << endl;
	//cout << "Clear (without N) reads: " << mps.size() << endl;
	cout << "seconds: " << (time(NULL) - now) << endl;
	fqp->close();

	return 0;
}

#include "cute.h"
#include <set>
#include "seq.hpp"
#include "sequence.hpp"
#include "nucl.hpp"
#include "debruijn.hpp"
#include "graphVisualizer.hpp"
#include <fstream>
#include <iostream>

using namespace std;

//void TestAddNode() {
//	DeBruijn<5> g;
//	Seq<5, int> seq1("ACAAA");
//	Seq<5, int> seq2("CAAAC");
//	Seq<5, int> seq3("CAAAA");
//	g.addNode(seq1);
//	g.addNode(seq2);
//	g.addNode(seq3);
//	int c = 0;
//	for (DeBruijn<5>::kmer_iterator it = g.kmer_begin(); it != g.kmer_end(); it++) {
//		c++;
//		ASSERT(*it == seq1 || *it == seq2 || *it == seq3);
//	}
//	ASSERT_EQUAL(3, c);
//}
namespace de_bruijn {

void TestAddEdge() {
//	DeBruijn<5> g;
//	Seq<5> seq1("ACAAA");
//	Seq<5> seq2("CAAAG");
//	Seq<5> seq3("CAAAA");
//	g.addEdge(seq1, seq2);
//	g.addEdge(seq1, seq3);
//	int c = 0;
//	for (DeBruijn<5>::kmer_iterator it = g.begin(); it != g.end(); ++it) {
//		c++;
//		ASSERT(*it == seq1 || *it == seq2 || *it == seq3);
//	}
//	ASSERT_EQUAL(3, c);
//	ASSERT_EQUAL(2, g.OutgoingEdgeCount(seq1));
//	DeBruijn<5>::edge_iterator n_it = g.OutgoingEdges(seq1);
//	ASSERT_EQUAL(seq3, (*n_it).end<5>());
//	++n_it;
//	ASSERT_EQUAL(seq2, (*n_it).end<5>());
}

void TestAddEdge2() {
//	DeBruijn<5> g;
//	Seq<5> seq1("ACAAA");
//	Seq<5> seq2("CCAAA");
//	Seq<5> seq3("CAAAC");
//	g.addEdge(seq1, seq3);
//	g.addEdge(seq2, seq3);
//	int c = 0;
//	for (DeBruijn<5>::kmer_iterator it = g.begin(); it != g.end(); it++) {
//		c++;
//		ASSERT(*it == seq1 || *it == seq2 || *it == seq3);
//	}
//	ASSERT_EQUAL(3, c);
//	ASSERT_EQUAL(2, g.IncomingEdgeCount(seq3));
//	DeBruijn<5>::edge_iterator n_it = g.IncomingEdges(seq3);
//	ASSERT_EQUAL(seq1, (*n_it).start<5>());
//	++n_it;
//	ASSERT_EQUAL(seq2, (*n_it).start<5>());
}

void TestSimpleConstruction() {
//	string ss[] = { "CGAAACCAC", "CGAAAACAC", "AACCACACC", "AAACACACC" };
//	vector<Read> input;
//	for (int i = 0; i < 4; ++i) {
//		input.push_back(Read("noname", ss[i], "noqual"));
//	}
//	DeBruijn<5> g;
//	g.ConstructGraph(input);
//	int c = 0;
//	Seq<5> seq("CGAAA");
//	Seq<5> seq2("CACAC");
//	for (DeBruijn<5>::kmer_iterator it = g.begin(); it != g.end(); it++) {
//		c++;
//	}
//	ASSERT_EQUAL(26, c);
//	ASSERT_EQUAL(2, g.OutgoingEdgeCount(seq));
//	DeBruijn<5>::edge_iterator n_it = g.OutgoingEdges(seq);
//	ASSERT_EQUAL(Seq<5>("GAAAA"), (*n_it).end<5>());
//	++n_it;
//	ASSERT_EQUAL(Seq<5>("GAAAC"), (*n_it).end<5>());
//
//	ASSERT_EQUAL(2, g.IncomingEdgeCount(seq2));
//	DeBruijn<5>::edge_iterator n_it2 = g.IncomingEdges(seq2);
//	ASSERT(Seq<6>("ACACAC") == *n_it2);
//	++n_it2;
//	ASSERT_EQUAL(Seq<6>("CCACAC"), *n_it2);
}

cute::suite DeBruijnGraphSuite() {
	cute::suite s;
//	s.push_back(CUTE(TestAddNode));
	s.push_back(CUTE(TestAddEdge));
	s.push_back(CUTE(TestAddEdge2));
	s.push_back(CUTE(TestSimpleConstruction));
	return s;
}

}

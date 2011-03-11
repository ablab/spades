#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include "abruijngraph.hpp"
#include "hash.hpp"
#include "graphBuilder.hpp"
#include "parameters.hpp"
#include <ext/hash_map>
#include "logging.hpp"

LOGGER("a.graphBuilder");

using namespace std;
using namespace __gnu_cxx;
using namespace abruijn;

typedef hash_map< Sequence, int, HashSym<Sequence>, EqSym<Sequence> > SeqCount;
SeqCount seqCount;

abruijn::Graph graph;

LoggerPtr logger(Logger::getLogger("a.graphBuilder"));
HashSym<Sequence> hashSym;

void processReadA(Seq<MPSIZE> r) {
	// Processing in O(length * k).
	// Can be done in O(length), provided the hash is polynomial.
	unsigned int h1 = -1;
	int pos1 = -1;
	unsigned int h2 = -1;
	int pos2 = -1;
	for (int i = 0; i + K <= MPSIZE; i++) {
		unsigned int h = hashSym(Sequence(r).Subseq(i, i + K));
		if (h < h1) {
			h2 = h1;
			pos2 = pos1;
			h1 = h;
			pos1 = i;
		} else if (h < h2) {
			h2 = h;
			pos2 = i;
		}
	}
	if (pos1 > pos2) {
		int t = pos1; pos1 = pos2; pos2 = t;
		t = h1; h1 = h2; h2 = t;
	}
	seqCount[Sequence(r).Subseq(pos1, pos1 + K)]++;
	seqCount[Sequence(r).Subseq(pos2, pos2 + K)]++;
}

void GraphBuilder::selectGood()
{
    for (SeqCount::iterator p = seqCount.begin(); p != seqCount.end(); ++p) {
    	graph.addVertex(&(p->first));
	}
}

void processReadB(Seq<MPSIZE> r) {
	vector<abruijn::Vertex*> vs;
	vector<int> index;
	for (int i = 0; i + K <= MPSIZE; i++) {
		Sequence s(Sequence(r).Subseq(i, i + K));
		if (graph.hasVertex(&s)) {
			vs.push_back(graph.getVertex(&s));
			index.push_back(i);
		}
	}
	DEBUG(vs.size());
	for (size_t i = 0; i + 1 < vs.size(); i++) {
		graph.addEdge(*vs[i], *vs[i + 1], index[i + 1] - index[i]);
	}
}

void condense() {

}

void GraphBuilder::build() {
	INFO("Building graph...");
	mate_read<MPSIZE>::type mp;
	for (int i = 0; !irs.eof() && i < CUT; i++) {
		irs >> mp;
		processReadA(mp[0]);
		processReadA(mp[1]);
	}
	INFO("processReadA done: " << seqCount.size() << " vertices");
	irs.reset();
    selectGood();
	INFO("selectGood done: " << graph.vertexArray.size() << " vertices");
	for (int i = 0; !irs.eof() && i < CUT; i++) {
		irs >> mp;
		processReadB(mp[0]);
		processReadB(mp[1]);
	}
	irs.close();
	condense();
    graph.output("output\\z.dot");
}

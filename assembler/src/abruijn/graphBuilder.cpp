#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <list>
#include <set>
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

//typedef hash_map< Sequence, int, HashSym<Sequence>, EqSym<Sequence> > SeqCount;
//SeqCount seqCount;
set<hash_t> goodHashes;

abruijn::Graph graph;

LoggerPtr logger(Logger::getLogger("a.graphBuilder"));
HashSym<Sequence> hashSym;

static hash_t H;

void initH() {
	H = 1;
	for (int i = 0; i < K; i++) {
		H = HASH_X(H);
	}
}

hash_t ha[MPSIZE];
hash_t hb[MPSIZE];

void countHashes(Sequence s) {
	size_t sz = s.size();
	hash_t h = 0;
	for (int i = 0; i < K; i++) {
		h = HASH_X(h) + s[i];
	}
	ha[0] = h;
	for (size_t i = 0; i + K < sz; i++) {
		ha[i + 1] = HASH_X(ha[i]) + s[i + K] - s[i] * H;
	}

	h = 0;
	for (size_t i = sz - 1; i + K >= sz; i--) {
		h = HASH_X(h) + (s[i] ^ 3);
	}
	hb[sz - K] = h;
	for (int i = sz - K; i > 0; i--) {
		hb[i - 1] = HASH_X(hb[i]) + (s[i - 1] ^ 3) - (s[i + K - 1] ^ 3) * H;
		ha[i] ^= hb[i] ^ HASH_XOR;
	}
	ha[0] ^= hb[0] ^ HASH_XOR;
}

void processReadA(Seq<MPSIZE> r) {
	countHashes(Sequence(r));
	hash_t h1 = -1;
//	int pos1 = -1;
	hash_t h2 = -1;
//	int pos2 = -1;
	for (int i = 0; i + K <= MPSIZE; i++) {
		if (ha[i] < h1) {
			h2 = h1;
//			pos2 = pos1;
			h1 = ha[i];
//			pos1 = i;
		} else if (ha[i] < h2) {
			h2 = ha[i];
//			pos2 = i;
		}
	}
	goodHashes.insert(h1);
	goodHashes.insert(h2);
}

void selectGood() {
//    for (SeqCount::iterator p = seqCount.begin(); p != seqCount.end(); ++p) {
//    	graph.createVertex(&(p->first));
//	}
}

void processReadB(Seq<MPSIZE> r) {
	vector<abruijn::Vertex*> vs;
	vector<int> index;
	Sequence s(r);
	countHashes(s);
	for (int i = 0; i + K <= MPSIZE; i++) {
		if (goodHashes.find(ha[i]) != goodHashes.end()) {
			//			Sequence* ss = new Sequence(s.Subseq(i, i + K));
			vs.push_back(graph.getVertex(s.Subseq(i, i + K)));
			index.push_back(i);
		}
	}
	for (size_t i = 0; i + 1 < vs.size(); i++) {
		graph.addEdge(vs[i], vs[i + 1], index[i + 1] - index[i] + K);
	}
}


void condenseA() {
	int condensations = 0;
	list<Vertex*> vect(graph.vertices.begin(), graph.vertices.end());
	for (list<Vertex*>::iterator v = vect.begin(); v != vect.end(); ++v) {
		if (graph.vertices.find(*v) == graph.vertices.end()) {
			continue;
		}
		if ((*v)->degree() == 1) {
			if ((*v)->edges_.begin()->first->complement_->degree() == 1) {
				VERBOSE(condensations++, " condensations");
				TRACE("CondenseA " << (*v)->size() << " " << (*v)->edges_.begin()->first->size());
				Vertex* u = graph.condense(*v);
				if (u != NULL) {
					vect.push_back(u);
				}
			}
		}
	}
}

void GraphBuilder::build() {
	initH();

	INFO("Processing-A...");
	mate_read<MPSIZE>::type mp;
	for (int i = 0; !irs.eof() && i < CUT; i++) {
		irs >> mp;
		processReadA(mp[0]);
		processReadA(mp[1]);
		VERBOSE(i, " reads read");
	}
	irs.reset();
	INFO("processReadA done: " << goodHashes.size() << " vertex-pairs");

	INFO("Selecting good kmers...");
	selectGood();

	INFO("Processing-B...");
	for (int i = 0; !irs.eof() && i < CUT; i++) {
		irs >> mp;
		processReadB(mp[0]);
		processReadB(mp[1]);
		VERBOSE(i, " reads processed");
	}
	irs.close();

	INFO("Condensing graph...");
	condenseA();

	INFO("Outputting graph to " << OUTPUT_FILE);
	graph.output(OUTPUT_FILE + ".dot");
}

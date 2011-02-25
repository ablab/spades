#include <iostream>
#include "graph.h"
#include "parameters.h"
#include "../parser.hpp"
#include <ext/hash_map>

using namespace std;
using namespace __gnu_cxx;

// seq[0]*31^(n-1) + seq[1]*31^(n-2) + ... + seq[n-1]*31^0
template< typename T >
struct Hash {
	unsigned int operator() (const T* seq) const {
	    unsigned int h = HASH_SEED;
		for (int i = 0; i < seq->len(); i++) {
			h = ((h << 5) - h) + (*seq)[i];
		}
		return h;
	}
};

template< typename T >
struct HashSym {
	unsigned int operator() (const T* seq) const {
		Hash<T> h;
		return h(seq) ^ h(!seq);
	}
};

struct HashSymWeighted {
  unsigned int operator() (const Kmer* seq) {
    return 0;
    // Will use frequency of seq
  }
};

hash_map<Sequence, CVertex, HashSym<Kmer> > kmers_map;
CGraph graph(1000000);

void processRead(Seq<MPSIZE> r) {
	// Processing in O(length * k).
	// Can be done in O(length), provided the hash is polynomial.
	unsigned int h1 = -1;
	int pos1 = -1;
	unsigned int h2 = -1;
	int pos2 = -1;
	for (int i = 0; i + K <= MPSIZE; i++) {
		Hash<Sequence> hh;
		Sequence* kmer = &r.substring(i, i + K);
		unsigned int h = hh(kmer);
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
//	Kmer k1 = r.substring(pos1, pos1 + K);
	CVertex v1 = kmers_map[r.substring(pos1, pos1 + K)];
	CVertex v2 = kmers_map[r.substring(pos2, pos2 + K)];
	graph.AddVertex(v1);
	graph.AddVertex(v2);
	CEdge e(v2, r, pos1, pos2);
	v1.AddEdge(e);
}

CGraph build() {
	FASTQParser<MPSIZE>* fqp = new FASTQParser<MPSIZE>();
	fqp->open(filenames.first, filenames.second);
	while (!fqp->eof()) {
		MatePair<MPSIZE> r = fqp->read(); // is it copy? :)
		if (r.id == -1) { // don't have 'N' in reads
			continue;
		}
		processRead(r.seq1);
		processRead(r.seq2);
		//cout <<  mp.id << endl << mp.seq1.str() << endl <<  mp.seq2.str() << endl;
	}
	return graph;
}

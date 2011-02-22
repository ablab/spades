#include <iostream>
#include "graph.h"
#include "parameters.h"
#include "../parser.hpp"
#include <ext/hash_map>
namespace std { using namespace __gnu_cxx; }

using namespace std;

// seq[0]*31^(n-1) + seq[1]*31^(n-2) + ... + seq[n-1]*31^0
struct Hash {
	unsigned int operator() (const Kmer* seq) const {
	    unsigned int h = HASH_SEED;
		for (int i = 0; i < seq->len(); i++) {
			h = ((h << 5) - h) + (*seq)[i];
		}
		return h;
	}
};

struct HashSym {
	unsigned int operator() (const Kmer* seq) const {
		return Hash()(seq);// ^ Hash()(!seq); TODO
	}
};

struct HashSymWeighted {
	unsigned int operator() (const Kmer* seq) {
		return 0;
		// Will use frequency of seq
	}
};

int main() {
	cout << "Hello!";
	hash_map<Kmer, CVertex, HashSym> kmers_map;
	FASTQParser<MPSIZE>* fqp = new FASTQParser<MPSIZE>();
	fqp->open(filenames.first, filenames.second);
	int cnt = 0;
	while (!fqp->eof()) {
		MatePair<MPSIZE> r = fqp->read(); // is it copy? :)
		if (r.id == -1) { // don't have 'N' in reads
			continue;
		}
		//cout <<  mp.id << endl << mp.seq1.str() << endl <<  mp.seq2.str() << endl;
		cnt++;
		// Processing in O(length * k).
		// Can be done in O(length), provided the hash is polynomial.
		unsigned int h1 = -1;
		int pos1 = -1;
		unsigned int h2 = -1;
		int pos2 = -1;
		for (int i = 0; i + K <= MPSIZE; i++) {
			unsigned int h = 0;// TODO Hash(r.substring(i, i + K));
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
		//TODO
		//v1 = kmers_map[r.substring(pos1, pos1 + K)];
		//v2 = kmers_map[r.substring(pos2, pos2 + K)];
		//v1.add_edge(v2, r, pos1, pos2);
	}
	return 0;
}

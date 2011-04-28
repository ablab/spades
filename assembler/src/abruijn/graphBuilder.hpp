#ifndef GRAPHBUILDER_H_
#define GRAPHBUILDER_H_

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <list>
#include <set>
#include "hash.hpp"
#include "parameters.hpp"
#include "logging.hpp"
#include "abruijngraph.hpp"

namespace abruijn {

using namespace std;
using namespace __gnu_cxx;
using hashing::hash_t;

class GraphBuilder
{
	//typedef hash_map< Sequence, int, HashSym<Sequence>, EqSym<Sequence> > SeqCount;
	//SeqCount seqCount;
	/**
	 * Stores bitmask:
	 * 1 = the lexicographically lesser k-mer has been seen as non-rightmost k-mer in read
	 * 2 = the greater k-mer ...
	 */
	map<hash_t, char> has_right;
	map<hash_t, char> tips;
	set<hash_t> earmarked_hashes;
	map<hash_t, set<hash_t> > tip_extensions;

	Graph graph;

	hashing::HashSym<Sequence> hashSym;
	typedef vector<hash_t> hash_vector;
	hash_vector ha;
	hash_vector hbest;

	size_t htake;
	void findMinimizers(Sequence s);
	void findLocalMinimizers(Sequence s, size_t window_size);
	void findSecondMinimizer(Sequence s);
	void revealTips(Sequence s);
	void findTipExtensions(Sequence s);
	void lookRight(Sequence s);
	void addToGraph(Sequence s);

public:
	GraphBuilder(size_t htake) : htake(htake) {
		hbest.reserve(htake);
	};
	Graph build();
};

}

#endif /* GRAPHBUILDER_H_ */

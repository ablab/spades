#ifndef _abruijngraph_h_
#define _abruijngraph_h_

#include <cassert>
#include <ostream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <ext/hash_map>
#include "seq.hpp"
#include "sequence.hpp"
#include "parameters.hpp"
#include "hash.hpp"
#include "logging.hpp"

using namespace std;
using namespace __gnu_cxx;

namespace abruijn {

class Edge {
public:
	map<int, int> lengths_;

	void addLength(int len) {
		lengths_[len]++;
	}

	string toString() {
		stringstream ss;
		for (map<int, int>::iterator it = lengths_.begin(); it != lengths_.end(); ++it) {
			ss << " " << it->first << " {" << it->second << "};";
		}
		return ss.str();
	}
};

class Vertex {
public:
	const Sequence* kmer_;
	Vertex* complement_;
	typedef map<Vertex*, Edge> Edges;
	Edges edges_;
	Vertex(const Sequence* kmer) : kmer_(kmer) {};

	void addEdge(Vertex* to, int len) {
		edges_[to].addLength(len);
	}
	int degree() {
		return edges_.size();
	}
};

class Graph {
public:
	typedef hash_map < Sequence, Vertex*, HashSym<Sequence>, EqSym<Sequence> > SeqVertice;
	SeqVertice seqVertice;

	typedef set<Vertex*> Vertices;
	Vertices vertices;

	Graph() {}
	Vertex* createVertex(const Sequence* kmer);
	void addEdge(Vertex* from, Vertex* to, int len);
	void removeVertex(Vertex* v);
	void removeVertex_single(Vertex* v);
	bool hasVertex(const Sequence* kmer);
	Vertex* getVertex(const Sequence* kmer);

	Vertex* condense(Vertex* v);

	void output(std::ofstream &out);
	void output(string filename);
};

}

#endif // _abruijngraph_h_

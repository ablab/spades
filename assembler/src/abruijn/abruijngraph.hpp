#ifndef _abruijngraph_h_
#define _abruijngraph_h_

#include <vector>
#include <ostream>
#include <cassert>
#include <algorithm>
#include <ext/hash_map>
#include "seq.hpp"
#include "sequence.hpp"
#include "parameters.hpp"
#include "hash.hpp"
#include "logging.hpp"

using namespace std;
using namespace __gnu_cxx;

namespace abruijn {

class Edge;

class Vertex {
public:
	const Sequence* kmer_;
	Vertex* complement_;
	typedef vector<Edge> EdgeArray;
	EdgeArray edges_;
	Vertex(const Sequence* kmer) : kmer_(kmer) {};
	void addEdge(Edge* e) {
		edges_.push_back(*e);
	}
};

class Edge {
public:
	Vertex* to_;
	int len_;
	Edge (Vertex* to, int len) : to_(to), len_(len) {};
};

class Graph {
public:
	typedef hash_map < Sequence, Vertex*, HashSym<Sequence>, EqSym<Sequence> > SeqVertice;
	SeqVertice seqVertice;

	typedef vector <Vertex*> VertexArray;
	VertexArray vertexArray;

	Graph() {}
	void addVertex(const Sequence* kmer);
	void addEdge(Vertex &from, Vertex &to, int len);
	bool hasVertex(const Sequence* kmer);
	Vertex* getVertex(const Sequence* kmer);
	void output(const char* filename);
};

}

#endif // _abruijngraph_h_

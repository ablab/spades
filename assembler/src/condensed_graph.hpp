/**
 * condensed_graph.h
 *
 *  Created on: Feb 21, 2011
 *      Author: sergey
 */
#include <vector>
#include <ext/hash_map>
#include <cstring>
#include "seq.hpp"

using namespace std;
using namespace __gnu_cxx;

#ifndef CONDENSED_GRAPH_H_
#define CONDENSED_GRAPH_H_

typedef int _v_idx;

namespace assembler {

#define K 25
#define HASH_SEED 1845724623

typedef Seq<K> Kmer;

static const int k = 25;

char toIndex(char c) {
	switch (c) {
	case 'A':
		return 0;
	case 'C':
		return 1;
	case 'G':
		return 2;
	case 'T':
		return 3;
	default:
		break;//throw exception
	}
}

class Vertex;

class Vertex {
	Sequence _nucls;
	Vertex* _arcs[4];
	Vertex* _complement;

	int _coverage;
	int _arc_coverage[4];
public:
	//todo talk with Kolya about problem with Sequence copying!!!
	Vertex(Sequence nucls, Vertex** arcs, Vertex* complement) :
		_nucls(nucls), _complement(complement) {
		memcpy(arcs, _arcs, 4 * sizeof(Vertex*));
		//		_arcs = arcs;
	}
	;

	//static Vertex AbsentVertex = Vertex(0, 0, NULL, true, 0, NULL);
	int nucl_count() {
		return _nucls.len();
	}
	;

	char operator[](const int &index) const {
		return _nucls[index];
	}

	int coverage() {
		return _coverage;
	}
	;

	int arc_count() {
		int c = 0;
		for (int i = 0; i < 4; ++i)
			if (_arcs[i] != NULL)
				c++;
		return c;
	}
	;

	Vertex** arcs() {
		return _arcs;
	}
	;

	Vertex* Arc(char nucl) {
		return _arcs[(int)nucl];
	}

	Vertex* complement() {
		return _complement;
	}
	;
};

class Graph {
	vector<Vertex*> _component_roots;
public:
	vector<Vertex*> component_roots();
};

/*class KmerPos {
	//vertex, offset and direction of the canonical
	Vertex* _v;
	int _offset;
public:
	KmerPos(Vertex* v, int offset) :
		_v(v), _offset(offset) {
	}
};*/

class SimpleHash {
public:
	unsigned int operator() (const Kmer& seq) const {
	    unsigned int h = HASH_SEED;
		for (int i = 0; i < seq.len(); i++) {
			h = ((h << 5) - h) + seq[i];
		}
		return h;
	}
};

class MySimpleHashTable {
	//pair<Vertex*, int> - vertex and offset
//	hash_map<const Kmer, pair <Vertex*, int> , SimpleHash> h;
	hash_map<const Kmer, pair<Vertex*, int>, SimpleHash > h;
	//vector<V> array[size];
public:
	//todo think of using references
	void put(Kmer k, pair <Vertex*, int> v) {

	}

	const pair <Vertex*, int>& get(Kmer k);
	void remove(Kmer k);
};

/*class VertexPool {
 _v_idx _size;
 bool* _free;
 Vertex* _vertices;
 _v_idx _max_idx;
 public:
 VertexPool(_v_idx size);
 ~VertexPool();
 Vertex& operator[](_v_idx index) const;
 bool IsFree(_v_idx index) const;
 void Free(_v_idx index);
 _v_idx AllocatePair();
 };*/

}
#endif /* CONDENSED_GRAPH_H_ */

/**
 * condensed_graph.h
 *
 *  Created on: Feb 21, 2011
 *      Author: sergey
 */
#include <vector>
#include "seq.hpp"

using namespace std;

#ifndef CONDENSED_GRAPH_H_
#define CONDENSED_GRAPH_H_

namespace assembler {
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

/*char ComplementNucl(char c) {
 switch (c) {
 case 'A':
 return 'T';
 case 'C':
 return 'G';
 case 'G':
 return 'C';
 case 'T':
 return 'A';
 default:
 return 'N';
 }
 }*/

class Vertex;

class Arc {
	int _coverage;
	Vertex* _head;
public:
	Arc(int coverage, Vertex* head) :
		_coverage(coverage), _head(head) {
	}

	//static
	int coverage() {
		return _coverage;
	}
	;
	Vertex* head() {
		return _head;
	}
	;
};

class Vertex {
	int _coverage;
	Sequence* _nucls;
	bool _direction;
	int _arc_count;
	Arc* _arcs;
	//delete this pointer later
	Vertex* _complement;
public:
	Vertex(int coverage, int nucl_count, Sequence* nucls, bool direction,
			int arc_count, Arc* arcs, Vertex* complement) :
		_coverage(coverage), _nucls(nucls),
				_direction(direction), _arc_count(_arc_count), _arcs(arcs), _complement(complement) {
	}

	//static Vertex AbsentVertex = Vertex(0, 0, NULL, true, 0, NULL);
	int nucl_count() {
		return _nucls->len();
	}
	;

	/*Nucl operator[](const int &index) const {
		return (*_nucls)[index];
	}
	;*/

	int coverage() {
		return _coverage;
	}
	;

	int arc_count() {
		return _arc_count;
	}
	;

	Arc* arcs() {
		return _arcs;
	}
	;

	Arc* FindArc(char nucl) {
		return _arcs + toIndex(nucl);
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

class KmerPos {
	//vertex, offset and direction of the canonical
	Vertex* _v;
	int _offset;
	bool _direction;
public:
	KmerPos(Vertex* v, int offset, bool direction): _v(v), _offset(offset), _direction(direction)
	{}
};

class MyHashTable {
	int _array_size;
	//pair<Vertex*, int> - vertex and offset
	vector<KmerPos>* array;
};

}
#endif /* CONDENSED_GRAPH_H_ */

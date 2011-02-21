/**
 * condensed_graph.h
 *
 *  Created on: Feb 21, 2011
 *      Author: sergey
 */
#include <vector>

using namespace std;

#ifndef CONDENSED_GRAPH_H_
#define CONDENSED_GRAPH_H_

namespace assembler {
static const int k = 25;

/*enum Nucl {
	A, T, G, C
};*/

char ComplementNucl(char c) {
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
}

class Vertex;

class Arc {
	int _coverage;
	Vertex* _head;
public:
	Arc(int coverage, Vertex* head) : _coverage(coverage), _head(head)
	{}

	//static
	int coverage() {return _coverage;};
	Vertex* head() {return _head;};
};

class Vertex {
	int _coverage;
	//int _nucl_count;
	SeqVarLen* _nucls;
	bool _direction;
	int _arc_count;
	Arc* _arcs;
public:
	Vertex(int coverage, int nucl_count, char* nucls, bool direction
			, int arc_count, Arc* arcs):
	_coverage(coverage), _nucl_count(nucl_count), _nucls(nucls), _direction(direction)
	, _arc_count(arc_count), _arcs(arcs)
	{}

	//static Vertex AbsentVertex = Vertex(0, 0, NULL, true, 0, NULL);
	int nucl_count() {return _nucls->count();};

	char operator[](const int &index) const {
		if (_direction) {
			return _nucls[index];
		} else {
			return ComplementNucl(_nucls[_nucls->len() - 1 - index]);
		}
	};
	int coverage() {return _coverage;};
	int arc_count() {return _arc_count;};
	Arc* arcs() {return _arcs;};
	Arc* FindArc(char nucl);
	Vertex* Complement() {};
};

class Graph {
	vector<Vertex*> _component_roots;
public:
	vector<Vertex*> component_roots();
};
}
#endif /* CONDENSED_GRAPH_H_ */

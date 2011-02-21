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

static const int k = 25;

enum Nucl {
	A, T, G, C
};

class Vertex {
	int _coverage;
	int _nucl_count;
	char* _nucls;
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
	int nucl_count();
	Nucl operator[](const int &index) const;
	int coverage();
	int arc_count();
	void ArcBounds(Arc* start, Arc* end);
	Arc* FindArc(Nucl& nucl);
	Vertex* Complement();
};

class Arc {
	int _coverage;
	Vertex* _head;
public:
	Arc(int coverage, Vertex* head) : _coverage(coverage), _head(head)
	{}

	//static
	int coverage();
	Vertex* head();
};

class Graph {
	vector<Vertex*> component_roots;
public:
	vector<Vertex*> component_roots();
};

#endif /* CONDENSED_GRAPH_H_ */

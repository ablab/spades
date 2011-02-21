/*
 * condensed_graph.h
 *
 *  Created on: Feb 21, 2011
 *      Author: sergey
 */

#ifndef CONDENSED_GRAPH_H_
#define CONDENSED_GRAPH_H_

enum Nucl {

};

class NuclIterator {
	int nuclCount;
	int currCharPos;
	int currNuclCount;
public:
	void next
};

class Nucls {
	int nuclCount;
	char* nucls;
public:
	NuclIterator iter(bool direction);
	NuclIterator iter(bool direction, int offset);
};

class Vertex {
	int coverage;
	Nucls& nucls;
	bool direction;
	int arcCount;
	Arc* arcs;
public:
	int getCoverage();
	int getArcCount();
	void getArcBounds(Arc* start, Arc* end);
	Arc& getArc(Nucl& nucl);
	NuclIterator nuclIter();
	NuclIterator nucKOffsetIter();
	Vertex* complement();
};

class Arc {
	int coverage;
	Vertex* v;
public:
	int getCoverage();
	Vertex* getV();
};

class Graph {
	int rootCount;
	Vertex* roots;
};

#endif /* CONDENSED_GRAPH_H_ */

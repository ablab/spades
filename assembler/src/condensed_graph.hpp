/*
 * condensed_graph.h
 *
 *  Created on: Feb 21, 2011
 *      Author: sergey
 */

#ifndef CONDENSED_GRAPH_H_
#define CONDENSED_GRAPH_H_

class Nucl {

};

class nucl_iterator {

};

class Nucls {
	int nuclCount;
	char* nucls;
public:
	nucl_iterator iter(bool direction);
	nucl_iterator iter(bool direction, int offset);
};

class Vertex {
	int coverage;
	Nucls& nucls;
	bool direction;
	int arcCount;
	Arc* arcs;
public:
	int getCoverage();
	void getArcCount();
	void getArcBounds(Arc* start, Arc* end);
	Arc& getArc(Nucl& nucl);
	nucl_iterator nuclIter();
	nucl_iterator nuclChainIter();
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

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

class Vertex {
	int coverage;
	bool straight;
	int nuclCount;
	char* nucls;
	int arcCount;
	Arc* arcs;
public:
	Arc* getArcs();
	Arc& getArc(Nucl& nucl);
};

class Arc {
	int coverage;
	Vertex* v;
public:
	int getCoverage();
	Vertex* getV();
};

#endif /* CONDENSED_GRAPH_H_ */

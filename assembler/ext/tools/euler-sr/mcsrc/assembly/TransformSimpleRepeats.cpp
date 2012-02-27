/***************************************************************************
 * Title:          TransformSimpleRepeats.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "IntervalGraph.h"
#include <map>
#include "IntegralTupleStatic.h"



ssize_t IsSimpleRepeat(IntervalGraph &g, ssize_t e);

int main(int argc, char* argv[]) {

	std::string graphBase = argv[1];
	int vertexSize = atoi(argv[2]);
	IntervalGraph graph;
	ReadIntervalGraph(graphBase, graph, vertexSize);
	
	ssize_t e;
	for (e = 0; e < graph.edges.size(); e++) {
		if (IsSimpleRepeat(graph,e)) {
			std::cout << e << " " << graph.edges[e].length 
								<< " " << graph.vertices[graph.edges[e].src].InDegree()
								<< " " << graph.vertices[graph.edges[e].dest].OutDegree() 
								<< std::endl;
		}
	}
	return 0;
}


ssize_t IsSimpleRepeat(IntervalGraph &g, ssize_t e) {
	ssize_t inEdgeIndex, inEdge;
	ssize_t outEdgeIndex, outEdge;
	ssize_t srcVertex, destVertex;

	srcVertex = g.edges[e].src;
	destVertex = g.edges[e].dest;
	
	ssize_t isSimple = 1;
	for(inEdgeIndex = g.vertices[srcVertex].FirstIn();
			inEdgeIndex != g.vertices[srcVertex].EndIn() and isSimple;
			inEdgeIndex = g.vertices[srcVertex].NextIn(inEdgeIndex)) {
		inEdge = g.vertices[srcVertex].in[inEdgeIndex];
		if (g.vertices[g.edges[inEdge].src].OutDegree() != 1)
			isSimple = 0;
	}

	for(outEdgeIndex = g.vertices[destVertex].FirstOut();
			outEdgeIndex != g.vertices[destVertex].EndOut() and isSimple;
			outEdgeIndex = g.vertices[destVertex].NextOut(outEdgeIndex)) {
		outEdge = g.vertices[destVertex].out[outEdgeIndex];
		if (g.vertices[g.edges[outEdge].src].InDegree() != 1)
			isSimple = 0;
	}

	return isSimple and 
		g.vertices[srcVertex].InDegree() == g.vertices[destVertex].OutDegree();
}

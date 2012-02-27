/***************************************************************************
 * Title:          PrintDisjointEdges.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "IntervalGraph.h"
#include "IntegralTupleStatic.h"

int main(int argc, char* argv[]) {

	std::string graphName = argv[1];
	IntervalGraph graph;
	int vertexSize;
	ReadIntervalGraph(graphName, graph, vertexSize);

	ssize_t e;
	for (e = 0; e < graph.edges.size(); e++ ){ 
		ssize_t src = graph.edges[e].src;
		if (graph.vertices[src].InDegree() != 0 and
				graph.IsEdgeInDisjoint(e,2)) {
			std::cout << e << " " << graph.edges[e].length << " " 
								<< graph.vertices[graph.edges[e].src].InDegree() << " "
								<< graph.vertices[graph.edges[e].dest].OutDegree() << std::endl;
		}
	}
}

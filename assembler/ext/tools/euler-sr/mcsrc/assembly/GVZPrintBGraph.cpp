/***************************************************************************
 * Title:          GVZPrintBGraph.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "DeBruijnGraph.h"
#include "BVertex.h"
#include "BEdge.h"


int main(int argc, char* argv[]) {
	if (argc != 3) {
		std::cout << "usage: gvzPrintBGraph bgraph gvzgraph" << std::endl;
		exit(0);
	}

	std::string bGraphName, gvzGraphName;
	bGraphName = argv[1];
	gvzGraphName = argv[2];
	BVertexList vertices;
	BEdgeList edges;

	ReadBGraph(bGraphName, vertices, edges);
	GVZPrintBGraph(vertices, edges, gvzGraphName);
	return 0;
}
	

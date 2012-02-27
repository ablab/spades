/***************************************************************************
 * Title:          PrintRedRepeatEdges.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  04/13/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <string>
#include "IntervalGraph.h"


int main(int argc, char* argv[]) {

	string graphFileName, gvzGraphFileName;

	if (argc < 3) {
		cout << "usage: printRedRepeatEdges graphBase graphvizFileName" << endl;
		exit(0);
	}

	graphFileName = argv[1];
	gvzGraphFileName = argv[2];

	int vertexSize;
	IntervalGraph g;
	ReadIntervalGraph(graphFileName, g, vertexSize);
	ssize_t e;
	for (e = 0; e < g.edges.size(); e++) {
		if (g.edges[e].intervals->size() > 1) {
			g.edges[e].flagged = GraphEdge::Marked;
		}
	}

	GVZPrintBGraph(g.vertices, g.edges, gvzGraphFileName);

	return 0;
}

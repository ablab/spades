/***************************************************************************
 * Title:          PrintEdgeSpanningReads.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "IntervalGraph.h"
#include "SeqReader.h"
#include "ReadPaths.h"


int main(int argc, char* argv[]) {

	std::string graphFileName, readsName, readsOut;
	graphFileName = argv[1];
	int vertexSize = atoi(argv[2]);
	readsName = argv[3];
	readsOut  = argv[4];

	IntervalGraph graph;
	ReadIntervalGraph(graphFileName, graph, vertexSize);

	ssize_t p;
	for (p = 0; p < graph.paths.size(); p++) {
		if (graph.pathLengths[p] > 1) {


		}
	}

}

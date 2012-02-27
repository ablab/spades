/***************************************************************************
 * Title:          PrintPath.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  10/24/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "IntervalGraph.h"
#include "DeBruijnGraph.h"
#include "SimpleSequence.h"

#include "IntegralTupleStatic.h"


int main(int argc, const char* argv[]) {

	std::string graphName;
	IntervalGraph graph;

	if (argc < 2) {
		std::cout << "printPaths  A program to print paths corresponding to reads in the graph.\n\n";
		std::cout << " usage: printPaths grpahBase path1 [path2...] [-vertexSize v]\n";
		std::cout << " vertexSize v  (20) The graph was constructed using a vertex size of 'v'" 
							<< std::endl;
		exit(1);
	}
	graphName = argv[1];
	std::vector<ssize_t> pathsToPrint;
	int argi;
	argi = 2;
	int vertexSize = 20;
	while (argi < argc) {
		if (strcmp(argv[argi], "-vertexSize") == 0) {
			vertexSize = atoi(argv[++argi]);
		}
		else {
			pathsToPrint.push_back(atoi(argv[argi]));
		}
		argi++;
	}

	ReadIntervalGraph(graphName, graph, vertexSize);


	ssize_t p;
	for (p = 0; p < pathsToPrint.size(); p++) {
		graph.PrintPath(pathsToPrint[p], std::cout);
		std::cout << std::endl;
		graph.PrintPathReverse(pathsToPrint[p] + 1, std::cout);
		std::cout << std::endl;
	}
	return 0;
}

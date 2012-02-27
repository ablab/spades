/***************************************************************************
 * Title:          ReorderIntervals.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  10/24/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "DeBruijnGraph.h"
#include "IntervalGraph.h"
#include "IntegralTupleStatic.h"


int main(int argc, char *argv[]) {
	
	std::string baseInName, baseOutName;

	int argi = 1;
	if (argc < 4) {
		std::cout << "usage: reorderIntervals baseInName baseOutName numReads [-vertexSize v]"
							<< std::endl;
		std::cout << "Transforms the read1,comp1 ... readN,compN paired ordering"
							<< " to read1...readN " << std::endl
							<< "format used by euler_et. " << std::endl;
		exit(0);
	}
	baseInName = argv[argi++];
	baseOutName = argv[argi++];
	

	std::string intvFileName, bgraphFileName, pathFileName;
	ssize_t numReads;
	std::string intvOutName, bGraphOutName;
  int vertexSize = 20;

	numReads = atosz(argv[argi++]);
	while (argi < argc) {
		if (strcmp(argv[argi], "-vertexSize") == 0) {
			vertexSize = atoi(argv[++argi]);
		}
		++argi;
	}
	intvFileName = baseInName + ".intv";
	bgraphFileName = baseInName + ".bgraph";
	pathFileName = baseInName + ".path";
	intvOutName  = baseOutName + ".intv";
	bGraphOutName = baseOutName + ".bgraph";


	ssize_t e, i;
	IntervalGraph graph;
  graph.vertexSize = vertexSize;

	graph.ReadIntervalGraph(bgraphFileName, intvFileName, pathFileName);

	for (e = 0; e < graph.edges.size(); e++) {
		for (i = 0; i < graph.edges[e].intervals->size(); i++) {
			if (graph.IsForwardRead((*graph.edges[e].intervals)[i].read)) {
				(*graph.edges[e].intervals)[i].read /= 2;
			}
			else {
				(*graph.edges[e].intervals)[i].read = numReads + ((*graph.edges[e].intervals)[i].read-1)/2;
			}
		}
	}

	graph.PrintIntervalGraph(bGraphOutName, intvOutName);



}

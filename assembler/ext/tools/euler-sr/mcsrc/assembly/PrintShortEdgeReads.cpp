/***************************************************************************
 * Title:          PrintShortEdgeReads.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  09/06/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "IntervalGraph.h"
#include "DNASequence.h"
#include "IntegralTupleStatic.h"

int main(int argc, char* argv[]) {

	std::string graphBaseName, readsName, shortEdgeReadsName;
	ssize_t shortEdgeLength;

	std::set<ssize_t> shortEdgeIndices;

	int argi = 1;
	if (argc < 5) {
		std::cout << "usage: printShortEdgeReads graph reads length shortEdgeReads"
							<<std::endl;
		exit(1);
	}
	graphBaseName = argv[argi++];
	readsName     = argv[argi++];
	shortEdgeLength = atoi(argv[argi++]);
	shortEdgeReadsName = argv[argi++];

	
	IntervalGraph graph;
	std::ifstream readsIn;
	std::ofstream readsOut;
	openck(readsName, readsIn, std::ios::in);
	openck(shortEdgeReadsName, readsOut, std::ios::out);
	int vertexSize;
	ReadIntervalGraph(graphBaseName, graph, vertexSize);
	
	ssize_t e, i;
	for (e = 0; e < graph.edges.size(); e++) {
		if (graph.edges[e].length <= shortEdgeLength) {
			for (i = 0; i < graph.edges[e].intervals->size();i++) {
				shortEdgeIndices.insert((*graph.edges[e].intervals)[i].read/2);
			}
		}
	}

	ssize_t printSeq = 0;
	std::string line;
	ssize_t readIndex = 0;
	while(readsIn) {
		if (!std::getline(readsIn, line)) {
			break;
		}
		if (line.size() > 0 and line[0] == '>') {
			if (shortEdgeIndices.find(readIndex) !=
					shortEdgeIndices.end()) {
				printSeq = 1;
				readsOut << line << std::endl;
			}
			else {
				printSeq = 0;
			}
			++readIndex;
		}
		else {
			if (printSeq) {
				readsOut << line << std::endl;
			}
		}
	}
	return 0;
}



/***************************************************************************
 * Title:          GenerateEdgeCoveringReads.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "IntervalGraph.h"
#include "DeBruijnGraph.h"
#include "RepeatSearch.h"
#include <string>
#include <fstream>
#include "utils.h"


void PrintUsage() {
	std::cout << "usage: generateEdgeCoveringReads graphBase graphVertexSize edgeCoverSize readsOut" << std::endl;
}

int main (int argc, char* argv[]) {

	std::string graphBase, readsOutName;
	int vertexSize, edgeCoverSize;
	if (argc < 5) {
		std::cout << "usage: generateEdgeCoveringReads graphBase vertexSize edgeCoverSize readsOutName" << std::endl;
		exit(1);
	}
	int argi = 1;
	graphBase     = argv[argi++];
	vertexSize    = atoi(argv[argi++]);
	edgeCoverSize = atoi(argv[argi++]);
	readsOutName  = argv[argi++];
	
	IntervalGraph graph;
	ReadIntervalGraph(graphBase, graph, vertexSize);
	ssize_t e;
	std::ofstream readsOut;
	openck(readsOutName, readsOut, std::ios::out);
	DNASequence coverRead;
	ssize_t p;
	coverRead.length = edgeCoverSize;
	std::stringstream title;
	for (e = 0; e < graph.edges.size(); e++) {
		if (graph.edges[e].length > edgeCoverSize) {
			for (p = 0; p < graph.edges[e].length - edgeCoverSize + 1; p++ ){
				coverRead.seq = &(graph.edges[e].seq.seq[p]);
				title.str("");
				title << e << "_" << p;
				coverRead.namestr = title.str();
				coverRead.PrintlnSeq(readsOut);
			}
		}
	}
	readsOut.close();
}

	

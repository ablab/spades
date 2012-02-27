/***************************************************************************
 * Title:          GrowEdges.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "IntervalGraph.h"
#include "DNASequence.h"


int main(int argc, char* argv[]) {

	std::string graphBaseName, newGraphName, originalReadsName;
	IntervalGraph graph;
	if (argc < 4) {
		std::cout << "usage: printEdgeEnds graphBase originalReads newGraph" << std::endl;
		exit(1);
	}
	assert("not done");
	ssize_t edgeEndLength;
	graphBaseName     = argv[1];
	int vertexSize    = argv[2];
	originalReadsName = argv[3];
	newGraphName      = argv[4];

	ReadIntervalGraph(graphBaseName, graph, );

	std::ofstream edgeEndsOut;
	openck(edgeFileName, edgeEndsOut, std::ios::out);
	DNASequence edgeEnd;
	ssize_t e;
	std::stringstream namestrm;

	for (e = 0; e < graph.edges.size(); e++ ) {
		if (graph.vertices[graph.edges[e].dest].OutDegree() == 0 and
				graph.vertices[graph.edges[e].dest].InDegree() == 1) {
			// this is a sink edge, print the suffix
			if (graph.edges[e].length >= edgeEndLength) {
				edgeEnd.seq = &(graph.edges[e].seq.seq[graph.edges[e].length - edgeEndLength]);
				edgeEnd.length = edgeEndLength;
				namestrm.str("");
				namestrm << e << "_" << graph.edges[e].length - edgeEndLength;
				edgeEnd.namestr = namestrm.str();
				edgeEnd.PrintlnSeq(edgeEndsOut);
			}
		}
		else if (graph.vertices[graph.edges[e].src].InDegree() == 0 and
						 graph.vertices[graph.edges[e].src].OutDegree() == 1) {
			if (graph.edges[e].length >= edgeEndLength) {
				edgeEnd.seq = &(graph.edges[e].seq.seq[0]);
				edgeEnd.length = edgeEndLength;
				namestrm.str("");
				namestrm << e << "_0";
				edgeEnd.namestr = namestrm.str();
				edgeEnd.PrintlnSeq(edgeEndsOut);
			}
		}
	}
	return 0;
}
			
				
			

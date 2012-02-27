/***************************************************************************
 * Title:          MarkErroneousEdges.cpp 
 * Author:         Mark Chaisson
 * Created:        2010
 * Last modified:  02/26/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "bbbwt/BBBWTQuery.h"
#include "../IntervalGraph.h"

#include "../IntegralTupleStatic.h"

int main(int argc, char* argv[]) {

	std::string graphBase, csaFileName;
	std::string bgraphFile, edgeFile, intvFile;
	std::string graphOutName;
	std::string pathFile;
	
	if (argc != 4) {
		std::cout <<" usage: markErroneousEdges graphBase genomeCSA graphOut.dot" <<std::endl;
		exit(1);
	}
	graphBase = argv[1];
	csaFileName = argv[2];
	graphOutName = argv[3];
	BBBWT csa;
	BW::Read(csaFileName, csa);


	bgraphFile = graphBase + ".bgraph";
	intvFile   = graphBase + ".intv";
	edgeFile  = graphBase + ".edge";
	pathFile  = graphBase + ".path";
	IntervalGraph graph;
	graph.ReadIntervalGraph(bgraphFile, intvFile, pathFile);

  ReadSequences(edgeFile, graph.edges);

	ssize_t e;
	for (e = 0; e < graph.edges.size(); e++) {
		ssize_t low, high;
		BW::Query(graph.edges[e].seq, csa, low, high);
		if (high - low  == 0) {
			graph.edges[e].flagged = GraphEdge::Marked;
		}
	}

	
	GVZPrintBGraph(graph.vertices, graph.edges, graphOutName);
	return 0;
}

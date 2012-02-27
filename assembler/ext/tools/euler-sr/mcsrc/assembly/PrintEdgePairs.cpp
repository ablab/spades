/***************************************************************************
 * Title:          PrintEdgePairs.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  10/24/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "IntervalGraph.h"
#include "DNASequence.h"
#include <iostream>
#include <ios>
#include "ParseTitle.h"
#include "IntegralTupleStatic.h"

int main(int argc, char* argv[]) {

	std::string graphBaseName, pairFileName;
	IntervalGraph graph;
	if (argc < 4) {
		std::cout << "usage: printEdgePairs graphBase vertexSize minLength pairFileName " << std::endl;
		std::cout << "       Prints pairs of adjacent edges from the repeat graph." << std::endl;
		std::cout << "       The parameter minLength is the minimum length of two pairs of edges to print."
							<< std::endl;
		exit(1);
	}
	ssize_t minLength;
	int vertexSize;
	graphBaseName = argv[1];
	vertexSize = atoi(argv[2]);
	minLength  = atoi(argv[3]);
	pairFileName =  argv[4];
	ReadIntervalGraph(graphBaseName, graph, vertexSize, 0);

	std::ofstream edgePairOut;
	openck(pairFileName, edgePairOut, std::ios::out);
	DNASequence edgePair;

	//UNUSED+// ssize_t i;
	ssize_t e ;
	std::stringstream namestrm;
	//UNUSED// ssize_t curEndLength;
	ssize_t destVertex;
	ssize_t outEdge, outEdgeIndex;
	ssize_t curEdgeLength, outEdgeLength;
	
	for (e = 0; e < graph.edges.size(); e++ ) {
		destVertex = graph.edges[e].dest;
		if (graph.edges[e].seq.length > minLength) 
			curEdgeLength = minLength;
		else 
			curEdgeLength = graph.edges[e].seq.length;

		for (outEdgeIndex = graph.vertices[destVertex].FirstOut();
				 outEdgeIndex < graph.vertices[destVertex].EndOut();
				 outEdgeIndex = graph.vertices[destVertex].NextOut(outEdgeIndex)) {
			outEdge = graph.vertices[destVertex].out[outEdgeIndex];
			
			if (graph.edges[outEdge].seq.length > minLength) 
				outEdgeLength = minLength - graph.vertices[destVertex].vertexSize;
			else 
				outEdgeLength = graph.edges[outEdge].seq.length - graph.vertices[destVertex].vertexSize;
			// Don't bother printing super short edges

			ssize_t edgeLength = curEdgeLength + outEdgeLength;

			if (edgeLength < minLength) 
				continue;

			edgePair.Reset(edgeLength);
			memcpy((char*) edgePair.seq,
						 &graph.edges[e].seq.seq[graph.edges[e].seq.length - curEdgeLength],
						 curEdgeLength);

			assert (outEdgeLength > 0);
			// copy the second pair edge
			memcpy((char*) &(edgePair.seq[curEdgeLength]),
						 &(graph.edges[outEdge].seq.seq[graph.vertices[destVertex].vertexSize]),
						 outEdgeLength);
			std::stringstream namestrm;
			namestrm << e << "_" << outEdge;
			edgePair.namestr = namestrm.str();
			edgePair.PrintlnSeq(edgePairOut);
		}
	}
	return 0;
}
			
				
			

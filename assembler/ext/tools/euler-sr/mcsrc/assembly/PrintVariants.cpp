/***************************************************************************
 * Title:          PrintVariants.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  10/24/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "IntervalGraph.h"
#include <sstream>
using namespace std;
#include "IntervalGraph.h"
#include <iostream>
#include <fstream>
#include "IntegralTupleStatic.h"

void PrintUsage() {
	cout << "usage: printVariants graph variants [-minAnchorLength m]" << endl << endl;
	cout << "  -minAnchorLength m (200) Only consider variants after edges" << endl;
	cout << "                         that are at least 'm' nucleotides."  << endl;
	cout << "  -maxBranches b" << endl; // TODO: explain what this means
}


void CollectVariants(IntervalGraph &g, ssize_t edgeIndex, ssize_t anchorLength, ssize_t depth, ssize_t maxDepth,
										string curSeq,	list<ssize_t> &edgeIndices, ostream &varOut) {

	if (depth > maxDepth) {
		return;
	}
	edgeIndices.push_back(edgeIndex);
	if (depth == 0) {
		// Get the 'anchor length' suffix of the first edge.
		assert(g.edges[edgeIndex].length >= anchorLength);
		ssize_t edgeLength = g.edges[edgeIndex].length;
		curSeq.assign((const char*) g.edges[edgeIndex].seq.seq, edgeLength-anchorLength, anchorLength);
	}
	else if (g.edges[edgeIndex].length >= anchorLength) {
		if (depth > 1) {
			// This is the tail end of a sequence.
			curSeq.append((const char*) &g.edges[edgeIndex].seq.seq[g.vertexSize], anchorLength - g.vertexSize);
			stringstream titleStrm;
			list<ssize_t>::iterator indexIt;
			titleStrm << ">";
			for (indexIt = edgeIndices.begin(); indexIt != edgeIndices.end(); ++indexIt) 
				titleStrm << *indexIt << "(" << g.edges[*indexIt].length << ", " 
									<< g.edges[*indexIt].multiplicity << ") ";
			
			varOut << titleStrm.str() << endl;
			varOut << curSeq << endl;
		}
	}
	else {
		curSeq.append((const char*) &g.edges[edgeIndex].seq.seq[g.vertexSize], g.edges[edgeIndex].length - g.vertexSize);
	}

	ssize_t dest;
	ssize_t outEdge, outEdgeIndex;
	dest = g.edges[edgeIndex].dest;
	

	for (outEdgeIndex = g.vertices[dest].FirstOut();
			 outEdgeIndex != g.vertices[dest].EndOut();		 
			 outEdgeIndex = g.vertices[dest].NextOut(outEdgeIndex)) {
		outEdge = g.vertices[dest].out[outEdgeIndex];
		CollectVariants(g, outEdge, anchorLength, depth + 1, maxDepth, curSeq, edgeIndices, varOut);
	}
	edgeIndices.pop_back();
}

int main(int argc, char* argv[]) {
	string graphName, varName, edgeName;
	if (argc < 3) {
		PrintUsage();
		exit(0);
	}
	
	graphName = argv[1];
	edgeName  = graphName + ".edge";
	varName   = argv[2];
	ssize_t minAnchorLength = 200;
	ssize_t maxBranches = 2;
	int argi = 4;
	while (argi < argc){ 
		if (strcmp(argv[argi], "-minAnchorLength") == 0) {
			minAnchorLength = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-maxBranches") == 0) {
			maxBranches = atoi(argv[++argi]);
		}
		++argi;
	}

	IntervalGraph graph;
	std::ofstream variantsOut;
	openck(varName, variantsOut, std::ios::out);
	int vertexSize;
	ReadIntervalGraph(graphName, graph, vertexSize);
	ReadSequences(edgeName, graph.edges);
	ssize_t e;
	for (e = 0; e < graph.edges.size(); e++) {
		//UNUSED+// ssize_t  src;
		ssize_t dest;
		dest = graph.edges[e].dest;
		if (graph.edges[e].length < minAnchorLength)
			continue;
		if (graph.vertices[dest].OutDegree() > 1) {
			// possible variant.
			list<ssize_t> ei;
			std::cout << "collecting variants for " << e << " " << graph.edges[e].length << std::endl;
			CollectVariants(graph, e, minAnchorLength, 0, maxBranches, "", ei, variantsOut);
		}
	}

}

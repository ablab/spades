/***************************************************************************
 * Title:          LastChanceReads.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  10/03/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "compatibility.h"
#include "IntervalGraph.h"
#include "SeqReader.h"
#include "SimpleSequence.h"
#include "DNASequence.h"
#include <set>


#include "IntegralTupleStatic.h"

void StoreComponentReads(IntervalGraph &graph, ssize_t numAnchors,
												 ssize_t edge,
												 std::set<ssize_t> &readIndices);


int main(int argc, char* argv[]) {
	std::string graphBaseName, readFileName;
	ssize_t numAnchors;
	std::string lastChanceName;
	int vertexSize;

	if (argc < 6) {
		std::cout << "usage: lastChanceReads graphFile readsFile numAnchors vertexSize lastChanceName" << std::endl;
		exit(1);
	}
	int argi = 1;
	
	graphBaseName = argv[argi++];
	readFileName  = argv[argi++];
	numAnchors    = atosz(argv[argi++]);
	vertexSize    = atoi(argv[argi++]);
	lastChanceName= argv[argi++];

	IntervalGraph graph;
	ReadIntervalGraph(graphBaseName, graph, vertexSize); 

	// Open the last chance reads for writing before
	// processing all the reads so that the user will
	// know if opening failed before running some long check
	std::ofstream lastChanceReads;
	openck(lastChanceName, lastChanceReads, std::ios::out);
	
	std::ifstream readsIn;
	openck(readFileName,readsIn, std::ios::in);

	ssize_t e, i;
	std::set<ssize_t> lastChanceIndices;
	std::set<ssize_t> firstComponentIndices;
	ssize_t isFirstComponent = 1;
	for (e = 0; e < graph.edges.size(); e++ ) {
		// Don't try to process this component if it 
		// has already been traversed
		if (graph.edges[e].traversed == GraphEdge::Marked)
			continue;

		for (i = 0; i < graph.edges[e].intervals->size(); i++) {
			if ((*graph.edges[e].intervals)[i].read/2 < numAnchors) {
				// Found an anchoring read.  Output all reads 
				StoreComponentReads(graph, numAnchors, e, lastChanceIndices);
				// Done processing this component, don't look through
				// other intervals.
				if (isFirstComponent) 
					firstComponentIndices = lastChanceIndices;
				isFirstComponent = 0;
				break;
			}
		}
	}
	
	ssize_t r;
	DNASequence read;
	std::cout << "Reclaiming " << lastChanceIndices.size() << " reads." << std::endl;
	std::cout << "first component: " << firstComponentIndices.size() << std::endl;
	// Read in all the anchors
	for (r = 0; r < numAnchors; r++) {
		SeqReader::GetSeq(readsIn, read);
	}
	
	// Output all the last chance reads that were stored
	ssize_t readIndex = numAnchors;
	while (SeqReader::GetSeq(readsIn, read)) {
		if (lastChanceIndices.find(readIndex) != lastChanceIndices.end()) {
			read.PrintlnSeq(lastChanceReads);
		}
		//		if (firstComponentIndices.find(readIndex) != firstComponentIndices.end()) {
		//			read.PrintlnSeq(std::cout);
		//		}
		++readIndex;
	}
	return 0;
}

void StoreComponentReads(IntervalGraph &graph, ssize_t numAnchors,
												 ssize_t edge,
												 std::set<ssize_t> &readIndices) {
	if (graph.edges[edge].traversed == GraphEdge::Marked)
		return;

	graph.edges[edge].traversed = GraphEdge::Marked;

	ssize_t i;
	for (i = 0; i < graph.edges[edge].intervals->size(); i++) {
		if ((*graph.edges[edge].intervals)[i].read/2 >= numAnchors) {
			readIndices.insert((*graph.edges[edge].intervals)[i].read/2);
		}
	}

	ssize_t destVertex, srcVertex;
	destVertex = graph.edges[edge].dest;
	srcVertex  = graph.edges[edge].src;

	ssize_t outEdgeIndex, outEdge, inEdgeIndex, inEdge;
	for (outEdgeIndex = graph.vertices[destVertex].FirstOut();
			 outEdgeIndex != graph.vertices[destVertex].EndOut();
			 outEdgeIndex = graph.vertices[destVertex].NextOut(outEdgeIndex)) {
		outEdge = graph.vertices[destVertex].out[outEdgeIndex];
		StoreComponentReads(graph, numAnchors, outEdge, readIndices);
	}

	for (inEdgeIndex = graph.vertices[srcVertex].FirstIn();
			 inEdgeIndex != graph.vertices[srcVertex].EndIn();
			 inEdgeIndex = graph.vertices[srcVertex].NextIn(inEdgeIndex)) {
		inEdge = graph.vertices[srcVertex].in[inEdgeIndex];
		StoreComponentReads(graph, numAnchors, inEdge, readIndices);
	}
}

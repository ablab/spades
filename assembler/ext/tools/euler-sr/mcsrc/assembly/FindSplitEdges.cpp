/***************************************************************************
 * Title:          FindSplitEdges.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "ReadPos.h"
#include "SeqReader.h"
#include "IntervalGraph.h"
#include "utils.h"
#include "IntegralTupleStatic.h"

ssize_t DFSForDest(IntervalGraph g, ssize_t src, ssize_t dest, ssize_t maxSearch);

int main(int argc, char* argv[]) {
	if (argc < 5) {
		std::cout << "usage: findSplitEdges graphFile vertexSize refEdgeName minEdgeLength" << std::endl;
		exit(0);
	}
	std::string graphFileName, refEdgeFileName;
	int vertexSize;
	ssize_t minEdgeLength;
	graphFileName   = argv[1];
	vertexSize      = atoi(argv[2]);
	refEdgeFileName = argv[3];
	minEdgeLength   = atoi(argv[4]);

	IntervalGraph graph;
	ReadIntervalGraph(graphFileName, graph, vertexSize);

	// Creat an edge position database
	int tupleSize = vertexSize + 10;
	
	std::cout << "making edge database." << std::endl;
	// Make a list of sequences
	SimpleSequenceList edges;
	edges.resize(graph.edges.size());
	ssize_t e, p;
	for(e = 0; e < graph.edges.size(); e++ ) {
		edges[e].seq = graph.edges[e].seq.seq;
		edges[e].length = graph.edges[e].seq.length;
	}
	ReadPos::sequences  = &edges;
	ReadPos::hashLength = tupleSize;
	ReadPositions readPosList;

	for (e = 0; e < graph.edges.size(); e++ ) {
		for (p = 0; p < graph.edges[e].length - tupleSize + 1; p++ ) 
			readPosList.push_back(ReadPos(e, p));
	}
	std::cout << "sorting it." << std::endl;
	// Make the list query-able.
	//UNUSED// CompareReadPos comp;
	std::sort(readPosList.begin(), readPosList.end());

	DNASequence edge;
	std::ifstream in;
	openck(refEdgeFileName, in, std::ios::in);
	while (SeqReader::GetSeq(in, edge, SeqReader::noConvert)) {
		// Now look for the tuples in edge.
		if (edge.length < minEdgeLength)
			continue;

		ssize_t firstPos = 0;
		ssize_t lastPos = edge.length - tupleSize ;
		// Find the first pos
		ReadPositions::iterator firstIt, lastIt;
		CompareReadPos readPosCompStr;
		//		std::cout << "searching for match on edge of length: " << edge.length << std::endl;
		firstIt = std::lower_bound(readPosList.begin(), readPosList.end(),
															 (const char*) &(edge.seq[firstPos]), readPosCompStr );
		while (firstPos < edge.length - tupleSize + 1 and 
					 (firstIt == readPosList.end() or
						readPosCompStr(*firstIt, (const char*) &(edge.seq[firstPos])) != 0)) {
			++firstPos;
			firstIt = std::lower_bound(readPosList.begin(), readPosList.end(),
																 (const char*) &(edge.seq[firstPos]), readPosCompStr );
		}
		
		lastIt = std::lower_bound(readPosList.begin(), readPosList.end(),
															(const char*) &(edge.seq[lastPos]), readPosCompStr);
		while (lastPos > firstPos and
					 (lastIt == readPosList.end() or
						readPosCompStr(*lastIt, (const char*) &(edge.seq[lastPos])) != 0)) {
			--lastPos;
			lastIt = std::lower_bound(readPosList.begin(), readPosList.end(),
																(const char*) &(edge.seq[lastPos]), readPosCompStr );
		}
		
		if (lastPos != firstPos) {
			std::list<ssize_t> edgeStack;
			if ((*firstIt).read != (*lastIt).read) {
				//UNUSED// ssize_t maxSearch = 100;
				if (graph.vertices[graph.edges[(*firstIt).read].dest].OutDegree() > 0) {
					std::cout << "Split edge: " << edge.length << std::endl;
					std::cout << "edge: " << (*firstIt).read << " "
										<< graph.edges[(*firstIt).read].length << " "
										<< (*lastIt).read << " " 
										<< graph.edges[(*lastIt).read].length << std::endl;
					// try and find the next edge.
					edgeStack.push_back((*firstIt).read);
					ssize_t srcEdge = (*firstIt).read;
					ssize_t destEdge = (*lastIt).read;
					//UNUSED// ssize_t destFound = 0;


					DFSForDest(graph, srcEdge, destEdge, 5);
				/*				while (edgeStack.size() > 0 and maxSearch >= 0 and !destFound) {
					ssize_t outEdge, outEdgeIndex;

					edgeStack.pop_front();
					ssize_t dest = graph.edges[edgeIndex].dest;
					for (outEdgeIndex = graph.vertices[dest].FirstOut();
							 outEdgeIndex != graph.vertices[dest].EndOut();
							 outEdgeIndex = graph.vertices[dest].NextOut(outEdgeIndex)) {
						outEdge = graph.vertices[dest].out[outEdgeIndex];
						if (outEdge == destEdge) {
							destFound = 1;
							break;
						}
						edgeStack.push_back(outEdge);
						--maxSearch;
					}
				}
				if (destFound) {
					std::cout << "found a route to dest." << std::endl;
				}
				*/
				}
			}
		}
	}
}

ssize_t DFSForDest(IntervalGraph g, ssize_t src, ssize_t dest, ssize_t maxSearch) {
	if (src == dest) {
		return 1;
	}
	else {
		if (maxSearch == 0)
			return 0;
		ssize_t destVertex = g.edges[src].dest;
		ssize_t out, outIndex;
		for (outIndex = g.vertices[destVertex].FirstOut();
				 outIndex != g.vertices[destVertex].EndOut();
				 outIndex = g.vertices[destVertex].NextOut(outIndex)) {
			out = g.vertices[destVertex].out[outIndex];
			if (DFSForDest(g, out, dest, maxSearch-1)) {
				std::cout << out << " (" << g.edges[out].length << ") " << std::endl;
				return 1;
			}
		}
		return 0;
	}
}


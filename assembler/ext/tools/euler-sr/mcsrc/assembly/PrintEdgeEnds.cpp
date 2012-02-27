/***************************************************************************
 * Title:          PrintEdgeEnds.cpp 
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
#include <iostream>
#include <ios>
#include "ParseTitle.h"
#include "IntegralTupleStatic.h"

int main(int argc, char* argv[]) {

	std::string graphBaseName, edgeFileName, fixedReadsFile ;
	std::string fixedReadsName, originalReadsName, endReadsFileName;
	IntervalGraph graph;
	if (argc < 6) {
		std::cout << "usage: printEdgeEnds graphBase edgeLength fixedReads originalReads endFileName endReadsFileName" << std::endl;
		exit(1);
	}
	ssize_t edgeEndLength;
	graphBaseName = argv[1];
	edgeEndLength = atoi(argv[2]);
	fixedReadsName= argv[3];
	originalReadsName=argv[4];
	edgeFileName  = argv[5];
	endReadsFileName = argv[6];
	int vertexSize;
	ReadIntervalGraph(graphBaseName, graph, vertexSize);
	std::set<ssize_t> edgeEndFixedIndices, edgeEndOriginalIndices;

	std::ofstream edgeEndsOut;
	openck(edgeFileName, edgeEndsOut, std::ios::out);
	DNASequence edgeEnd;
	ssize_t e, i;
	std::stringstream namestrm;
	ssize_t curEndLength;
	for (e = 0; e < graph.edges.size(); e++ ) {
		if (graph.vertices[graph.edges[e].dest].OutDegree() == 0 and
				graph.vertices[graph.edges[e].dest].InDegree() == 1) {
			// this is a sink edge, print the suffix
			if (graph.edges[e].length >= edgeEndLength) {
				curEndLength = edgeEndLength;
			}
			else {
				curEndLength = graph.edges[e].length;
			}
			edgeEnd.seq = &(graph.edges[e].seq.seq[graph.edges[e].length - 
																						 curEndLength]);
			edgeEnd.length = curEndLength;
			namestrm.str("");
			namestrm << e << "_" << graph.edges[e].length - curEndLength;
			edgeEnd.namestr = namestrm.str();
			edgeEnd.PrintlnSeq(edgeEndsOut);
			
			// store all the indices of the reads that map to the end of this edge
			ssize_t endEdgeStart  = graph.edges[e].length - curEndLength;
			
			for (i = 0; i < graph.edges[e].intervals->size(); i++) { 
				if ((*graph.edges[e].intervals)[i].edgePos > endEdgeStart) {
					edgeEndFixedIndices.insert((*graph.edges[e].intervals)[i].read/2);
				}
			}
		}
		else if (graph.vertices[graph.edges[e].src].InDegree() == 0 and
						 graph.vertices[graph.edges[e].src].OutDegree() == 1) {
			if (graph.edges[e].length >= edgeEndLength) {
				curEndLength = edgeEndLength;
			}
			else {
				curEndLength = graph.edges[e].length;
			}
			edgeEnd.seq = &(graph.edges[e].seq.seq[0]);
			edgeEnd.length = curEndLength;
			namestrm.str("");
			namestrm << e << "_0";
			edgeEnd.namestr = namestrm.str();
			edgeEnd.PrintlnSeq(edgeEndsOut);
			for (i = 0; i < graph.edges[e].intervals->size(); i++) { 
				if ((*graph.edges[e].intervals)[i].edgePos < curEndLength) {
					edgeEndFixedIndices.insert((*graph.edges[e].intervals)[i].read/2);
				}
			}
		}
	}
	std::cout << "looking for " << edgeEndFixedIndices.size() << " reads." << std::endl;
	std::ifstream fixedIn; 
	openck(fixedReadsName, fixedIn, std::ios::in);
	std::string line;
	ssize_t fixedIndex = 0;
	ssize_t origIndex;
	std::set<ssize_t>::iterator it;
	while(fixedIn) {
		std::getline(fixedIn, line);
		if (line.size() > 0 and line[0] == '>') {
			if ((it = edgeEndFixedIndices.find(fixedIndex)) != edgeEndFixedIndices.end()) {
				if (ParseKeyword(line, "index", origIndex)) {
					edgeEndOriginalIndices.insert(origIndex);
				}
				edgeEndFixedIndices.erase(it);
			}
			fixedIndex++;
		}
		// otherwise the line is a read, and we don't care about it
	}
	for (it = edgeEndFixedIndices.begin(); 
			 it != edgeEndFixedIndices.end(); it++)
		std::cout << "missing: " << *it << std::endl;
	

	// Now print all the original reads.
	ssize_t printSeq = 0;
	origIndex = 0;
	std::ifstream origIn;
	std::ofstream endReadsOut;
	std::cout << "writing " << edgeEndOriginalIndices.size() << " reads." << std::endl;
	openck(originalReadsName, origIn, std::ios::in);
	openck(endReadsFileName, endReadsOut, std::ios::out);
	while (origIn) {
		std::getline(origIn, line);
		if (line.size() > 0 and line[0] == '>') {
			if ((it = edgeEndOriginalIndices.find(origIndex)) != 
					edgeEndOriginalIndices.end()) {
				printSeq = 1;
				endReadsOut << line << std::endl;
			}
			else {
				printSeq = 0;
			}
			origIndex++;
		}
		else {
			if (printSeq) {
				endReadsOut << line << std::endl;
			}
		}
	}
	return 0;
}
			
				
			

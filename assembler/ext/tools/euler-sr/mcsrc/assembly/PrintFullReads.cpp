/***************************************************************************
 * Title:          PrintFullReads.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "IntervalGraph.h"
#include "IntegralTupleStatic.h"



int main(int argc, char* argv[]) {

	if (argc < 5) {
		std::cout << "usage: printFullReads reads graph vertexSize fullReadLength contained not" << std::endl;
		exit(0);
	}

	std::string readsName = argv[1];
	std::string baseName = argv[2];
	int vertexSize = atoi(argv[3]);
	ssize_t fullReadLength  = atoi(argv[4]);
	std::string contained = argv[5];
	std::string notContained = argv[6];
	IntervalGraph graph;
	DNASequenceList reads;
	ReadDNASequences(readsName,reads);
	ReadIntervalGraph(baseName, graph, vertexSize);

	ssize_t p;
	ssize_t i;
	ssize_t edge, index;
	//UNUSED// ssize_t readIsFull = 0;
	ssize_t readLength;
	ssize_t destVertex;
	std::ofstream contFile, notContFile;
	openck(contained, contFile, std::ios::out);
	openck(notContained, notContFile, std::ios::out);
	for (p = 0; p < graph.paths.size(); p+=2) {
		readLength = 0;
		if (graph.pathLengths[p] == 0) {
			reads[p/2].PrintSeq(notContFile);
			notContFile << std::endl;
			continue;
		}
		edge  = graph.paths[p][0].edge;
		index = graph.paths[p][0].index;
		
		if ((*graph.edges[edge].intervals)[index].readPos != 0) {
			reads[p/2].PrintSeq(notContFile);
			notContFile << std::endl;
		}
		else {
			for (i =0 ; i < graph.pathLengths[p]-1; i++) {
				destVertex = graph.edges[edge].dest;
				readLength += (*graph.edges[edge].intervals)[index].length - graph.vertices[destVertex].vertexSize;
			}
			readLength += (*graph.edges[edge].intervals)[index].length;
			
			if (readLength == fullReadLength) {
				reads[p/2].PrintSeq(contFile);
				contFile << std::endl;
			}
			else {
				reads[p/2].PrintSeq(notContFile);
				notContFile << std::endl;
			}
		}
		
	}
	return 0;
}
	 
			


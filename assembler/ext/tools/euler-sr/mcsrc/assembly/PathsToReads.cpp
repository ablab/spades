/***************************************************************************
 * Title:          PathsToReads.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  11/28/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "IntervalGraph.h"
#include "DeBruijnGraph.h"
#include "RepeatSearch.h"
#include "ThreadUtils.h"
#include "utils.h"
#include <string>
#include <fstream>
#include "IntegralTupleStatic.h"




void PrintUsage() {
	std::cout << "usage: pathsToReads graphBase graphVertexSize readsOut [-onlyMultiEdge]" << std::endl;
}


int main (int argc, char* argv[]) {

	std::string graphBase, readsOutName;
	//UNUSED+// int edgeCoverSize;
	int vertexSize ;
	if (argc < 2) {
		PrintUsage();
		exit(1);
	}
	int argi = 1;
	graphBase    = argv[argi++];
	readsOutName = argv[argi++];
	ssize_t onlyMultiEdge = 0;
	while (argi < argc) {
		if (strcmp(argv[argi], "-onlyMultiEdge") == 0) {
			onlyMultiEdge = 1;
		}
		++argi;
	}
	
	IntervalGraph graph;
	ThreadPath path;
	std::ofstream readsOut;
	openck(readsOutName, readsOut, std::ios::out);
	
	ssize_t p, i;
	ssize_t pathEdge, pathInterval;
	std::string pathSeq;
	DNASequence pathRead;
	std::stringstream titlestrm;
	ssize_t intvLength, intvPos;
	ssize_t destVertex, destVertexLength;
	//UNUSED+// ssize_t  pathLength;
	ssize_t readLength;

	ReadIntervalGraph(graphBase, graph, vertexSize);
	ssize_t minPathLength = 0;
	// adjust for multi edge reads
	if (onlyMultiEdge )
		minPathLength = 2;

	for (p = 0; p < graph.paths.size(); p+=2 ) {
		path.clear();
		if (graph.pathLengths[p] >= minPathLength) {
			readLength = 0;
			cout << "pl: " << graph.pathLengths[p] << endl;
			for (i = 0; i < graph.pathLengths[p] - 1; i++) {
				// Look up the edge, interval for this path interval
				pathEdge     = graph.paths[p][i].edge;
				pathInterval = graph.paths[p][i].index;
				
				// Find the length along this path interval
				destVertex       = graph.edges[pathEdge].dest;
				destVertexLength = graph.vertices[destVertex].vertexSize;
				intvLength       = (*graph.edges[pathEdge].intervals)[pathInterval].length -
					destVertexLength;
				intvPos = (*graph.edges[pathEdge].intervals)[pathInterval].edgePos;
				path.push_back(ThreadPathInterval(pathEdge, intvLength, intvPos));
				readLength += intvLength;
				std::cout << " " << intvLength;
			}
			// The last interval includes the entire vertex length
			pathEdge     = graph.paths[p][i].edge;
			pathInterval = graph.paths[p][i].index;
			intvLength   = (*graph.edges[pathEdge].intervals)[pathInterval].length;
			intvPos      = (*graph.edges[pathEdge].intervals)[pathInterval].edgePos;
			//			std::cout << " " << intvLength << std::endl;
			readLength += intvLength;
			//			std::cout << graph.pathLengths[p] << " " << readLength << " ";
			path.push_back(ThreadPathInterval(pathEdge, intvLength, intvPos));
			
			ThreadToSeq(graph, path, pathSeq);
			titlestrm.str("");
			titlestrm << pathEdge << "_" << intvLength << "_" << intvPos;
			pathRead.namestr = titlestrm.str();
			pathRead.Reset(pathSeq.size());
			memcpy(pathRead.seq, pathSeq.c_str(), pathSeq.size());
			pathRead.PrintlnSeq(readsOut);
		}
	}
	return 0;
}

/***************************************************************************
 * Title:          PrintMateLengthDistribution.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "IntervalGraph.h"
#include "MateLibrary.h"
#include "MateTable.h"
#include "IntegralTupleStatic.h"

void PrintUsage() {
	std::cout << "usage: printMateLengthDistribution graphName mateTableName" << std::endl;
}
int main(int argc, char* argv[]) {

	std::string graphFileName, mateTableName;
	if (argc < 3) {
		PrintUsage();
		exit(0);
	}
	graphFileName = argv[1];
	mateTableName = argv[2];

	IntervalGraph graph;
	int vertexSize;
	ReadIntervalGraph(graphFileName, graph, vertexSize);

	ReadMateList  mateTable;
	ReadMateTable(mateTableName, mateTable);

	ssize_t e, i;
	std::map<ssize_t,ssize_t> lengthCount;
	TEdgeList &edges = graph.edges;
	PathLengthList &pathLengths = graph.pathLengths;
	PathIntervalList &paths = graph.paths;

	for (e = 0; e < edges.size(); e++ ){
		for (i = 0; i < (*edges[e].intervals).size(); i++) {
			ssize_t pathPos;
			ssize_t path;
			path = (*edges[e].intervals)[i].read;
			pathPos = (*edges[e].intervals)[i].pathPos;

			if (path % 2 == 0) {
				// This is the foward read.
				ssize_t readIndex = path / 2;
				ssize_t mateIndex = mateTable[readIndex].mateIndex;
				ssize_t mateReadIndex = mateIndex * 2 + 1;

				if (pathLengths[path] == 0 or pathLengths[mateReadIndex] == 0) {
					continue;
				}
			
				ssize_t pathLength = pathLengths[path];
				ssize_t readEndEdge = paths[path][pathLength-1].edge;
				ssize_t readEndIntv = paths[path][pathLength-1].index;
				ssize_t readEndPos  = (*edges[readEndEdge].intervals)[readEndIntv].edgePos + 
					(*edges[readEndEdge].intervals)[readEndIntv].length;
			
				//UNUSED// ssize_t readEndEdgeLength = edges[readEndEdge].length;
			
				ssize_t mateBeginEdge = paths[mateReadIndex][0].edge;
				ssize_t mateBeginIntv = paths[mateReadIndex][0].index;
				ssize_t mateBeginPos  = (*edges[mateBeginEdge].intervals)[mateBeginIntv].edgePos;

				if (readEndEdge == mateBeginEdge ) {
					lengthCount[mateBeginPos - readEndPos]++;
				}
			}
		}
	}

	std::map<ssize_t, ssize_t>::iterator lenIt;
	for (lenIt = lengthCount.begin(); lenIt != lengthCount.end(); ++lenIt) {
		std::cout << (*lenIt).first << " " << (*lenIt).second << std::endl;
	}
}

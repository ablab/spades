/***************************************************************************
 * Title:          CountMatesOnSameEdge.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "IntervalGraph.h"
#include "PathLib.h"
#include "MateLibrary.h"
#include "RuleList.h"
#include "IntegralTupleStatic.h"

using namespace std;
int main(int argc, char* argv[]) {


	string graphFileName;
	
	string mateTableName;
	graphFileName = argv[1];
	mateTableName = argv[2];
	int vertexSize;
	IntervalGraph graph;
	ReadIntervalGraph(graphFileName, graph, vertexSize);
	ReadMateList mateList;
	ReadMateTable(mateTableName, mateList);
	
	ssize_t p;
	ssize_t numSameEdge = 0;
	ssize_t numSingle = 0;
	for (p = 0; p < mateList.size(); p++) {
		if (mateList[p].mateIndex == -1) {
			numSingle++;
			continue;
		}
		ssize_t pathIndex = p*2;
		ssize_t matePathIndex = mateList[p].mateIndex * 2 + 1;
		ssize_t pathLen = graph.pathLengths[pathIndex];
		ssize_t endPathEdge, beginMatePathEdge;
		endPathEdge = beginMatePathEdge = -1;
		if (pathLen == 0 and graph.pathLengths[matePathIndex] == 0) 
			continue;
		
		if (pathLen == 0 or graph.pathLengths[matePathIndex] == 0) 
			numSingle++;

		if (pathLen > 0)
			endPathEdge = graph.paths[pathIndex][pathLen-1].edge;
		if (graph.pathLengths[matePathIndex] > 0)
			beginMatePathEdge = graph.paths[matePathIndex][0].edge;
		
		if (endPathEdge != -1 and beginMatePathEdge != -1 and
				endPathEdge == beginMatePathEdge) {
			// the two mates are on the same edge... maybe they should have been marked?
			numSameEdge++;
		}
	}
	cout << numSameEdge << " mates are on the same edge." << endl;
	cout << numSingle << " mates are missing a pair." << endl;

}

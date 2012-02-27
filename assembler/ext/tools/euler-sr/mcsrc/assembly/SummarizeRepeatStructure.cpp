/***************************************************************************
 * Title:          SummarizeRepeatStructure.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  11/26/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "IntervalGraph.h"
#include <string>

using namespace std;


ssize_t GetUniqueReadCoverage(IntervalGraph &g, ssize_t e) {
	ReadIntervalList *intv;
	ssize_t numIntv;
	intv = g.edges[e].intervals;
	numIntv = intv->size();

	ssize_t readCoverage = 0;
	ssize_t i = 0;
	while (i < numIntv - 1) {
		if ((*intv)[i].read != (*intv)[i+1].read) 
			++readCoverage;
		++i;
	}
	return readCoverage;
}


int main(int argc, char *argv[]) {

	string graphFileName = argv[1];
	
	IntervalGraph graph;
	int vertexSize;

	ReadIntervalGraph(graphFileName, graph, vertexSize);
	graph.SortAllEdgeIntervalsByReadPos();

  ssize_t e;
	vector<double> avgMultiplicity;
	avgMultiplicity.resize(graph.edges.size());
	double maxMult = 0;
	ssize_t maxMultIndex = -1;
	for (e = 0; e < graph.edges.size(); e++ ){
		avgMultiplicity[e] = (((1.0)*GetUniqueReadCoverage(graph,e))/ graph.edges[e].length);
		if (avgMultiplicity[e] > maxMult) {
			maxMult = avgMultiplicity[e];
			maxMultIndex = e;
		}
		cout << graph.edges[e].length << " " << avgMultiplicity[e] << endl;
	}

	return 0;
}
	

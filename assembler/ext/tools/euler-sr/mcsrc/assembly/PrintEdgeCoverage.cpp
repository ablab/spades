/***************************************************************************
 * Title:          PrintEdgeCoverage.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "IntervalGraph.h"
#include <iostream>
#include <fstream>
#include "IntegralTupleStatic.h"

using namespace std;

void PrintUsage() {
	cout << "usage: printEdgeCoverage graph edgeIndex coverageFile [-matlab]"<<endl << endl;
	cout << "           edgeIndex - A file of indices of edges in the graph.  It is "<<endl
			 << "                       of any format really, just a list of numbers."<<endl
			 << "           coverageFile- The output file." << endl
			 << "       [-matlab]       Print commands to generate plots using a batch matlab script."<<endl;
}

int main(int argc, char* argv[]) {
	string graphName;
	string edgeIndexFileName;
	string coverageFile;
	if (argc < 4) {
		PrintUsage();
		exit(0);
	}

	graphName = argv[1];
	edgeIndexFileName = argv[2];
	coverageFile = argv[3];
	int argi = 4;
	ssize_t printMatlab = 0;
	while (argi < argc) {
		if (strcmp(argv[argi], "-matlab") == 0) {
			printMatlab = 1;
		}
		else {
			PrintUsage();
			cout << "bad option: " << argv[argi] << endl;
			exit(0);
		}
		++argi;
	}
	IntervalGraph g;
	std::ifstream edgeIndexIn;
	std::ofstream coverageOut;
	openck(edgeIndexFileName, edgeIndexIn, std::ios::in);
	openck(coverageFile, coverageOut, std::ios::out);
	int vertexSize;
	ReadIntervalGraph(graphName, g, vertexSize);
	
	ssize_t edgeIndex = 0;
	while(edgeIndexIn) {
		if (!(edgeIndexIn >> edgeIndex)) break;
		
		vector<ssize_t> edgeCov;
		edgeCov.resize(g.edges[edgeIndex].length);
		fill(edgeCov.begin(), edgeCov.end(), 0);
		
		ssize_t i;
		ssize_t readPos, edgePos;
		for (i = 0; i < (*g.edges[edgeIndex].intervals).size(); i++) {
			edgePos = (*g.edges[edgeIndex].intervals)[i].edgePos;
			for (readPos = 0; readPos < (*g.edges[edgeIndex].intervals)[i].length; ++readPos) {
				edgeCov[edgePos+readPos]++;
			}
		}
		for (edgePos = 0; edgePos < edgeCov.size(); edgePos++){
			coverageOut << edgeCov[edgePos] << endl;
		}
	}
	coverageOut.close();
}


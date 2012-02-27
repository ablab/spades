/***************************************************************************
 * Title:          FilterValidMatePaths.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  10/24/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "IntervalGraph.h"
#include "MateLibrary.h"
#include "utils.h"
#include "IntegralTupleStatic.h"


void PrintUsage() {
	std::cout << "usage: filterValidMatePaths graphName vertexSize mateTable ruleFile"
						<< std::endl;
}

int main(int argc, char*argv[]) {
	std::string graphName;
	std::string mateTableName;
	std::string ruleFileName;
	int vertexSize;
	if (argc < 4) {
		PrintUsage();
		exit(1);
	}
	graphName = argv[1];
	vertexSize = atoi(argv[2]);
	mateTableName = argv[3];
	ruleFileName = argv[4];

	IntervalGraph graph;
	std::cout << "reading the interval graph." << std::endl;
	ReadIntervalGraph(graphName, graph, vertexSize);
	
	RuleList rules;
	ParseRuleFile(ruleFileName, rules);

	ReadMateList mateList;
	ReadMateTable(mateTableName, mateList);
	std::map<ssize_t,ssize_t> matePathCount;
	std::map<ssize_t,ssize_t> mateLengths;
	ssize_t p;
	ssize_t readIndex, mateIndex;
	ssize_t mateType;
	ssize_t lastPathIntv;
	ssize_t lastEdgeIntv, mateEdgeIntv;
	ssize_t lastEdge, mateEdge, lastEdgePos, mateStartEdgePos; 
	ssize_t mp;
	std::cout << "done reading mate table." << std::endl;
	std::cout << "Here is a sample:"<< std::endl;
	ssize_t i;
	for (i = 0; i< 10; i++) {
		std::cout << mateList[i].mateIndex << " " << mateList[i].mateType << std::endl;
	}
	for (p = 0; p < graph.paths.size(); p+=2 ) {
		readIndex = p/2;
		if (mateList[readIndex].mateIndex == -1)
			continue;

		if (graph.pathLengths[p] == 0 )
			continue;

		mateIndex = mateList[readIndex].mateIndex;
		mateType  = mateList[readIndex].mateType;

		lastPathIntv = graph.pathLengths[p] - 1;
		
		lastEdge = graph.paths[p][lastPathIntv].edge;
		lastEdgeIntv = graph.paths[p][lastPathIntv].index;

		lastEdgePos = (*graph.edges[lastEdge].intervals)[lastEdgeIntv].edgePos + 
			(*graph.edges[lastEdge].intervals)[lastEdgeIntv].length;

		// use the reverse complment path of the mate since
		// the mate is sequenced in the opposite direction
		mp = mateIndex * 2 + 1;
		
		if (graph.pathLengths[mp] == 0) 
			continue;

		mateEdge     = graph.paths[mp][0].edge;
		mateEdgeIntv = graph.paths[mp][0].index;

		mateStartEdgePos = (*graph.edges[mateEdge].intervals)[mateEdgeIntv].edgePos +
			(*graph.edges[mateEdge].intervals)[mateEdgeIntv].length;

		//UNUSED// ssize_t foundAValidPath;
		std::list<ssize_t> path;
		/*
		foundAValidPath = SearchValidMatePath(graph,
																			 lastEdge, lastEdgePos,
																			 mateEdge, mateStartEdgePos,
																			 0, 
																			 rules[mateType].cloneLength - 
																			 rules[mateType].cloneVar,
																			 rules[mateType].cloneLength + 
																					rules[mateType].cloneVar, rules[mateType].cloneLength, path);
		*/
		//		if (foundAValidPath) {
			// This could take much more time.
			ssize_t numValidPaths = 0;
			ssize_t foundValidPaths = 0;
			MatePathList matePath;
			ssize_t storeMatePath = 1;
			ssize_t totalPathLength = 0;
			map<ssize_t,ssize_t> visited;
			foundValidPaths = CountValidMatePaths(graph,
																						lastEdge, lastEdgePos, // cur edge pos
																						mateEdge, mateStartEdgePos,
																						0, 
																						rules[mateType].cloneLength - rules[mateType].cloneVar,
																						rules[mateType].cloneLength + rules[mateType].cloneVar, 
																						//																						rules[mateType].cloneLength,
																						30,
																						4, numValidPaths,
																						storeMatePath, matePath, 1, totalPathLength, visited);
			if (numValidPaths >= 1) {
				//				std::cout << "mate pair: " << readIndex << " has " << numValidPaths << " paths." << std::endl;
				if (matePath.size() == 0) {
					ssize_t mateLength = mateStartEdgePos - lastEdgePos;
					if (mateLength > 0) {
						if (mateLengths.find(mateLength) == mateLengths.end())
							mateLengths[mateLength] = 1;
						else
							mateLengths[mateLength]++;
					}
				}
			}
			//		}
			
		if (numValidPaths == 0)  {
			std::cout << "mate pair: " << readIndex << " is bad." << std::endl;
			std::cout << lastEdge << "." 
								<< lastEdgePos << "(" << graph.edges[lastEdge].length << ") -> " 
								<< mateEdge << "." << mateStartEdgePos 
								<< " (" << graph.edges[mateEdge].length << ")" << std::endl ;
			ssize_t destVertex;
			destVertex = graph.edges[lastEdge].dest;
			std::cout << "dest degree: " << graph.vertices[destVertex].OutDegree() << std::endl;
		}
		if (lastEdge != mateEdge) {
			/*			std::cout << "found edge spanning pair:"  << lastEdge << "." 
								<< lastEdgePos << "(" << graph.edges[lastEdge].length << ") -> " 
								<< mateEdge << "." << mateStartEdgePos << std::endl ;
			std::list<ssize_t>::iterator pathIt;
			for (pathIt = validPath.begin(); pathIt!= validPath.end(); ++pathIt) {
				std::cout << *pathIt << "(" << graph.edges[*pathIt].length << ") ";
			}
			std::cout << std::endl;
			*/
			if (matePathCount.find(numValidPaths) == matePathCount.end()) 
				matePathCount[numValidPaths] = 1;
			else
				matePathCount[numValidPaths]++;
		}
	}
	std::map<ssize_t,ssize_t>::iterator mapIt;
	std::cout << "map histogram: " << std::endl;
	for (mapIt = matePathCount.begin(); mapIt != matePathCount.end(); ++mapIt) {
		std::cout << mapIt->first << " " << mapIt->second << std::endl;
	}
	std::cout << " mate lengths: " << std::endl;
	for (mapIt = mateLengths.begin(); mapIt != mateLengths.end(); ++mapIt) {
		std::cout << mapIt->first << " " << mapIt->second << std::endl;
	}
	
}

/***************************************************************************
 * Title:          MateLibrary.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "MateLibrary.h"
#include "utils.h"
#include "ParseTitle.h"
#include <fstream>
#include <iostream>
#include <ext/functional>
#include <iterator>
#include "IntervalGraph.h"
#include <regex.h>
#include "compatibility.h"

using namespace std;

#define MAX_EDGE_TRAVERSALS 3

void ReadMateTable(std::string &mateTableName,
									 ReadMateList &matePairs,
									 std::ostream &report) {
	std::ifstream tableIn;
	openck(mateTableName, tableIn, std::ios::in, report);
	ReadMate mate;
	ssize_t numPaired = 0;
	while(tableIn) {
		if (!(tableIn >> mate.mateIndex >> mate.mateType)) {
			break;
		}
		if (mate.mateIndex != -1) ++numPaired;
		matePairs.push_back(mate);
	}
	cout << "read: " << numPaired << " mate pairs." << endl;
}



// This is a mroe 'lightweight' path search utility that simply marks 
// edges that are on a valid mate path from 
void MarkEdgesOnValidMatePaths(IntervalGraph &g,
															 ssize_t curEdge, ssize_t curEdgePos,
															 ssize_t endEdge, ssize_t endEdgePos,
															 ssize_t curLength, ssize_t minLength, ssize_t maxLength,
															 ssize_t maxDepth,	std::vector<ssize_t> &path, ssize_t curDepth, 
															 std::vector<ssize_t> &edgeOnPathCount, 
															 ssize_t &numPaths, ssize_t &totalLength, ssize_t &totalEdges,
															 std::vector<ssize_t> &pathEdgeCount) {
	if (curDepth > maxDepth)
		return;
	
	++pathEdgeCount[curEdge];

	// Look to see if this is the last edge on the path
	if (curEdge == endEdge && curEdgePos < endEdgePos) {
		ssize_t curEdgeLength = curLength + (endEdgePos - curEdgePos + 1);
		
		if ( curEdgeLength >= minLength and
				 curEdgeLength <= maxLength) {
			ssize_t e;
			// path[0] is the first edge, so we don't 
			// need to mark that as being part of a path since
			// we are looking for paths between edges.
			for (e = 1; e < curDepth - 1; e++) {
				edgeOnPathCount[path[e]]++;
			}
			numPaths++;
			totalLength += curEdgeLength;
			totalEdges  += curDepth;
		}
	}
	// Look to see if this is the last edge on the mate path.
	if (curLength <= maxLength) {
		// Add the current edge to the path.
		path[curDepth] = curEdge;
		ssize_t destVertex;
		destVertex = g.edges[curEdge].dest;
		// advance past this edge.
		curLength += (g.edges[curEdge].length - curEdgePos);
		ssize_t destVertexLength;

		destVertexLength = g.vertices[g.edges[curEdge].dest].vertexSize;

		ssize_t outEdge, outEdgeIndex;
		// search for more valid mate paths
		for(outEdgeIndex = g.vertices[destVertex].FirstOut();
				outEdgeIndex != g.vertices[destVertex].EndOut();
				outEdgeIndex = g.vertices[destVertex].NextOut(outEdgeIndex)) {
			outEdge = g.vertices[destVertex].out[outEdgeIndex];
			destVertexLength = g.vertices[destVertex].vertexSize;
			if (pathEdgeCount[outEdge] < 5) {
				MarkEdgesOnValidMatePaths(g, outEdge, destVertexLength,
																	endEdge, endEdgePos, 
																	curLength, minLength, maxLength,	
																	maxDepth, path, curDepth + 1, edgeOnPathCount, 
																	numPaths, totalLength, totalEdges, pathEdgeCount);
			}
		}
	}
	--pathEdgeCount[curEdge];

	// returns 0 if no valid paths are found.
}

ssize_t CountValidMatePaths(IntervalGraph &g,
												ssize_t curEdge, ssize_t curEdgePos,
												ssize_t endEdge, ssize_t endEdgePos,
												ssize_t curLength,
												ssize_t minLength, ssize_t maxLength,
												ssize_t maxDepth, ssize_t maxValid, ssize_t &numValid, 
												ssize_t &storePath,
												MatePathList &matePath,
												ssize_t print, ssize_t &totalPathLength, map<ssize_t,ssize_t> &visited) {
	if (maxDepth <= 0)
		return 0;

	if (maxValid <= 0) 
		return 0;

	if (numValid >= maxValid)
		return 1;

	ssize_t foundValidPath = 0;

	// Look to see if this is the last edge on the path
	if (curEdge == endEdge && curEdgePos < endEdgePos) {
		ssize_t curEdgeLength = curLength + (endEdgePos - curEdgePos + 1);
		//		cout << "found possible mate-path of length: " << curEdgeLength << endl;
		
		if ( curEdgeLength >= minLength and
				 curEdgeLength <= maxLength) {
			numValid++;
			// We don't want to store any alternative paths, just count them 
			// from now on.
			MatePathList::iterator listIt;
		
			if (print) {
				storePath = 1;

				for (listIt = matePath.begin(); listIt != matePath.end(); ++listIt) {
					//					std::cout << (*listIt).edge << "(" << g.edges[(*listIt).edge].length << ") ";
				}
				//				std::cout << std::endl;
			}
			else {
				storePath = 0;
			}
			// Compute the total path length.
			ssize_t src;
			MatePathList::iterator endIt;
			endIt = matePath.end();
			if (matePath.size() > 1)
				--endIt;

			for (listIt = matePath.begin(); listIt != endIt; ++listIt) {
				src = g.edges[(*listIt).edge].src;
				totalPathLength += g.edges[(*listIt).edge].length - g.vertices[src].vertexSize;
			}
			return 1;
		}		
	}
	
	// Should keep searching for a path, but check to see 
	// if this part of the graph has been visited yet.

	

	// Look to see if this is the last edge on the mate path.
	if (curLength <= maxLength) {

		// Don't allow cycles when looking for unique paths since
		// they create too long of a search, and are often not unique.
		if (visited.find(g.edges[curEdge].dest) != visited.end() and
				visited[g.edges[curEdge].dest] > 1)
			return 0;
		
		//		visited.insert(g.edges[curEdge].dest);
		visited[g.edges[curEdge].dest]++;


		ssize_t outEdge, outEdgeIndex;
		ssize_t destVertex;
		destVertex = g.edges[curEdge].dest;
		// advance past this edge.
		curLength += (g.edges[curEdge].length - curEdgePos);
		//		assert(g.edges[curEdge].length - curEdgePos >= 0);
		ssize_t destVertexLength;
		if (storePath) {
			destVertexLength = g.vertices[g.edges[curEdge].dest].vertexSize;
			matePath.push_back(MatePathInterval(curEdge));
		}
		
		// search for more valid mate paths
		for(outEdgeIndex = g.vertices[destVertex].FirstOut();
				outEdgeIndex != g.vertices[destVertex].EndOut();
				outEdgeIndex = g.vertices[destVertex].NextOut(outEdgeIndex)) {
			outEdge = g.vertices[destVertex].out[outEdgeIndex];
			destVertexLength = g.vertices[destVertex].vertexSize;
			if (CountValidMatePaths(g, outEdge, destVertexLength,
															endEdge, endEdgePos, curLength,
															minLength, maxLength,	maxDepth - 1, maxValid, numValid, 
															storePath, matePath, print, totalPathLength, visited)) {
				foundValidPath = 1;
			}
		}
		if (storePath) {
			assert(matePath.size() > 0);
			matePath.pop_back();
		}
	}
	// returns 0 if no valid paths are found.
	return foundValidPath;
}												

void AssignMateOrder(ssize_t p, ssize_t mateIndex, ssize_t &mp, ssize_t &firstMate, ssize_t& secondMate) {
	if (p % 2 == 0) {
		mp = mateIndex * 2 + 1;
		firstMate = p;
		secondMate = mp;
	}
	else {
		mp = mateIndex * 2;
		firstMate = mp;
		secondMate = p;
	}
}

ssize_t StoreEdgePairMap(IntervalGraph &graph, 
										 ReadMateList  &matePairs,
										 EdgePairMap   &edgePairs,
										 ssize_t ruleType) { 
	//UNUSED+// ssize_t mp;
	ssize_t p ;
	ssize_t readIndex, mateIndex, mateType;
	for (p = 0; p < graph.paths.size(); p++ ) {
		// Find the read that is paired with read 'p', and the
		// type of mate that pairs it.
		readIndex = p/2;
		mateIndex = matePairs[readIndex].mateIndex;
		mateType  = matePairs[readIndex].mateType;
			
		// Dont try and find a mate-path if no mates 
		// are mapped to this read
		if (mateIndex == -1)
			continue;
			
		//
		// Possibly store edge pairs for a single rule type.
		//
		if (ruleType != -1 and mateType != ruleType) 
			continue;
		// Determine the path corresponding to the mate read
			
		ssize_t mp;
		ssize_t firstMate, secondMate;
		AssignMateOrder(p, mateIndex, mp, firstMate, secondMate);
		
		// If one of these paths is removed, don't try and 
		// find a mate-path.
		if (graph.pathLengths[firstMate] == 0 ||
				graph.pathLengths[secondMate] == 0) {
			continue;
		}
		
		ssize_t lastPathIntv, lastEdge, lastEdgeIntv, lastEdgePos;
			
		lastPathIntv = graph.pathLengths[firstMate] - 1;
		lastEdge = graph.paths[firstMate][lastPathIntv].edge;
		lastEdgeIntv = graph.paths[firstMate][lastPathIntv].index;
	 
		lastEdgePos = (*graph.edges[lastEdge].intervals)[lastEdgeIntv].edgePos + 
			(*graph.edges[lastEdge].intervals)[lastEdgeIntv].length;

		// use the reverse complment path of the mate since
		// the mate is sequenced in the opposite direction

		
		ssize_t mateEdge, mateEdgeIntv, mateStartEdgePos;
		mateEdge     = graph.paths[secondMate][0].edge;
		mateEdgeIntv = graph.paths[secondMate][0].index;
			
		/*
			mateStartEdgePos = (*graph.edges[mateEdge].intervals)[mateEdgeIntv].edgePos +
			(*graph.edges[mateEdge].intervals)[mateEdgeIntv].length;
		*/
		mateStartEdgePos = (*graph.edges[mateEdge].intervals)[mateEdgeIntv].edgePos;

		EdgePair edgePair;
		edgePair.edge1 = lastEdge;
		edgePair.edge2 = mateEdge;
		edgePair.mateType = mateType;
		EdgePairMap::iterator epIt;
		ssize_t read1Length = graph.CalculateReadLength(firstMate);
		ssize_t read2Length = graph.CalculateReadLength(secondMate);
		if ((epIt = edgePairs.find(edgePair)) == edgePairs.end()) {
			edgePairs[edgePair].count = 1;
			edgePairs[edgePair].meanEdge1End = lastEdgePos;
			edgePairs[edgePair].meanEdge2Start = mateStartEdgePos;
			edgePairs[edgePair].edge1 = lastEdge;
			edgePairs[edgePair].edge2 = mateEdge;
			edgePairs[edgePair].read1Length = read1Length;
			edgePairs[edgePair].read2Length = read2Length;
		}
		else {
			(*epIt).second.count++;
			(*epIt).second.meanEdge1End += lastEdgePos;
			(*epIt).second.meanEdge2Start += mateStartEdgePos;
			(*epIt).second.read1Length += read1Length;
			(*epIt).second.read2Length += read2Length;
		}
	}
	StoreMeanEdgePosition(edgePairs);
	return edgePairs.size();
}


void RemoveLowFrequencyEdgePairs(IntervalGraph &graph, 
																 EdgePairMap &edgePairs,
																 ReadMateList &matePairs, 
																 ssize_t minMatePairCount, 
																 ssize_t ruleType) {
	//UNUSED+// ssize_t mp;
	ssize_t p ;
	ssize_t readIndex, mateIndex, mateType;
	//	int lastPathIntv, lastEdge, lastEdgeIntv, lastEdgePos;

	for (p = 0; p < graph.paths.size(); p++ ) {
		// Find the read that is paired with read 'p', and the
		// type of mate that pairs it.
		readIndex = p/2;
		mateIndex = matePairs[readIndex].mateIndex;
		mateType  = matePairs[readIndex].mateType;

		// Dont try and find a mate-path if no mates 
		// are mapped to this read
		if (mateIndex == -1)
			continue;
		
		if (ruleType != -1 and mateType != ruleType)
			continue;

		// Determine the path corresponding to the mate read

		ssize_t firstMate, secondMate;
		ssize_t mp;
		AssignMateOrder(p, mateIndex, mp, firstMate, secondMate);

			
		// If one of these paths is removed, don't try and 
		// find a mate-path.
		if (graph.pathLengths[firstMate] == 0 ||
				graph.pathLengths[secondMate] == 0) {
			continue;
		}
		ssize_t lastPathIntv = graph.pathLengths[firstMate] - 1;
		ssize_t lastEdge, lastEdgeIndex;
		lastEdge = graph.paths[firstMate][lastPathIntv].edge;
		lastEdgeIndex = graph.paths[firstMate][lastPathIntv].index;


		// use the reverse complment path of the mate since
		// the mate is sequenced in the opposite direction
		
		//UNUSED+// ssize_t mateStartEdgePos;
		ssize_t mateEdge, mateEdgeIntv ;

		//UNUSED// ssize_t mateReadPathLength = graph.pathLengths[secondMate];
	
		mateEdge     = graph.paths[secondMate][0].edge;
		mateEdgeIntv = graph.paths[secondMate][0].index;

		EdgePair edgePair;
		// Make a query-able matepair object.
		edgePair.edge1    = lastEdge;
		edgePair.edge2    = mateEdge;
		edgePair.mateType = mateType;
		EdgePairMap::iterator epIt;

		epIt = edgePairs.find(edgePair);
		assert(epIt != edgePairs.end());
		
		if ((*epIt).second.count < minMatePairCount) {
			matePairs[readIndex].mateIndex = -1;
			matePairs[mateIndex].mateIndex = -1;
		}
	}		
	EdgePairMap::iterator epIt, toErase;
	// Now remvoe the low coverage pairs.
	// Check the balance of the mate pairs.
	epIt = edgePairs.begin();
	while (epIt != edgePairs.end()) {
		if ((*epIt).second.count < minMatePairCount) {
			toErase = epIt;
			++epIt;
			edgePairs.erase(toErase);
		}
		else 
			++epIt;
	}
			

	// Check the balance of the mate pairs.
	for (epIt = edgePairs.begin(); epIt != edgePairs.end(); ++epIt) {
		ssize_t edge1Bal, edge2Bal;
		edge1Bal = graph.edges[(*epIt).first.edge1].balancedEdge;
		edge2Bal = graph.edges[(*epIt).first.edge2].balancedEdge;
		EdgePair ep;
		ep.edge1 = edge2Bal;
		ep.edge2 = edge1Bal;
		ep.mateType = (*epIt).first.mateType;
		if (edgePairs.find(ep) == edgePairs.end()) {
			std::cout << "pair: " << (*epIt).first.edge1 << " " << (*epIt).first.edge2
								<< " has no bal pair: "
								<<  edge2Bal << " " << edge1Bal << std::endl;
			//			exit(0);
		}
	}
}

void StoreMeanEdgePosition(EdgePairMap &edgePairs) {
	EdgePairMap::iterator epIt;
	for (epIt = edgePairs.begin(); epIt != edgePairs.end(); ++epIt) {
		ssize_t lastEdgePos = (*epIt).second.meanEdge1End / (*epIt).second.count;
		ssize_t mateStartEdgePos = (*epIt).second.meanEdge2Start / (*epIt).second.count;
		(*epIt).second.meanEdge1End = lastEdgePos;
		(*epIt).second.meanEdge2Start = mateStartEdgePos;

		// TODO: Why are these converting to floating arithmetic when they reduce to integer multiplication?
		// Is there supposed to be parentheses around (1.0 * (*epIt).second.count) ?
		// But how would that be any different from lastEdgePos, mateStartEdgePos computed above?
		// I'm fixing it to division, but it looks like the value is never used at all. -- GT
		//		(*epIt).second.read1Length = (int) ((*epIt).second.read1Length / 1.0* (*epIt).second.count);
		//		(*epIt).second.read2Length = (int) ((*epIt).second.read2Length / 1.0* (*epIt).second.count);

		(*epIt).second.read1Length = (*epIt).second.read1Length / (*epIt).second.count;
		(*epIt).second.read2Length = (*epIt).second.read2Length / (*epIt).second.count;
	}
}


ssize_t  PathsMayOverlap(Path &forPath, ssize_t forPathLength,
										 Path &revPath, ssize_t revPathLength,
										 ssize_t minCloneSize, ssize_t maxCloneSize ) {

	cout << "PathsMayOverlap isn't written yet" << endl;
	exit(0);
}

ssize_t  FindOverlappingMatePaths(Path &forPath, ssize_t forPathLength,
															Path &revPath, ssize_t revPathLength,
															Path &overlappingPath, ssize_t &overlappingPathLength) {
	cout << "FindOverlappingMatePaths isn't written yet." << endl;
	exit(1);


}

void StoreUniqueMatePaths(IntervalGraph &graph, TVertexList &vertices, TEdgeList &edges,
													RuleList &rules, 
													EdgePairMap &edgePairs, PathBranch &pathTree) {

	// Some output files for reporting statistics.

	// Hold and later print some statistics on mate-pair paths.
	std::vector<MatePathStatistics>  stats;	
	ssize_t numDeadEnd = 0;
	ssize_t numNoPath  = 0;
	EdgePairMap::iterator epIt;
	std::map<ssize_t,ssize_t> counts;
	for (epIt = edgePairs.begin(); epIt != edgePairs.end(); ++epIt) {

		ssize_t lastEdgePos = (*epIt).second.meanEdge1End;
		ssize_t mateStartEdgePos = (*epIt).second.meanEdge2Start;
		/*		std::out << "mate pair " << (*epIt).first.edge1 << " " << (*epIt).first.edge2 << " " 
							<< edges[(*epIt).first.edge1].length << " "
							<< edges[(*epIt).first.edge2].length << " "
							<< (*epIt).second.count << " " << lastEdgePos << " "
							<< mateStartEdgePos << std::endl;
		*/
		ssize_t count = (*epIt).second.count;
		ssize_t lastEdge, mateEdge;
		lastEdge = (*epIt).first.edge1;
		mateEdge = (*epIt).first.edge2;

		if (count > 2 and lastEdge != mateEdge ) {

			/*				std::cout << "checking for paths between joined mates: " << lastEdge << " (" << edges[lastEdge].length
									<< ") " << mateEdge << " (" << edges[mateEdge].length << ")  count: " << count << std::endl;
			*/
			// use the reverse complment path of the mate since
			// the mate is sequenced in the opposite direction

			// This could take much more time.
			ssize_t numValidPaths = 0;
			ssize_t foundValidPaths = 0;
			MatePathList matePath;
			ssize_t storeMatePath = 1;
			//UNUSED// ssize_t doPrint = 0;
			ssize_t mateType = (*epIt).first.mateType;
			ssize_t totalPathLength = 0;
			//UNUSED// int totalReadLength = (*epIt).second.read1Length + (*epIt).second.read2Length;
			//			int cloneSep = rules[mateType].cloneLength + totalReadLength;
			ssize_t cloneSep = rules[mateType].cloneLength;
			map<ssize_t,ssize_t> visited;
			foundValidPaths = CountValidMatePaths(graph,
																						lastEdge, lastEdgePos, // cur edge pos
																						mateEdge, mateStartEdgePos,
																						0,  // starting length is 0
																						cloneSep - rules[mateType].cloneVar, 
																						cloneSep + rules[mateType].cloneVar, 
																						10, // max depth to search
																						40, // max paths to find
																						numValidPaths, 
																						storeMatePath, matePath, 0, totalPathLength, visited);
			(*epIt).second.numPaths = numValidPaths;
			counts[numValidPaths]++;
			//			std::cout << "found " << numValidPaths << " paths." << std::endl;

			if (numValidPaths == 1 and (matePath.size() > 0 or lastEdge != mateEdge) ) {
				//				cout << "found one path for" << std::endl;
				MatePathList::iterator pathIt;
				PathInterval *path;
				ssize_t matePathLength = matePath.size() + 1;
				path = new PathInterval[matePathLength];
				//				std::cout << "mate path: " << lastEdge << " " << mateEdge << " : ";
				path[matePath.size()].edge = mateEdge;
				path[matePath.size()].index = -1;
				ssize_t pathPos = 0;
				for (pathIt = matePath.begin(); pathIt != matePath.end(); ++ pathIt) {
					//					cout << (*pathIt).edge << ", ";
					path[pathPos].edge = (*pathIt).edge;
					path[pathPos].index = -1;
					++pathPos;
				}
				//				std::cout << std::endl;

				CollectPathTreeOnPath(path, matePathLength, pathTree);
			
				// Add the balance of this path.
				PathInterval *rcPath;
				rcPath = new PathInterval[matePathLength];
				ssize_t i;
				for(i = 0; i < matePathLength; i++ ){ 
					rcPath[i].edge = edges[path[matePathLength - i - 1].edge].balancedEdge;
					rcPath[i].index = -1;
				}
				CollectPathTreeOnPath(rcPath, matePathLength, pathTree);
				delete[] path;
				delete[] rcPath;
			}
			else {
				if (numValidPaths <= 1) {
					if (vertices[edges[lastEdge].dest].OutDegree() == 0) {
						numDeadEnd++;
						//						disPathStats << "0 " << mateEdge << std::endl;
					}
					else {
						numNoPath++;
						//						disPathStats << (*epIt).second.count << " " << lastEdge << " " << mateEdge << " " << vertices[edges[lastEdge].dest].OutDegree() << std::endl;
					}
				}
			}
		}
	}
	//UNUSED// ssize_t s;
	/*
		for (s = 0; s < stats.size(); s++ ){
		std::cout << stats[s].numPaths << " " << stats[s].meanPathLength << std::endl;
		}
	*/
	//	std::cout << "num dead end: " << numDeadEnd << " no path: " << numNoPath << std::endl;
	std::map<ssize_t,ssize_t>::iterator countIt;
	/*
		std::cout << "mate path count is: " << std::endl;
		for (countIt = counts.begin(); countIt != counts.end(); ++countIt) {
		std::cout << (*countIt).first << " " << (*countIt).second << std::endl;
		}
	*/
}

ssize_t FindPairedScaffoldEdges(IntervalGraph &g, ReadMateList &readMates, ssize_t edge, 
														std::set<ssize_t> &pairedScaffoldEdges) {

	ssize_t pathIndex, readIndex, mateIndex;
	ssize_t i;

	for (i = 0; i < g.edges[edge].intervals->size(); i++ ){
		pathIndex = (*g.edges[edge].intervals)[i].read;
		readIndex = pathIndex / 2;
		if (readMates[readIndex].mateIndex != -1) {
			if (pathIndex % 2 == 0) {
				mateIndex = readMates[readIndex].mateIndex * 2 + 1;
				ssize_t pi;
				if (g.pathLengths[mateIndex] <= 0)
					continue;
				
				for (pi = 0; pi < g.pathLengths[mateIndex]; pi++) {
					ssize_t edgeIndex= g.paths[mateIndex][pi].edge;
					//UNUSED// ssize_t pathIndex= g.paths[mateIndex][pi].index;
					if (g.edges[edgeIndex].marked == GraphEdge::Marked) {
						// This edge is marked as a scaffolded edge.
						pairedScaffoldEdges.insert(edgeIndex);
					}
				}
			}
		}
	}
	return pairedScaffoldEdges.size();
}

ssize_t ComputeMatePairLengthDistribution(IntervalGraph g, ReadMateList &readMates,
																			 ssize_t mateType, ssize_t &meanSep, double &stddevSep) {

	ssize_t nMates = 0;
	ssize_t	maxSamples = 10000;
	ssize_t e;
	//UNUSED// ssize_t totalLength = 0;
	meanSep = 0;
	stddevSep   = 0.0;
	double sumSq = 0.0;
	for (e = 0; e < g.edges.size() and nMates < maxSamples; e++) {
		ssize_t pathIndex, readIndex, mateIndex;
		ssize_t i;
		for (i = 0; i < g.edges[e].intervals->size(); i++) {
			pathIndex = (*g.edges[e].intervals)[i].read;
			readIndex = pathIndex / 2;
			if (readIndex < 0 or readIndex >= readMates.size())
				continue;
			// Don't process all the cases that are not informative.
			if (g.pathLengths[pathIndex] == 0)
				continue;
			if (readMates[readIndex].mateIndex == -1)
				continue;
			if (readMates[readIndex].mateType != mateType)
				continue;
			if (pathIndex % 2 == 1)
				continue;

			mateIndex = readMates[readIndex].mateIndex * 2 + 1;
			if (g.pathLengths[mateIndex] == 0)
				continue;
			
			ssize_t mateIntv = g.paths[mateIndex][0].index;
			ssize_t mateEdge  = g.paths[mateIndex][0].edge;

			if (mateEdge != e) 
				continue;

			ssize_t matePos, readPos, readLength;

			readPos = (*g.edges[e].intervals)[i].edgePos;
			readLength = (*g.edges[e].intervals)[i].length;

			matePos = (*g.edges[e].intervals)[mateIntv].edgePos;
			if (readPos + readLength > matePos)
				continue;
				
			ssize_t mateSep;
			mateSep = (matePos - (readPos + readLength));

			sumSq += (mateSep * mateSep);
			meanSep += mateSep;
			nMates++;
		}
	}
	if (nMates == 0) {
		return 0;
	}
	meanSep /= nMates;
	double varSep;
	varSep = (sumSq / nMates - (meanSep * meanSep));
	stddevSep = sqrt(varSep);
	return 1;
}

void CollectMateEdges(IntervalGraph &g, ReadMateList &readMates, ssize_t edge, 
											MateEdgeMap &mateEdges, ssize_t mateType, ssize_t dir) {

	ssize_t i;
	ssize_t pathIndex, readIndex, mateIndex;

	for (i = 0; i < g.edges[edge].intervals->size(); i++ ){
		pathIndex = (*g.edges[edge].intervals)[i].read;
		readIndex = pathIndex / 2;
		/*
		if (edge == 5) 
			cout << "edge: " << edge << " path: " << pathIndex << " " << (*g.edges[edge].intervals)[i].pathPos 
					 << " " << readMates[readIndex].mateType << " " << mateType << endl;
		*/
		if (mateType != -1 and readMates[readIndex].mateType != mateType)
			continue;
		if (readMates[readIndex].mateIndex != -1) {
			if (dir == 0 and pathIndex % 2 == 0) {
				//
				// Collect edges in the forward direction.
				//

				// Either the read is in the forward direction
				// it is the first in a clone,
				//  -- or --
				// the read is a reverse complement read, and 
				// it is at the end of a mate-pair.
				// If this is not the 
				
				mateIndex = readMates[readIndex].mateIndex * 2 + 1;
			}
			else if (dir == 1 and pathIndex % 2 == 1) {
				mateIndex = readMates[readIndex].mateIndex * 2;
			}
			else {
				continue;
			}
			
			ssize_t pi;
			if (g.pathLengths[mateIndex] <= 0)
				continue;

			// 
			// Collect the edges paired to 'edge' simply by extending 
			// the read passing through 'edge'.
			//
			
			for (pi = pathIndex+ 1; pi < g.pathLengths[readIndex]; pi++) {
				ssize_t edgeIndex= g.paths[mateIndex][pi].edge;
				//UNUSED// ssize_t matePathIndex= g.paths[mateIndex][pi].index;
				if (mateEdges.find(edgeIndex) == mateEdges.end()) {
					mateEdges[edgeIndex].count = 1;
				}
				else {
					mateEdges[edgeIndex].count++;
				}
			}

			// 
			// Collect all the edges paired to 'edge' by mate-pairs.
			// This includes all edges along the paths that the mates
			// map to.
			//
			for (pi = 0; pi < g.pathLengths[mateIndex]; pi++ ){
				ssize_t edgeIndex= g.paths[mateIndex][pi].edge;
				ssize_t matePathIndex= g.paths[mateIndex][pi].index;
				if (edgeIndex != -1 and matePathIndex != -1) {
					if (mateEdges.find(edgeIndex) == mateEdges.end()) {
						mateEdges[edgeIndex].count = 1;
						mateEdges[edgeIndex].avgStartPos = (*g.edges[edge].intervals)[i].edgePos;
						mateEdges[edgeIndex].avgEndPos = (*g.edges[edgeIndex].intervals)[matePathIndex].edgePos +
							(*g.edges[edgeIndex].intervals)[matePathIndex].length;
					}
					else {
						mateEdges[edgeIndex].count++;
						mateEdges[edgeIndex].avgStartPos += (*g.edges[edge].intervals)[i].edgePos;
						mateEdges[edgeIndex].avgEndPos += (*g.edges[edgeIndex].intervals)[matePathIndex].edgePos +
							(*g.edges[edgeIndex].intervals)[matePathIndex].length;
					}/*
					if (edge == 5 and edgeIndex == 1) {
						cout << "5,1 pair " << dir << " " << pathIndex << " " << edge << " " << i << " " 
								 << (*g.edges[edge].intervals)[i].edgePos << " " 
								 << (*g.edges[edgeIndex].intervals)[matePathIndex].edgePos << " + " 
								 << (*g.edges[edgeIndex].intervals)[matePathIndex].length << endl;
								 }*/
				}
			} // end adding this path to the map.
				
			ssize_t beginPathLength;
			beginPathLength = g.pathLengths[pathIndex];

			// Count the number and positions that this edge maps to others.
			if (g.paths[mateIndex][0].edge != -1)
				mateEdges[g.paths[mateIndex][0].edge].cloneCount++;
		}
	}

	MateEdgeMap::iterator mateIt, mateEnd;
	mateEnd = mateEdges.end();
	for(mateIt = mateEdges.begin(); mateIt!=  mateEnd; ++mateIt) {
		if ((*mateIt).second.count > 0) {
			(*mateIt).second.avgEndPos /= (*mateIt).second.count;
			//		if ((*mateIt).second.cloneCount > 0) {
			(*mateIt).second.avgStartPos /= (*mateIt).second.count; //(*mateIt).second.cloneCount;
			
			//		}
			/*		else {
			// This no mate-pairs have started on this edge, so
			// we can't record an average start position.
			(*mateIt).second.avgStartPos = -1;
			}
			*/
		}
	}
}

void GetLastStartPosition(MateEdgeMap &mateEdgeMap, ssize_t &startPos, ssize_t &pairedEdge) {
	MateEdgeMap::iterator mapIt;
	startPos = -1;
	pairedEdge = -1;
	for (mapIt = mateEdgeMap.begin(); mapIt != mateEdgeMap.end(); ++mapIt) {
		if ((*mapIt).second.avgStartPos > startPos) {
			startPos = (*mapIt).second.avgStartPos;
			pairedEdge  = (*mapIt).first;
		}
	}
}


ssize_t GetAverageMateStartPos(IntervalGraph &g, ReadMateList &readMates, 
													 ssize_t srcEdge, ssize_t destEdge, ssize_t &avgSrcEnd, ssize_t &avgDestBegin) {

	ssize_t pathIndex, readIndex, mateIndex;
	ssize_t numSupported = 0;
	ssize_t i;
	avgSrcEnd = avgDestBegin = 0;
	for (i = 0; i < g.edges[srcEdge].intervals->size(); i++ ){
		pathIndex = (*g.edges[srcEdge].intervals)[i].read;
		readIndex = pathIndex / 2;
		if (readMates[readIndex].mateIndex != -1) {
			if (pathIndex % 2 == 0) {
				mateIndex = readMates[readIndex].mateIndex * 2 + 1;
				ssize_t mateContainsDest = 0;
				ssize_t pi;
				for (pi = 0; pi < g.pathLengths[mateIndex];  pi++) {
					if (g.paths[mateIndex][pi].edge == destEdge) {
						mateContainsDest = 1;
						avgDestBegin += (*g.edges[destEdge].intervals)[g.paths[mateIndex][pi].index].edgePos;
						avgSrcEnd += ((*g.edges[srcEdge].intervals)[i].edgePos + 
													(*g.edges[srcEdge].intervals)[i].length);
						numSupported++;
					}
				}
			}
		}
	}
	if (numSupported > 0) {
		avgDestBegin /= numSupported;
		avgSrcEnd    /= numSupported;
	}
	return numSupported;
}

ssize_t SearchForMateEdge(IntervalGraph &g, ssize_t rootEdge, ssize_t srcEdge, 
											ssize_t maxSearchLength, ssize_t maxSearchDepth, 
											ssize_t mateEdge, std::list<ssize_t> &path) {
	MateEdgeMap mateEdgeMap;
	mateEdgeMap[mateEdge].count = 1;
	ssize_t altDestEdge;
	return SearchForMateEdge(g, rootEdge, srcEdge, 
													 maxSearchLength, maxSearchDepth,
													 mateEdgeMap, path, altDestEdge);
}

// return values:
//    0: no paths found
//    1: unique path found
//    2: multiple paths found
ssize_t SearchForMateEdge(IntervalGraph &g, ssize_t rootEdge, ssize_t srcEdge, 
											ssize_t maxSearchLength, ssize_t maxSearchDepth, 
											MateEdgeMap &mateEdges, std::list<ssize_t> &path,
											ssize_t &altDestEdge) {
	// If this is on an invalid path.
	if (maxSearchLength < 0 or maxSearchDepth < 0) {
		return 0;
	}

	// Found a path ending in a mate edge.
	// store which edge this was.
	if (srcEdge != rootEdge) {
		MateEdgeMap::iterator mateEdgeIt;
		mateEdgeIt = mateEdges.find(srcEdge);
		if (mateEdgeIt != mateEdges.end()) {
			path.push_front(srcEdge);
			altDestEdge = (*mateEdgeIt).first;
			return 1;
		}
	}

	if (maxSearchLength >= 0 and maxSearchDepth >= 0) {
		
		ssize_t outEdge, outEdgeIndex;
		ssize_t srcEdgeDest = g.edges[srcEdge].dest;

		if (srcEdge != rootEdge) {
			// Only decrement the length if we are not checking edges from the source
			// edge.
			maxSearchLength -= (g.edges[srcEdge].length - g.vertices[srcEdgeDest].vertexSize);
			maxSearchDepth--;
		}
		// Break if this edge pushes us over the limit.
		if (maxSearchLength < 0 or maxSearchDepth < 0) 
			return 0;

		// Otherwise, move forward in the quest to find a valid edge.
		ssize_t numRoutesToEnd = 0;
		ssize_t retVal = 0;
		ssize_t pathEdge = -1;
		for (outEdgeIndex = g.vertices[srcEdgeDest].FirstOut();
				 outEdgeIndex != g.vertices[srcEdgeDest].EndOut();
				 outEdgeIndex = g.vertices[srcEdgeDest].NextOut(outEdgeIndex)) {
			outEdge = g.vertices[srcEdgeDest].out[outEdgeIndex];
			retVal = SearchForMateEdge(g, rootEdge, outEdge, 
																 maxSearchLength, maxSearchDepth, 
																 mateEdges, path, altDestEdge);
			if (retVal == 1) {
				numRoutesToEnd++;
				pathEdge = srcEdge;
			}
			else if (retVal == 2) {
				return 2;
			}
		}
		
		if (numRoutesToEnd == 1) {
			//
			// Case 1: only one route was found from here to the dest edge.
			//         Store that, and return 1 to signal that just a path was found.
			if (srcEdge != rootEdge) 
				path.push_front(srcEdge);
			return 1;
		}
		else if (numRoutesToEnd > 1) {
			// Case 2: Multiple routes were found from here to the dest edge.
			//         Don't bother storing the path edge, since it will be considered
			//         an invalid path.
			return 2;
		}
		else {
			// 
			// Case 3. No valid paths were found from this vertex.  No paths 
			//         to store, just return 0 to signal this.
			return 0;
		}
		return 0;
	}
	assert(0);
	return 0; // To quiet compiler warningse
}


void UntraverseMateEdges(IntervalGraph &g, MateEdgeMap &readMates,
												 std::set<ssize_t> &extraEdges) {
	// Step 1. Mark edgs in the read mates and extra edges as not traversed.
	MateEdgeMap::iterator mateEdgeIt, mateEdgeEnd;
	mateEdgeEnd = readMates.end();
	for (mateEdgeIt = readMates.begin(); mateEdgeIt != mateEdgeEnd; ++mateEdgeIt) {
		g.edges[(*mateEdgeIt).first].traversed = GraphEdge::NotMarked;
	}
	std::set<ssize_t>::iterator extraEdgeIt, extraEdgeEnd;
	extraEdgeEnd = extraEdges.end();
	for (extraEdgeIt = extraEdges.begin(); extraEdgeIt != extraEdgeEnd; ++extraEdgeIt) {
		g.edges[*extraEdgeIt].traversed = GraphEdge::NotMarked;
	}
}

void UnmarkMateEdges(IntervalGraph &g, MateEdgeMap &readMates,
										 std::set<ssize_t> &extraEdges) {
	// Step 1. Mark edgs in the read mates and extra edges as not traversed.
	MateEdgeMap::iterator mateEdgeIt, mateEdgeEnd;
	mateEdgeEnd = readMates.end();
	for (mateEdgeIt = readMates.begin(); mateEdgeIt != mateEdgeEnd; ++mateEdgeIt) {
		g.edges[(*mateEdgeIt).first].marked = GraphEdge::NotMarked;
	}
	std::set<ssize_t>::iterator extraEdgeIt, extraEdgeEnd;
	extraEdgeEnd = extraEdges.end();
	for (extraEdgeIt = extraEdges.begin(); extraEdgeIt != extraEdgeEnd; ++extraEdgeIt) {
		g.edges[*extraEdgeIt].marked = GraphEdge::NotMarked;
	}
}

void ClearDistances(MateEdgeMap &readMates, std::set<ssize_t> &extraEdges,
										std::vector<ssize_t> distToSrc) {
	// Step 1. Mark edgs in the read mates and extra edges as not traversed.
	MateEdgeMap::iterator mateEdgeIt, mateEdgeEnd;
	mateEdgeEnd = readMates.end();
	for (mateEdgeIt = readMates.begin(); mateEdgeIt != mateEdgeEnd; ++mateEdgeIt) {
		distToSrc[(*mateEdgeIt).first] = 0;
	}
	std::set<ssize_t>::iterator extraEdgeIt, extraEdgeEnd;
	extraEdgeEnd = extraEdges.end();
	for (extraEdgeIt = extraEdges.begin(); extraEdgeIt != extraEdgeEnd; ++extraEdgeIt) {
		distToSrc[(*extraEdgeIt)] = 0;
	}
}

ssize_t StoreTreeDepth(IntervalGraph &g, ssize_t srcEdge, 
									 MateEdgeMap &mateEdges, std::vector<ssize_t>  &distToEnd) {
	
	ssize_t outEdge, outEdgeIndex;
	ssize_t dest;
	g.edges[srcEdge].traversed = GraphEdge::Marked;
	dest = g.edges[srcEdge].dest;
	ssize_t longestDist = -1;
	ssize_t pathDist = -1;
	for (outEdgeIndex = g.vertices[dest].FirstOut();
			 outEdgeIndex != g.vertices[dest].EndOut();
			 outEdgeIndex = g.vertices[dest].NextOut(outEdgeIndex)) {
		outEdge = g.vertices[dest].out[outEdgeIndex];
		// Don't compute cycles
		if (g.edges[outEdge].traversed == GraphEdge::Marked)
			continue;
		// Next, if the out edge is marked as a puttative path-edge
		// in this mate-set, look for the end of the path.
		
		if (g.edges[outEdge].marked == GraphEdge::Marked) {
			//			std::cout << "storing tree depth." << std::endl;
			pathDist = StoreTreeDepth(g,  outEdge, mateEdges, distToEnd);
		}
		if (pathDist > longestDist) {
			longestDist = pathDist;
		}
	}
	if (longestDist == -1) {
		// This path is the end of the line.
		
		if (mateEdges.find(srcEdge) != mateEdges.end()) {
			longestDist = mateEdges[srcEdge].avgEndPos;
			distToEnd[srcEdge] = longestDist;
		}
		/*		else {
		// This edge is the root edge of the tree.
		distToEnd[srcEdge] = 0;
		}
		*/
	}
	else {
		distToEnd[srcEdge] = (g.edges[srcEdge].length 
													- g.vertices[dest].vertexSize 
													+ longestDist );
	}
	return distToEnd[srcEdge];
}


void FindTreeDepth(IntervalGraph &g, ssize_t srcEdge, 
									 MateEdgeMap &mateEdges, std::set<ssize_t> &extraEdges,
									 std::vector<ssize_t> &distToEnd) {
	UntraverseMateEdges(g, mateEdges, extraEdges);
	StoreTreeDepth(g, srcEdge, mateEdges, distToEnd);
	UnmarkMateEdges(g, mateEdges, extraEdges);
	UntraverseMateEdges(g, mateEdges, extraEdges);
}


ssize_t ExtractMin(std::map<ssize_t, ssize_t> &distMap, std::set<ssize_t> &removed) {
	// Ok, this is slow, but I don't feel like implementing 
	// a fib. heap for now.  This is just exploring a *small* subset
	// of the graph (max subgraph size is limited), so this should run quickly. 
	//

	std::map<ssize_t, ssize_t>::iterator distMapIt, distMapEnd;
	distMapEnd = distMap.end();
	if (distMap.size() == 0) {
		return -1;
	}

	distMapIt = distMap.begin();
	ssize_t minDist = SSIZE_MAX;
	ssize_t minDistVertex = -1;

	for (distMapIt = distMap.begin(); distMapIt != distMap.end(); ++distMapIt) {
		if (removed.find((*distMapIt).first) == removed.end() and 
				minDist > (*distMapIt).second) {
			minDist = (*distMapIt).second;
			minDistVertex = (*distMapIt).first;
		}
	}
	return minDistVertex;
}

ssize_t FindClosestVertex(IntervalGraph &g, ssize_t startVertex, std::set<ssize_t> &destVertices) {
	
	std::map<ssize_t,ssize_t> distMap;
	std::set<ssize_t> removed;
	distMap[startVertex] = 0;
	
	//UNUSED// ssize_t found;
	// just hard wire a stopping point now.
	ssize_t maxSearchSize = 100;
	
	ssize_t outEdge, outEdgeIndex;
	ssize_t bestDestDistance = SSIZE_MAX;
	ssize_t bestDestVertex;
	ssize_t curVertex = startVertex;
	
	// Loop while not all nodes have been processed,
	// and while not too many have been processed.
	//

	while (distMap.size() != removed.size() and
				 distMap.size() < maxSearchSize) {

		curVertex = ExtractMin(distMap, removed);

		bestDestVertex = -1;
		
		for (outEdgeIndex = g.vertices[curVertex].FirstOut();
				 outEdgeIndex != g.vertices[curVertex].EndOut();
				 outEdgeIndex = g.vertices[curVertex].NextOut(outEdgeIndex)) {
			outEdge = g.vertices[curVertex].out[outEdgeIndex];

			//
			// Look to see if the dest is part of the destVertices set, if so,
			// we're done.
			if (destVertices.find(g.edges[outEdge].dest) != destVertices.end()) {
				if (bestDestVertex == -1) {
					bestDestVertex = g.edges[outEdge].dest;
					bestDestDistance = distMap[curVertex] + g.edges[outEdge].length;
				}
				else {
					if (distMap[curVertex] + g.edges[outEdge].length < bestDestDistance) {
						bestDestVertex = g.edges[outEdge].dest;
						bestDestDistance = distMap[curVertex] + g.edges[outEdge].length;
					}
				}
			}

			// Otherwise, do the distance relaxations.
			if (distMap.find(g.edges[outEdge].dest) != distMap.end()) {
				if (distMap[g.edges[outEdge].dest] > distMap[curVertex] + g.edges[outEdge].length) {
					distMap[g.edges[outEdge].dest] = distMap[curVertex] + g.edges[outEdge].length;
				}
			}
			else {
				distMap[g.edges[outEdge].dest] = distMap[curVertex] + g.edges[outEdge].length;
			}
		}
		if (bestDestVertex != -1) {
			return bestDestVertex;
		}
		//
		// Done relaxing dest vertices reachable from 'curVertex', remove
		// curVertex from the queue.
		//
		//		distMap.erase(curVertex);
		removed.insert(curVertex);
		/*		std::cout << "cur dist map: ";
					std::transform(distMap.begin(), distMap.end(), 
					std::ostream_iterator<ssize_t>(std::cout, " "),
					__gnu_cxx::select1st<std::map<ssize_t, ssize_t>::value_type>() );
					std::cout << std::endl;
		*/
		
	}
	return -1;
}

ssize_t CollectPairedOutEdges(IntervalGraph &g, MateEdgeMap &pairedEdges, ssize_t curVertex, 
												 set<ssize_t> &pairedOutEdges) {

	ssize_t numPairedDestEdges = 0;
	ssize_t destOutEdgeIndex, destOutEdge;
	for (destOutEdgeIndex = g.vertices[curVertex].FirstOut();
			 destOutEdgeIndex != g.vertices[curVertex].EndOut();
			 destOutEdgeIndex = g.vertices[curVertex].NextOut(destOutEdgeIndex)) {
		destOutEdge = g.vertices[curVertex].out[destOutEdgeIndex];
		if (pairedEdges.find(destOutEdge) != pairedEdges.end()) {
			numPairedDestEdges++;
			pairedOutEdges.insert(destOutEdge);
		}
	}
	return numPairedDestEdges;
}

ssize_t FindPairedOutEdge(IntervalGraph &g, MateEdgeMap &pairedEdges, ssize_t curVertex, ssize_t &pairedOutEdge) {
	// 
	// Given a vertex, lookup the out edge that is stored  marked as
	// a paired edge in pairedEdges.  Return the number of out edges
	// that are paired, since some methods need to follow unique paths.
	ssize_t destOutEdge, destOutEdgeIndex;
	ssize_t numPairedDestEdges = 0;
	pairedOutEdge = -1;
	for (destOutEdgeIndex = g.vertices[curVertex].FirstOut();
			 destOutEdgeIndex != g.vertices[curVertex].EndOut();
			 destOutEdgeIndex = g.vertices[curVertex].NextOut(destOutEdgeIndex)) {
		destOutEdge = g.vertices[curVertex].out[destOutEdgeIndex];
		if (pairedEdges.find(destOutEdge) != pairedEdges.end()) {
			numPairedDestEdges++;
			pairedOutEdge = destOutEdge;
		}
	}
	return numPairedDestEdges;
}


ssize_t FindMatePath(IntervalGraph &g, ssize_t srcEdge,
								 MateEdgeMap &srcMateEdges,
								 std::list<ssize_t> &srcEdgePath) {
	// Step 1. Mark edges that most likely reach the src edge.  It's possible
	//         that the correct path is not entirely covered by mate edges, furthermore,
	//         it's possible that there are erroneous edges in the srcMateEdges set.
	//         The goal of creating a path tree is to detect which edges that leave the
	//         src are invalid.
	srcEdgePath.push_back(srcEdge);
	std::set<ssize_t> pairedSrcVertices;
	MateEdgeMap::iterator mateEdgeIt, mateEdgeEnd;
	mateEdgeEnd = srcMateEdges.end();
	//	std::cout << "finding mate path from: " << srcEdge << " (" << g.edges[srcEdge].length << "): ";
	for (mateEdgeIt = srcMateEdges.begin(); mateEdgeIt != mateEdgeEnd; ++mateEdgeIt) {
		if ((*mateEdgeIt).first != srcEdge) {
			g.edges[(*mateEdgeIt).first].marked = GraphEdge::Marked;
			g.edges[(*mateEdgeIt).first].traversed = GraphEdge::NotMarked;
			pairedSrcVertices.insert(g.edges[(*mateEdgeIt).first].src);
			//			std::cout << (*mateEdgeIt).first << " (" << (*mateEdgeIt).second.count << ") ";
		}
	}
	//	std::cout << std::endl;
	
	ssize_t curEdge = srcEdge;
	ssize_t multiplePathsFound = 0;
	if (srcMateEdges.find(srcEdge) != srcMateEdges.end()) {
		srcMateEdges.erase(srcEdge);
	}

	while(srcMateEdges.size() > 0 and !multiplePathsFound and curEdge != -1) {
		//
		// Count how many out edges from 'curEdge' are paired to
		// the source edge.
		//
		
		ssize_t numPairedOutEdges = 0;
		//UNUSED// ssize_t outEdge, outEdgeIndex;
		ssize_t curEdgeDest = g.edges[curEdge].dest;
		ssize_t nextEdge = -1;		
		//
		// This loop attempts to update 'nextEdge' with 
		// the that follows 'curEdge' on the mate path.
		//

		MateEdgeMap::iterator pairedEdgeIt;

		numPairedOutEdges = FindPairedOutEdge(g, srcMateEdges, curEdgeDest, nextEdge);

		if (numPairedOutEdges > 1) {
			// There are multiple paths from this one, don't trying to create a
			// mate path from this since is is ambiguous.
			srcEdgePath.clear();
			
			// For now, mark the sentinal that multiple paths are found, but just
			// bail with a failed search.  It's possible that later on
			// more will be done to update the search 
			multiplePathsFound = 1;
			//			std::cout << " multiple (" << numPairedOutEdges << ") paired out edges." << std::endl;
			return 0;
		}

		else if (numPairedOutEdges == 1) {
			srcEdgePath.push_back(nextEdge);
		}
		else {
			// No paired out edges were found that are compatible with this one.
			// Try and find a path to one of the 
			
			ssize_t closestPairedVertex = FindClosestVertex(g, g.edges[curEdge].dest, pairedSrcVertices);

			if (closestPairedVertex == -1) {
				// There was no vertex that was found, so 
				// there is no hope of tracing a path out from this edge, so no mate
				// path may be found.
				//				std::cout << "no close paired vertex. " << std::endl;
				return 0;
			}
			
			//
			// Otherwise, it is possible that there is a path from curEdge to
			// the closest mate path.  Look for that path.
			//

			// First, we found the closest vertex, not the closest edge,
			// look at the closest vertex to find the closest edge.
			// If there are multiple edges from this vertex, don't try.

			ssize_t pairedDestEdge, numPairedDestEdges;
			numPairedDestEdges = FindPairedOutEdge(g, srcMateEdges, closestPairedVertex,  pairedDestEdge);
			if (numPairedDestEdges > 1 ) {
				// There are two paired edges from the same vertex, so 
				// the mate-path is ambiguous (there may just be a cycle, but
				// that's a bit complicated for now).
				//				std::cout << "multiple (" << numPairedDestEdges << ") paired dest edges." << std::endl;
				return 0;
			}
			
			// There should have been one dest edge that is paired
			// with the source. 
			if (numPairedDestEdges != 1) {
				std::set<ssize_t>::iterator setIt;
				/*				std::cout << "no good edge found at: " << std::endl;
				for (setIt = pairedSrcVertices.begin();
						 setIt != pairedSrcVertices.end();
						 ++setIt) {
					std::cout << (*setIt) << " ";
				}
				std::cout << std::endl;
				std::cout << "searched list: ";
				MateEdgeMap::iterator mapit;
				for (mapit = srcMateEdges.begin(); mapit != srcMateEdges.end(); ++mapit) {
					std::cout << (*mapit).first << " ";
				}
				std::cout << std::endl;
				*/
			}
					
			assert(numPairedDestEdges == 1);

			ssize_t searchRetVal;
			std::list<ssize_t> pathToPairedDest;
			searchRetVal = SearchForMateEdge(g, curEdge, curEdge, 200, 5, 
																			 pairedDestEdge, pathToPairedDest);
			
			if (searchRetVal == 0) {
				// This shouldn't happen.  We were able to find a paht to the 
				// vertex, but not the edge.
				//				assert(0);
			}
			else if (searchRetVal == 2) {
				// There was a gap in the mate-path, but
				// there are multiple ways to fill it, so the path
				// is ambiguous.
				//				std::cout << "multiple patch paths." << std::endl;
				srcEdgePath.clear();
				return 0;
			}
			else {
				assert(searchRetVal == 1);
				std::list<ssize_t>::iterator pathIt;
				std::cout << "found extra path from " << curEdge<< " to " << pairedDestEdge << std::endl;
				for (pathIt = pathToPairedDest.begin(); pathIt != pathToPairedDest.end(); ++pathIt) {
					std::cout << *pathIt << " ";
					srcEdgePath.push_back(*pathIt);
					nextEdge = *pathIt;
				}
				std::cout << std::endl;
			}
		} // Done searching for an alternative path.
		if (nextEdge != -1) {
			// The next edge should be removed from the list of candidate next edges
			// since it has been added to the path and loops aren't cosidered for now.
			std::cout << "done with next edge: " << nextEdge << std::endl;
			pairedSrcVertices.erase(g.edges[nextEdge].src);
			srcMateEdges.erase(nextEdge);
			
		}
		curEdge = nextEdge;			
	}

	if (srcMateEdges.size() == 0 and !multiplePathsFound) {
		// All edges in srcMateEdges have been processed, so a single path 
		// has been found 
		std::cout << "found src edge path: " << std::endl;
		std::list<ssize_t>::iterator pathIt;
		for (pathIt = srcEdgePath.begin(); pathIt != srcEdgePath.end(); ++pathIt) {
			std::cout << *pathIt << " ";
		}
		std::cout << std::endl;
		return 1;
	}
	return 0;
}



void RemoveLowCountMateEdges(MateEdgeMap &mateEdges, ssize_t minCount) {
	MateEdgeMap::iterator mateEdgeIt, deletedMateEdgeIt;
	mateEdgeIt = mateEdges.begin();
	while (mateEdgeIt != mateEdges.end()) {
		if ((*mateEdgeIt).second.count < minCount) {
			deletedMateEdgeIt = mateEdgeIt;
			++mateEdgeIt;
			mateEdges.erase(deletedMateEdgeIt);
		}
		else {
			++mateEdgeIt;
		}
	}
}

void FindScaffoldPaths(IntervalGraph &g, 
											 ReadMateList &mateTable,
											 RuleList &ruleList, ssize_t mateType, ssize_t scaffoldEdgeLength,
											 PathIntervalList &paths,
											 PathLengthList &pathLengths) {
	ssize_t e;
	PathInterval *path, *pathRC;
	ssize_t numOptimalPaths = 0;
	ssize_t numTotalPairs   = 0;
	for (e = 0; e < g.edges.size(); e++ ) {
		g.edges[e].traversed = GraphEdge::NotMarked;
	}

	for (e = 0; e < g.edges.size(); e++ ) {
		
		/*
		cout << "Finding paths for: " << e << " " 
				 << (ssize_t) g.edges[e].marked  << " not " << (ssize_t) GraphEdge::NotMarked 
				 << " " << g.edges[e].length << " <? " << scaffoldEdgeLength << endl;
		*/
		// Only attempt to scaffold edges that are marked as being possible
		// for scaffold edges.
		if (g.edges[e].marked == GraphEdge::NotMarked)
			continue;

		// 
		// An edge is set to traversed if the search for 
		// mate paths is conducted from the balanced
		// edge first.
		if (g.edges[e].traversed == GraphEdge::Marked)
			continue;

		MateEdgeMap::iterator mit, mpeIt;
		MateEdgeMap mateEdges;
		MateEdgeMap matePairEdges;

		//UNUSED// ssize_t destEdge;
		if (g.edges[e].length < scaffoldEdgeLength)
			continue;
		// collect mate-pairs from the src edge
		//		cout << "collecting edges for: " << e << endl;
		CollectMateEdges(g, mateTable, e, mateEdges, mateType, 0); 
		
		RemoveLowCountMateEdges(mateEdges, 3);

		//
		// Each edge that a src is paired with is a potential
		// dest edge.  Try and find 
		for (mit = mateEdges.begin(); mit != mateEdges.end(); ++mit) {
			//			cout << "  got mate edge: " << (*mit).first << endl;
			if ((*mit).first == e)
				continue;
			if (g.edges[(*mit).first].marked != GraphEdge::Marked)
				continue;
			

			std::list<ssize_t> srcEdgePath;
			matePairEdges = mateEdges;
			// Add in mate-pairs from the dest edge
			CollectMateEdges(g, mateTable, (*mit).first, matePairEdges, mateType, 1); 
			RemoveLowCountMateEdges(matePairEdges, 3);
			//			cout << endl << endl << endl;


			
			// Only scaffold on long edges.
			if (g.edges[(*mit).first].length < scaffoldEdgeLength)
				continue;
			
			//
			// Try and find mate-paths from src (e) to dest ((*mit).first).  
			// 
			ssize_t srcEdge  = e;
			ssize_t destEdge = (*mit).first;
			/*
			cout << "looking for path from " << srcEdge << " [" << g.edges[srcEdge].index << ", "
					 << g.edges[srcEdge].length << "] "
					 << " to " << destEdge << " [" << g.edges[destEdge].index << ", "
					 << g.edges[destEdge].length << "] using rule: " << mateType << endl;
			*/
			MateEdgeMap::iterator mapIt;
			/*
			for (mapIt = mateEdges.begin(); mapIt != mateEdges.end(); ++mapIt) {
				cout << (*mapIt).first << " " << (*mapIt).second.count << " " 
						 << (*mapIt).second.avgStartPos << " " << (*mapIt).second.avgEndPos << endl;
			}
			*/
			if (FindMaximallySupportedPath(g, mateTable, ruleList, mateType, srcEdge, destEdge, srcEdgePath)) {
				// 
				// A path from src to dest was found and stored in 
				// the list srcEdgePath.  Convert that to a mate-path.
				//
				path   = new PathInterval[srcEdgePath.size()];
				pathRC = new PathInterval[srcEdgePath.size()];
				ssize_t pathLength = srcEdgePath.size();
				std::list<ssize_t>::iterator edgeIt;
				ssize_t pathPos = 0;

				//cout << "found supported path: " << endl;
				for (edgeIt = srcEdgePath.begin();
						 edgeIt != srcEdgePath.end();
						 ++edgeIt) {
					path[pathPos].edge = *edgeIt;
					//					cout << *edgeIt << " ";
					pathRC[pathLength - pathPos - 1].edge = g.edges[(*edgeIt)].balancedEdge;
					++pathPos;
				}
				//				cout << endl;
				// 
				// Compute statistics about mate-pairs converted to this path.
				//
				ssize_t i;
				for (i = 0; i < g.edges[srcEdge].intervals->size(); i++) { 
					ssize_t pathIndex;
					ssize_t mateIndex;
					pathIndex = (*g.edges[srcEdge].intervals)[i].read;
					// process reads in forward strand here.
					if (pathIndex % 2 != 0)	continue;
					mateIndex = mateTable[pathIndex/2].mateIndex;
					if (mateIndex == -1) continue;
					
					ssize_t matePath = mateIndex * 2 + 1;
					ssize_t mp, pi;
					for (pi = 0; pi < pathLength; pi++) {
						for (mp = 0; mp < g.pathLengths[matePath]; mp++) {
							if (path[pi].edge == g.paths[matePath][mp].edge) {
								mateTable[pathIndex/2].marked = 1;
								mateTable[mateIndex].marked = 1;
							}
						}
					}
				}
				// process reads that are incorporated into the 
				// dest edge, which are RC edges.
				for (i = 0; i < g.edges[destEdge].intervals->size(); i++) { 
					ssize_t pathIndex;
					ssize_t mateIndex;
					pathIndex = (*g.edges[destEdge].intervals)[i].read;
					// process reads in forward strand here.
					if (pathIndex % 2 != 1)	continue;
					mateIndex = mateTable[pathIndex/2].mateIndex;
					if (mateIndex == -1) continue;

					ssize_t matePath = mateIndex * 2 ;
					ssize_t mp, pi;
					for (pi = 0; pi < pathLength; pi++) {
						for (mp = 0; mp < g.pathLengths[matePath]; mp++) {
							if (path[pi].edge == g.paths[matePath][mp].edge) {
								mateTable[pathIndex/2].marked = 1;
								mateTable[mateIndex].marked = 1;
							}
						}
					}
				}				
						
				
				paths.push_back(path);
				paths.push_back(pathRC);
				pathLengths.push_back(pathLength);
				pathLengths.push_back(pathLength);
				numOptimalPaths++;
			}

			// Unmark the source and the balanced edge of the dest so
			// that ambiguous paths are not defined.

			g.edges[g.edges[(*mit).first].balancedEdge].traversed = GraphEdge::Marked;
			numTotalPairs++;
		}
	}
	cout << "out of " << numTotalPairs << ", " << numOptimalPaths << " had connecting paths." << std::endl;
}

ssize_t AdvancePathAlongPairedEdges(IntervalGraph &g, ssize_t curEdge, MateEdgeMap &pairedEdges, 
																ssize_t stopEdge, 
																std::list<ssize_t> &pairedPath, ssize_t &pathSeqLength,
																std::map<ssize_t,ssize_t> &edgeTraversals,
																ssize_t &lastEdge, ssize_t &numPairedOutEdges) {
	
	// Move forward along a path that is marked by paired edges.
	// This assumes 'curEdge' is already on the path, and so the first 
	// edge that is added to the path is an edge following curEdge.
	ssize_t destVertex = g.edges[curEdge].dest;
	ssize_t pairedOutEdge;
	ssize_t numNewEdges = 0;

	numPairedOutEdges = 0;

	if (pairedEdges.find(curEdge) == pairedEdges.end())
		return 0;

	
	//
	// This path must have started off on a paired edge.
	// 
	//UNUSED// ssize_t pathLength = 1;
	pairedPath.push_back(curEdge);

	// Only one edge is paired to the curEdge.
	numPairedOutEdges = 1; 

	// Record how long the path is extended.
	++numNewEdges;
	edgeTraversals[curEdge]++;
	// The curEdge is paired.
	pairedOutEdge = curEdge;
	lastEdge      = pairedOutEdge;
	if (edgeTraversals[curEdge] > MAX_EDGE_TRAVERSALS) {
		return numNewEdges;
	}
	// Move forward as long as there is no ambiguous path.
	while (numPairedOutEdges == 1 and pairedOutEdge != stopEdge) {	 
		destVertex = g.edges[pairedOutEdge].dest;

		numPairedOutEdges = FindPairedOutEdge(g, pairedEdges, destVertex, pairedOutEdge);
		if (numPairedOutEdges != 1)
			break;

		if (pairedOutEdge != stopEdge)
			// This is an internal edge.  Increment the path length. 
			// The final edge length will be added when checking to 
			// see if this is a valid path because the 
			// position of the mate pairs is stored on the last edge.
			pathSeqLength += g.edges[pairedOutEdge].length - g.vertices[destVertex].vertexSize;			

		pairedPath.push_back(pairedOutEdge);
		lastEdge = pairedOutEdge;
		edgeTraversals[curEdge]++;
		if (edgeTraversals[curEdge] > MAX_EDGE_TRAVERSALS) {
			return numNewEdges;
		}
		++numNewEdges;
	}


	return numNewEdges;
}


void MarkScaffoldEdges(IntervalGraph &g, ssize_t minCoverage, ssize_t minLength) {
	ssize_t e;
	for (e = 0; e < g.edges.size(); e++ ){ 
		// Reset this just in case
		g.edges[e].marked = GraphEdge::NotMarked;
		if ((*g.edges[e].intervals).size() >= minCoverage and g.edges[e].length > minLength) {
			g.edges[e].marked = GraphEdge::Marked;
		}
	}
}


ssize_t IsLengthValid(ssize_t length, ssize_t min, ssize_t max) {
	return length >= min && length <= max;
}

ssize_t StoreMaximalPath(std::list<ssize_t> &curSupportedPath, ssize_t curSupportedPathScore,
										 std::list<ssize_t> &maxSupportedPath, ssize_t &maxSupportedPathScore) {
	// Replace maxspanningpath by curspanning path if possible.
	if (curSupportedPathScore > maxSupportedPathScore) {
		maxSupportedPath.clear();
		maxSupportedPath = curSupportedPath;
		maxSupportedPathScore = curSupportedPathScore;
		return maxSupportedPathScore;
	}
	else {
		return 0;
	}
}
 
ssize_t FindMaximallySupportedPath(IntervalGraph &g, MateEdgeMap &pairedEdges, 
															 ssize_t curEdge, ssize_t curEdgePos, ssize_t destEdge, ssize_t destEdgePos,
															 ssize_t curLength, ssize_t minPathLength, ssize_t maxPathLength,
															 ssize_t maxSearchDepth,
															 std::list<ssize_t> &curSupportedPath, ssize_t curSupportedPathScore,
															 std::list<ssize_t> &maxSupportedPath, ssize_t &maxSupportedPathScore,
															 std::map<ssize_t,ssize_t> &edgeTraversals,
															 ssize_t &numOptPaths) {
	// 
	// Given a graph, a source edge, and a dest edge, find a path from
	// source to dest that has the most support.
	// Support may be defined using mate-pairs from the source edge or the dest edge,
	// and may be counted as the number of edges along the path that provide support
	// for the path, or the number of reads that support the path.  I'm not
	// sure which will be the best.  
	//

	// Limit this from searching too much of the graph.
	if (maxSearchDepth == 0)
		return 0;

	if (curLength > maxPathLength) 
		return 0;
	/*
	cout << "FMSP: cur " << curEdge << " pos " << curEdgePos 
			 << " dest " << destEdge << " destpos: " << destEdgePos 
			 << " curlen " << curLength << " min: " << minPathLength <<  " max: " << maxPathLength
			 << " maxdep " << maxSearchDepth << endl;
	*/
	if (curEdge == destEdge and destEdgePos > curEdgePos) {
		// Reached the dest.  If the current path has more supported edges than the
		// current optimal path, replace opt with cur path.
		ssize_t curPathLength = destEdgePos - curEdgePos + curLength;
		if (IsLengthValid(curPathLength, minPathLength, maxPathLength)) {
			if (StoreMaximalPath(curSupportedPath, curSupportedPathScore,
													 maxSupportedPath, maxSupportedPathScore)) {
				numOptPaths = 1;
			}
			else if (curSupportedPathScore == maxSupportedPathScore) {
				numOptPaths++;
			}
			//			cout << " num opt paths: " << numOptPaths << endl;
			return maxSupportedPathScore;
		}
		else {
		}
		return 0;
	}
	
	
	ssize_t dest = g.edges[curEdge].dest;
	//UNUSED// ssize_t nextEdge;
	ssize_t pathSeqLength;
	ssize_t outEdge, outEdgeIndex;

	set<ssize_t> pairedOutEdges;
	CollectPairedOutEdges(g, pairedEdges, dest, pairedOutEdges);
	if (pairedOutEdges.size() == 0) {
		// Try all paired out edges.
		//		cout << "FMSP no paired out edges, trying all." << endl;
		for (outEdgeIndex = g.vertices[dest].FirstOut();
				 outEdgeIndex != g.vertices[dest].EndOut();
				 outEdgeIndex = g.vertices[dest].NextOut(outEdgeIndex)) {
			pairedOutEdges.insert(g.vertices[dest].out[outEdgeIndex]);
		}
	}

	set<ssize_t>::iterator poeIt;
	for (poeIt = pairedOutEdges.begin(); poeIt != pairedOutEdges.end(); ++poeIt) {
		outEdge = *poeIt;
		//		cout << "FMSP edge traversals: " <<edgeTraversals[outEdge] << endl;
		// If there are too many cycles, don't try and search.
		if (edgeTraversals[outEdge] > MAX_EDGE_TRAVERSALS)
			return 0;

		edgeTraversals[outEdge]++;

		pathSeqLength = 0;
		//UNUSED// ssize_t pairedPathLength;
		//UNUSED// ssize_t lastEdge;

		//
		// Update the statistics from the concatenated path.
		//
		ssize_t outEdgeScore = 0;
		if (pairedEdges.find(outEdge) != pairedEdges.end()) {
			// As a heuristic, give the search a new chance.
			//			cout << "unsupported edge. " << endl;
			maxSearchDepth = 6; 
			outEdgeScore = 1;
		}

		//		cout << "advancing along edge: " << outEdge << " " << pathSeqLength << endl;
		curSupportedPath.push_back(outEdge);
		
		if (pathSeqLength > maxPathLength)
			return 0;


		FindMaximallySupportedPath(g, pairedEdges, outEdge, 0, destEdge, destEdgePos,
															 (curLength 
																+ g.edges[curEdge].length 
																- curEdgePos
																- g.vertices[g.edges[curEdge].dest].vertexSize),
															 minPathLength, maxPathLength,
															 maxSearchDepth - 1,
															 curSupportedPath, curSupportedPathScore + outEdgeScore, 
															 maxSupportedPath, maxSupportedPathScore, 
															 edgeTraversals,
															 numOptPaths);
		// Get rid of the cur edge from the path stack.
		curSupportedPath.pop_back();
	}

	return maxSupportedPathScore; // TODO: check if we can use void return value
}

ssize_t FindMaximallySupportedPath(IntervalGraph &g, ReadMateList &mateTable, 
															 RuleList &mateRules, ssize_t mateType,
															 ssize_t srcEdge, ssize_t destEdge,
															 std::list<ssize_t> &optimalPath) {
	MateEdgeMap srcMateEdges;
	std::set<ssize_t> srcMateEdgeSet;
	/*
	cout << "Starting FMSP " << srcEdge << " " << destEdge << endl;
	*/
	CollectMateEdges(g, mateTable, srcEdge, srcMateEdges, mateType);
	RemoveLowCountMateEdges(srcMateEdges, 3);
	
	std::list<ssize_t> curPath;

	assert(srcMateEdges.find(destEdge) != srcMateEdges.end());
	MateEdgeMap::iterator mapIt;

	ssize_t srcEdgePos  = srcMateEdges[destEdge].avgStartPos;
	ssize_t destEdgePos = srcMateEdges[destEdge].avgEndPos;

	ssize_t numPathPairedEdges, optPathPairedEdges;
	numPathPairedEdges = optPathPairedEdges = 0;
	ssize_t numOptimalPaths = 0;
	ssize_t optimalPathScore;
	
	// The optimal path will always have src edge.
	curPath.push_back(srcEdge);
	map<ssize_t,ssize_t> edgeTraveralCount;
	optimalPathScore = FindMaximallySupportedPath(g, srcMateEdges,
																								srcEdge, srcEdgePos, destEdge, destEdgePos,
																								0, // cur length
																								mateRules[mateType].cloneLength - mateRules[mateType].cloneVar, 
																								mateRules[mateType].cloneLength + mateRules[mateType].cloneVar, 
																								10, // max search depth.  Don't allow too many gaps in the path.
																								curPath, numPathPairedEdges,
																								optimalPath, optPathPairedEdges, 
																								edgeTraveralCount, numOptimalPaths);
	ssize_t retval;
	if (numOptimalPaths == 1 ) {
		retval = 1;
	}
	else {
		retval = 0;
	}

	return retval;
}


ssize_t QueryUpstreamEdgesForDownstreamConnection(IntervalGraph &g, ReadMateList &mateTable,
																							ssize_t curVertex, ssize_t radius,
																							std::set<ssize_t> &dsEdgeSet, ssize_t minPairedEdges,
																							std::set<ssize_t> &traversedVertices) {
	if (traversedVertices.find(curVertex) !=
			traversedVertices.end()) {
		return 0;
	}
	traversedVertices.insert(curVertex);
		
	// This search has gone too far.
	if (radius <= 0)
		return 0;
	// Search edges backwards form cur vertex for edges that contain mate-pairs
	// somewhere in the dsEdgeSet.
	//UNUSED+// ssize_t  inEdge;
	ssize_t inEdgeIndex;
	for (inEdgeIndex = g.vertices[curVertex].FirstIn();
			 inEdgeIndex != g.vertices[curVertex].EndIn();
			 inEdgeIndex = g.vertices[curVertex].NextIn(inEdgeIndex)) {

		ssize_t intv;
		//UNUSED// ssize_t path;
		ssize_t pathIndex, readIndex, mateIndex;
		ssize_t matePathPos;
		ssize_t numPairedEdges = 0;
		ssize_t edgeLength = g.edges[inEdgeIndex].length;
		for (intv = 0; intv < g.edges[inEdgeIndex].intervals->size(); intv++) {
			pathIndex = (*g.edges[inEdgeIndex].intervals)[intv].read;
			readIndex = pathIndex / 2;
			mateIndex = mateTable[readIndex].mateIndex;
			ssize_t edgePos = (*g.edges[inEdgeIndex].intervals)[intv].edgePos;
			
			if (mateIndex == -1 or (edgeLength - edgePos < radius))
				continue;
			
			if (readIndex % 2 == 0) {
				ssize_t readMateIndex = mateIndex * 2 + 1;
				for (matePathPos = 0; matePathPos < g.pathLengths[readMateIndex]; ++matePathPos) {
					ssize_t mateEdge;
					mateEdge = g.paths[readMateIndex][matePathPos].edge;
					if (dsEdgeSet.find(mateEdge) != dsEdgeSet.end()) {
						++numPairedEdges;
						if (numPairedEdges >= minPairedEdges) {
							return 1;
						}
					} // end searching for paired edges
				} // end searching through mate paths for this read.
			} // end checking this read
		} // end checking all reads on this edge.
		ssize_t srcVertex = g.edges[inEdgeIndex].src;
		if (QueryUpstreamEdgesForDownstreamConnection(g, mateTable, curVertex, 
																									radius - (g.edges[inEdgeIndex].length -
																														 g.vertices[srcVertex].vertexSize),
																									dsEdgeSet, minPairedEdges, traversedVertices)) {
			return 1;
		}
	}
	return 0;
}


void StoreDownstreamEdges(IntervalGraph &g, ssize_t curVertex, ssize_t radius, 
													std::set<ssize_t> &dsEdgeSet) {
	if (radius < 0) 
		return;
	ssize_t outEdge, outEdgeIndex;
	for (outEdgeIndex = g.vertices[curVertex].FirstOut();
			 outEdgeIndex != g.vertices[curVertex].EndOut();
			 outEdgeIndex = g.vertices[curVertex].NextOut(outEdgeIndex)) {
		outEdge = g.vertices[curVertex].out[outEdgeIndex];
		// This edge may have already been traversed.
		if (dsEdgeSet.find(outEdge) != dsEdgeSet.end()) 
			continue;
		/*		std::cout << "storing downstream with radius: " << radius << " "
							<< outEdge << " "
							<< g.edges[outEdge].length << endl;
		*/
		dsEdgeSet.insert(outEdge);
		StoreDownstreamEdges(g, g.edges[outEdge].dest, 
												 radius - (g.edges[outEdge].length - 
																	 g.vertices[g.edges[outEdge].dest].vertexSize),
												 dsEdgeSet);
	}
}

ssize_t RemoveBadStartEndMatePairs(IntervalGraph &g, ReadMateList &mateTable, RuleList &rules, ssize_t ruleType) {
	// Now detect mate-pairs with invalid start/end positions on edges, and remove them.
	ssize_t i, e;
	ssize_t numRemovedMatePairs = 0;
	for (e = 0; e < g.edges.size(); e++ ){ 
		for (i = 0; i < (*g.edges[e].intervals).size(); i++ ) {
			// Determine 
			ssize_t path = (*g.edges[e].intervals)[i].read;
			//UNUSED// int pathPos = (*g.edges[e].intervals)[i].pathPos;
			
			if (path % 2 == 0) {
				// This is the forward read.
				ssize_t readIndex = path / 2;
				ssize_t mateIndex = mateTable[readIndex].mateIndex;
				ssize_t mateReadIndex = mateIndex * 2 + 1;

				if (mateIndex == -1) continue;
				if (g.pathLengths[path] == 0 or g.pathLengths[mateReadIndex] == 0) {
					continue;
				}
				if (ruleType != -1 and mateTable[readIndex].mateType != ruleType)
					continue;

				ssize_t pathLength  = g.pathLengths[path];
				ssize_t readEndEdge = g.paths[path][pathLength-1].edge;
				ssize_t readEndIntv = g.paths[path][pathLength-1].index;
				ssize_t readEndPos  = (*g.edges[readEndEdge].intervals)[readEndIntv].edgePos + 
					(*g.edges[readEndEdge].intervals)[readEndIntv].length;
				
				ssize_t readEndEdgeLength = g.edges[readEndEdge].length;

				ssize_t mateBeginEdge = g.paths[mateReadIndex][0].edge;
				ssize_t mateBeginIntv = g.paths[mateReadIndex][0].index;
				ssize_t mateBeginPos  = (*g.edges[mateBeginEdge].intervals)[mateBeginIntv].edgePos;

				ssize_t mateType = mateTable[readIndex].mateType;
				// The reads map to the same edge, definitely no problem.
				if (readEndEdge == mateBeginEdge) {
					continue;
				}
				
				if ((readEndEdgeLength - readEndPos) > rules[mateType].cloneLength + rules[mateType].cloneVar or
						(mateBeginPos) > rules[mateType].cloneLength + rules[mateType].cloneVar) {
					// The mate pairs are bad.  Remove the paths corresponding to the mate-pair
					// as well as the reverse complement.

					ssize_t pathsToRemove[4];
					
					pathsToRemove[0] = path;
					pathsToRemove[1] = mateReadIndex;
					pathsToRemove[2] = path + 1;
					pathsToRemove[3] = mateReadIndex - 1;
					
					// remove these 4 paths
					ssize_t r;
					for (r = 0; r < 4; r++) { 
						ssize_t pathToRemove = pathsToRemove[r];
						g.MarkPathForRemoval(pathToRemove);
						g.pathLengths[pathToRemove] = 0;
						delete[] g.paths[pathToRemove];
						g.paths[pathToRemove] = 0;
					}
					++numRemovedMatePairs;
				}
			}
		}
	}
	return numRemovedMatePairs;
}

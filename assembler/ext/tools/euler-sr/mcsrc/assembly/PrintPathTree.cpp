/***************************************************************************
 * Title:          PrintPathTree.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  09/09/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "IntervalGraph.h"
#include "PathLib.h"
#include "ReadPaths.h"

#include "IntegralTupleStatic.h"

int main(int argc, char* argv[]) {
	std::string graphFileName;
	if (argc <= 1) {
		std::cout << "usage: printPathTree graphName" << std::endl;
		exit(0);
	}

	graphFileName = argv[1];
	int argi = 2;
	ssize_t printPathSeq = 0;
	std::string pathSeqName = "";
	while (argi < argc ) {
		if (strcmp(argv[argi], "-pathSeq") == 0) {
			printPathSeq = 1;
			pathSeqName = argv[++argi];
		}
		++argi;
	}
	IntervalGraph graph;
	int vertexSize;
	ReadIntervalGraph(graphFileName, graph, vertexSize);

	PathIntervalList &paths = graph.paths;
	PathLengthList   &pathLengths = graph.pathLengths;
	TEdgeList &edges      = graph.edges;
	TVertexList &vertices = graph.vertices;
	//	ReadReadPaths(pathFileName, paths, pathLengths);

	PathBranch pathTree;
	CollectPathTree(paths, pathLengths, pathTree);
	//	PrintPathTree(pathTree);

	//	std::cout << "Done collecting paths.  Looking for redundant paths." << std::endl;
	// Make the list of paths more easy to search.
	PathTraceList pathTraces;
	PathTreeToPathList(pathTree, pathTraces);

	// Find paths that are subpaths of others.
	MarkEnclosedPaths(pathTraces);
	RemoveEnclosedPaths(pathTraces);

	ssize_t maxEdgeIndex = 0;
	ssize_t p, pi;
	for (p = 0; p < paths.size(); p++ ) {
		for (pi = 0; pi < pathLengths[p]; pi++) {
			if (maxEdgeIndex < paths[p][pi].edge) {
				maxEdgeIndex =paths[p][pi].edge; 
			}
		}
	}
	/*
		std::vector<ssize_t> intList, startList, endList;
	intList.resize(maxEdgeIndex+1);
	startList.resize(maxEdgeIndex+1);
	endList.resize(maxEdgeIndex+1);
	std::fill(intList.begin(), intList.end(), 0);
	std::fill(startList.begin(), startList.end(), 0);
	std::fill(endList.begin(), endList.end(), 0);
	*/
	//	MarkExInternal(pathTraces, startList, endList, intList);

	// Find edges that are both intList and startList. These
	// may not be resolved.
	ssize_t e;
	/*
	for (e =0 ; e < startList.size(); e++ ){ 
				if (startList[e] == 1 and intList[e] == 0) {
			std::cout << "possible src: " << e << std::endl;
		}
		if (endList[e] == 1 and intList[e] == 0) {
			std::cout << "possible dest: " << e << std::endl;
		}
	}
		*/
	
	TraceMapMatrix traceMaps;
	traceMaps.resize(edges.size());
	StoreTraceMaps(pathTraces, traceMaps);
	
	for (p = 0; p < pathTraces.size(); p++ ) {
		//UNUSED//		int traceEnd = pathTraces[p].edges->size() - 1;
		//UNUSED//		int firstEdge = (*pathTraces[p].edges)[0];
		//UNUSED//		int lastEdge =  (*pathTraces[p].edges)[traceEnd];
		ssize_t t;
		std::cout << p << " [ " << pathTraces[p].count << " ] ";
		for (t = 0; t < pathTraces[p].size(); t++ ) {
			std::cout << (*pathTraces[p].edges)[t] << " (" << edges[(*pathTraces[p].edges)[t]].length << ") ";
		}
		std::cout << std::endl;
	}
	return 0;
	/*
		if (traceMaps[firstEdge].numStart == 1 and 
				traceMaps[firstEdge].numInternal == 0 and 
				traceMaps[lastEdge].numEnd == 1 and
				traceMaps[lastEdge].numInternal == 0) {
			std::cout << "resolving path: " << p << " " << traceEnd + 1 << std::endl;
			if (AreInternalPathEdgesResolved(pathTraces, traceMaps, p)) {
				std::cout << "and the tangle is resolved as well." << std::endl;
			}
		}
	}
	*/
	return 0;
	for (e = 0; e < maxEdgeIndex + 1; e++ ) {
		std::cout << traceMaps[e].numStart << " " 
							<< traceMaps[e].numInternal << " " 
							<< traceMaps[e].numEnd << std::endl;
	}
	for (p = 0; p < pathTraces.size(); p++ ) {
		ssize_t traceEnd = pathTraces[p].edges->size() - 1;
		ssize_t firstEdge = (*pathTraces[p].edges)[0];
		ssize_t lastEdge =  (*pathTraces[p].edges)[traceEnd];
		ssize_t firstConsistentIndex;
		if (traceMaps[firstEdge].numStart == 1 and
				traceMaps[firstEdge].numInternal == 0 and 
				traceMaps[firstEdge].numEnd <= 1) {
			std::cout << " a resolving path of length: " << pathTraces[p].size() << " may begin from: " << firstEdge << " " 
								<< edges[firstEdge].length << " ";
			std::cout << vertices[edges[firstEdge].src].InDegree() << " " 
								<< vertices[edges[firstEdge].src].OutDegree() << std::endl;
			std::cout << " the traces that pass through this edge are: " << std::endl;
			ssize_t te;
			for (te = 0; te <= traceEnd; te++ ) {
				std::cout << (*pathTraces[p].edges)[te] << " (" 
									<< traceMaps[(*pathTraces[p].edges)[te]].numStart << ","
									<< traceMaps[(*pathTraces[p].edges)[te]].numInternal << ","
									<< traceMaps[(*pathTraces[p].edges)[te]].numEnd << ","
									<< edges[(*pathTraces[p].edges)[te]].length << ")";
			}
			std::cout << " " << edges[(*pathTraces[p].edges)[te-1]].length;
			std::cout << std::endl;
			//UNUSED// ssize_t t;
			ssize_t tracePos;
			// Look to see if traces that start on internal edges
			// are consistent with this trace
			//UNUSED// ssize_t edgeTrace;
			//UNUSED// ssize_t internalEdge;
			ssize_t firstConsistent = 0;
			ssize_t consistentTraceIndex;
			ssize_t numConsistent;
			for (tracePos = 1; tracePos < traceEnd ; tracePos++) {
				ssize_t curTraceEdge = (*pathTraces[p].edges)[tracePos];
				// Count the number of traces beginning in this edge that are consistent
				// with the trace beginning at firstEdge
				numConsistent = 
					CountForwardConsistent(pathTraces[p], tracePos, 2,
																 pathTraces, traceMaps,
																 curTraceEdge, consistentTraceIndex);
				
				if (numConsistent == 1 and tracePos == 1) {
					firstConsistent = 1;
					firstConsistentIndex = consistentTraceIndex;
				}
				ssize_t numStarting = 0;
				ssize_t tp;
				for (tp = 0; tp < traceMaps[curTraceEdge].traces.size(); tp++ ){
					if (traceMaps[curTraceEdge].traces[tp].pos == 0) {
						numStarting++;
					}
				}
				
				std::cout << "there are " << numConsistent << " / " << numStarting
									<< " paths starting in " 
									<< curTraceEdge << " consistent with path: " << p << std::endl;
									
			}
			ssize_t extNumber = 1;
			consistentTraceIndex = firstConsistentIndex;
			while (firstConsistent  ) {

				// March forward to see if there are more first consistent paths.
				// 0 <- consistentTracePos must begin in the edge 
				std::cout << "checking extension " << extNumber << " on path: ";
				++extNumber;
				ssize_t cp;
				ssize_t cpe;
				ssize_t prevTraceIndex= consistentTraceIndex;
				for (cp = 0; cp < pathTraces[consistentTraceIndex].size(); cp++ ) {
					cpe = (*pathTraces[consistentTraceIndex].edges)[cp];
					std::cout << cpe << "(" << traceMaps[cpe].numStart << "," << traceMaps[cpe].numInternal 
										<< "," << traceMaps[cpe].numEnd << "," << edges[cpe].length << ") ";
				}
				std::cout << std::endl;
				assert(pathTraces[consistentTraceIndex].size() > 0);
				numConsistent = 
					CountForwardConsistent(pathTraces[consistentTraceIndex], 1, 2,
																 pathTraces, traceMaps,
																 (*pathTraces[consistentTraceIndex].edges)[1], 
																 consistentTraceIndex);
				
				if (numConsistent == 1) {
					std::cout << "is consistent with 1 extension" << std::endl;
					firstConsistent = 1;
					ssize_t traceLength = pathTraces[prevTraceIndex].size();
					lastEdge = (*pathTraces[prevTraceIndex].edges)[traceLength-1];
				}
				else {
					std::cout << "there were " << numConsistent << " extensions to the extension." << std::endl;
					firstConsistent = 0;
				}
			}
		}

		/*		if (startList[firstEdge] == 1 and
				(endList[firstEdge] == 1 or
				 (endList[firstEdge] == 0 and 
					vertices[edges[firstEdge].src].InDegree() == 0)) and
				intList[firstEdge] == 0 and
				(startList[lastEdge] == 1 or
				 (startList[lastEdge] = 0 and
					vertices[edges[lastEdge].dest].OutDegree() == 0)) and
				endList[lastEdge] == 1 and
				intList[lastEdge] == 0) {
			std::cout << "path: " << p << " is resolving." << std::endl;
			for (pi = 0; pi <= traceEnd; pi++ ){ 
				std::cout << (*pathTraces[p].edges)[pi] << " ";
			}
			std::cout << std::endl;

			for (pi = 1; pi < traceEnd; pi++) {
				ssize_t intEdgeIndex = (*pathTraces[p].edges)[pi];
				std::cout << "internal edge: " << intEdgeIndex << " is in paths: ";
				ssize_t t;
				
				for (t = 0; t < traceMaps[intEdgeIndex].size(); t++ ) {
					std::cout << traceMaps[intEdgeIndex][t].trace << " ("
										<< traceMaps[intEdgeIndex][t].pos << ") ";
				}
				std::cout << std::endl;
			}
		}
		*/
	}
	/*
	for (e = 0; e < edges.size(); e++ ){ 
		if (IsEdgePathContained(pathTraces, traceMaps, e)) {
			std::cout << "Edge: " << e << " of length: " << edges[e].length << " is contained." << std::endl;
			std::set<EdgePair> edgePairs;
			
			if (AreEntranceAndExitEdgesPaired(pathTraces, traceMaps, vertices, edges, e, edgePairs)) {
				std::cout << "And they are paired!" << std::endl;
			}
		}
	}
	*/

	std::ofstream traceOut;
	openck("path.traces", traceOut, std::ios::out);
	PrintPathTraces(pathTraces, graph, traceOut);
	return 0;
}

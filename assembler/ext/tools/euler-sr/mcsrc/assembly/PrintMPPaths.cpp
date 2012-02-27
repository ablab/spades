/***************************************************************************
 * Title:          PrintMPPaths.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  08/13/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "IntervalGraph.h"
#include "PathLib.h"
#include "MateLibrary.h"



int main(int argc, char* argv[]) {
	std::string graphFileName, matePairOutName;
	int vertexSize;
	if (argc <= 3) {
		std::cout << "usage: mateTransformGraph graphName vertexSize mateTableName ruleFileName matePairPathname" << std::endl;
		exit(0);
	}
	std::string mateTableName, ruleFileName;

	graphFileName   = argv[1];
	vertexSize      = atoi(argv[2]);
	mateTableName   = argv[3];
	ruleFileName    = argv[4];
	matePairOutName = argv[5];

	std::ofstream matePairOut;
	openck(matePairOutName, matePairOut, std::ios::out);
	IntervalGraph graph;
	ReadIntervalGraph(graphFileName, graph, vertexSize);
	RuleList rules;
	ParseRuleFile(ruleFileName, rules);

	ReadMateList mateList;
	ReadMateTable(mateTableName, mateList);

	PathIntervalList &paths       = graph.paths;
	PathLengthList   &pathLengths = graph.pathLengths;
	TEdgeList        &edges       = graph.edges;
	TVertexList      &vertices    = graph.vertices;

	// Go from the mate paths to regular paths.

	PathIntervalList matePaths;
	PathLengthList   matePathLengths;
	ssize_t p, readIndex;

	for (p = 0; p < graph.paths.size(); p+=2 ) {
		readIndex = p/2;
		if (mateList[readIndex].mateIndex == -1)
			continue;

		if (graph.pathLengths[p] == 0 )
			continue;

		ssize_t mateType, mateIndex, lastPathIntv, lastEdge, lastEdgeIntv,
			lastEdgePos;

		mateIndex = mateList[readIndex].mateIndex;
		mateType  = mateList[readIndex].mateType;
		
		lastPathIntv = graph.pathLengths[p] - 1;
		lastEdge = graph.paths[p][lastPathIntv].edge;
		lastEdgeIntv = graph.paths[p][lastPathIntv].index;

		lastEdgePos = (*graph.edges[lastEdge].intervals)[lastEdgeIntv].edgePos + 
			(*graph.edges[lastEdge].intervals)[lastEdgeIntv].length;

		// use the reverse complment path of the mate since
		// the mate is sequenced in the opposite direction
		ssize_t mp;
		mp = mateIndex * 2 + 1;
		
		if (graph.pathLengths[mp] == 0) 
			continue;
		ssize_t mateEdge, mateEdgeIntv, mateStartEdgePos;

		mateEdge     = graph.paths[mp][0].edge;
		mateEdgeIntv = graph.paths[mp][0].index;

		mateStartEdgePos = (*graph.edges[mateEdge].intervals)[mateEdgeIntv].edgePos +
			(*graph.edges[mateEdge].intervals)[mateEdgeIntv].length;

		ssize_t foundAValidPath;

		// This could take much more time.
		ssize_t numValidPaths = 0;
		ssize_t foundValidPaths = 0;
		MatePathList matePath;
		ssize_t storeMatePath = 1;
		foundValidPaths = CountValidMatePaths(graph,
																					lastEdge, lastEdgePos, // cur edge pos
																					mateEdge, mateStartEdgePos,
																					0,  // starting length is 0
																					rules[mateType].cloneLength - rules[mateType].cloneVar, // lower end
																					rules[mateType].cloneLength + rules[mateType].cloneVar, // upper end
																					30, // max depth to search
																					4, numValidPaths, 
																					storeMatePath, matePath, 0);

		if (numValidPaths == 1 and (matePath.size() > 0 or lastEdge != mateEdge) ) {
			MatePathList::iterator pathIt;
			PathInterval *path;
			ssize_t matePathLength = matePath.size() + 1;
			path = new PathInterval[matePathLength];
			std::cout << "mate path: ";
			path[matePath.size()].edge = mateEdge;
			path[matePath.size()].index = -1;
			ssize_t pathPos = 0;
			for (pathIt = matePath.begin(); pathIt != matePath.end(); ++ pathIt) {
				path[pathPos].edge = (*pathIt).edge;
				path[pathPos].index = -1;
				std::cout << " " << (*pathIt).edge << " (" << edges[(*pathIt).edge].length << "), ";
				++pathPos;
			}
			std::cout << " " << mateEdge << " (" << edges[mateEdge].length << ")" << std::endl;
			matePaths.push_back(path);
			matePathLengths.push_back(matePathLength);


			// Add the balance of this path.
			PathInterval *rcPath;
			rcPath = new PathInterval[matePathLength];
			ssize_t i;
			for(i = 0; i < matePathLength; i++ ){ 
				rcPath[i].edge = edges[path[matePathLength - i - 1].edge].balancedEdge;
				rcPath[i].index = -1;
			}
			matePaths.push_back(rcPath);
			matePathLengths.push_back(matePathLength);
		}
	}
		

	

	PathBranch pathTree;
	CollectPathTree(matePaths, matePathLengths, pathTree);
	//	PrintPathTree(pathTree);

	std::cout << "Done collecting paths.  Looking for redundant paths." << std::endl;
	// Make the list of paths more easy to search.
	PathTraceList pathTraces;
	PathTreeToPathList(pathTree, pathTraces);

	// Find paths that are subpaths of others.
	MarkEnclosedPaths(pathTraces);
	RemoveEnclosedPaths(pathTraces);
	
	TraceMapMatrix traceMaps;
	traceMaps.resize(edges.size());
	StoreTraceMaps(pathTraces, traceMaps);


	MarkResolvedPaths(pathTraces, traceMaps);

	// Make sure each path has it's own balance, but don't bother 
	// doing this in debug mode.
	assert(CheckPathBalance(edges, pathTraces, traceMaps));

	// Check the balance of all traces.
	for (p = 0; p < pathTraces.size(); p++ ){ 
		ssize_t traceEnd  = pathTraces[p].edges->size() - 1;
		ssize_t firstEdge = (*pathTraces[p].edges)[0];
		ssize_t lastEdge  = (*pathTraces[p].edges)[traceEnd];
		ssize_t traceLength = (*pathTraces[p].edges).size();
		ssize_t balEdge   = edges[lastEdge].balancedEdge;
		
		// There must exist a path that starts in the balance of the last edge and is the exact
		// balance of this path.
		ssize_t balPathFound = 0;
		ssize_t t;
		for (t = 0; t < traceMaps[balEdge].traces.size() and !balPathFound; t++) {
			ssize_t bpi, bp;
			
			bp = traceMaps[balEdge].traces[t].trace;
			if (pathTraces[bp].edges->size() == traceLength) {
				ssize_t pathsAgree = 1;
				for (bpi = 0; bpi < pathTraces[bp].edges->size() and pathsAgree; bpi++ ) {
					if (edges[(*pathTraces[bp].edges)[bpi]].balancedEdge != (*pathTraces[p].edges)[traceLength - bpi - 1])
						pathsAgree = 0;
				}
				if (pathsAgree)
					balPathFound = 1;
			}
		}
		assert(balPathFound);
	}

	
	// Now check the balance of the resolved traces.
	ssize_t e;
	ssize_t resolvedTraceIndex;
	ssize_t endEdgeIndex, beginEdgeIndex;
	ssize_t traceMapIndex;
	for (e = 0; e < edges.size(); e++ ) {  
		if (traceMaps[e].resolved) {
			// Edge 'e' begins a resolved path. 
			std::cout << e << " is resolved. " << std::endl;
			traceMapIndex  = traceMaps[e].GetFirstStartTraceIndex();
			resolvedTraceIndex  = traceMaps[e].traces[traceMapIndex].trace;
			std::cout << "resolving trace: " << resolvedTraceIndex << std::endl;
			beginEdgeIndex = (*pathTraces[resolvedTraceIndex].edges)[0];
			endEdgeIndex   = (*pathTraces[resolvedTraceIndex].edges)[pathTraces[resolvedTraceIndex].edges->size() - 1];
			
			// The path that starts at the balanced edge should be resolved.
			ssize_t balEdge = edges[endEdgeIndex].balancedEdge;
			std::cout << "bal edge: " << balEdge << std::endl;
			if (!traceMaps[balEdge].resolved) {
				std::cout << "ERROR!  Edge " << beginEdgeIndex << " ends a trace but its balance " 
									<< balEdge << " does not begin one. " << beginEdgeIndex
									<< std::endl;
				assert(0);
			}
		}
	}
	// Now fix the balance of the paths.

 	for (e = 0; e < edges.size(); e++ ) {  
		if (traceMaps[e].resolved) {
			traceMapIndex       = traceMaps[e].GetFirstStartTraceIndex();
			resolvedTraceIndex  = traceMaps[e].traces[traceMapIndex].trace;
			//			assert((*pathTraces[resolvedTraceIndex].edges)[0] == traceMaps[e].startEdge);
			DetachPath(pathTraces, traceMaps, graph, vertices, edges,
								 paths, pathLengths,
								 pathTraces[resolvedTraceIndex], resolvedTraceIndex);
		}
	}

	graph.RemoveMarkedIntervalsNoPaths();
	graph.RemoveEmptyEdges();
	ssize_t nremoved = graph.RemoveEmptyVertices();
	std::cout << "removed: " << nremoved << " vertices." << std::endl;
	graph.CondenseSimplePaths();
	ssize_t pi;
	/*
	for (p = 0; p < paths.size(); p++) { 
		for (pi = 0; pi < pathLengths[p]; pi++) { 
			assert((*edges[paths[p][pi].edge].intervals)[paths[p][pi].index].read == p);
		}
		if (p % 2 == 0)
			assert(pathLengths[p] == pathLengths[p+1]);
	}
	*/
	std::string bGraphOutName = graphOutName + ".bgraph";
	std::string intvOutName = graphOutName + ".intv";
	std::string gvzOutName  = graphOutName + ".dot";
	std::string pathOutName = graphOutName + ".path";
	std::string edgeOutName = graphOutName + ".edge";
	std::string euGraphName = graphOutName + ".graph";
	CheckEdges(graph.vertices, graph.edges);

	graph.PrintIntervalGraph(bGraphOutName, intvOutName);


  PrintGraph(graph.vertices,graph.edges, euGraphOutName);
  PrintEdges(graph.vertices, graph.edges, edgeOutName);
  GVZPrintBGraph(graph.vertices, graph.edges, gvzOutName);
	WriteReadPaths(pathOutName, graph.paths, graph.pathLengths);

	return 0;

}



/***************************************************************************
 * Title:          TestFindMatePath.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  08/13/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "IntervalGraph.h"
#include "MateLibrary.h"
#include "PathLib.h"
#include <vector>
#include <iterator>
#include "IntegralTupleStatic.h"

void PrintUsage() {
	std::cout << " usage: testFindMatePath graphIn mateTable " << std::endl
						<< "  [-minPathCount c]  Remove paths with count less than 'c' " << std::endl
						<< "  [-minMatePairCount c] Remove paths with mate count less than 'c'" << std::endl;
}




int main(int argc, char* argv[]) {
	
	std::string graphFileName, mateTableName;
	if (argc < 3) {
		PrintUsage();
		exit(0);
	}
	graphFileName = argv[1];
	mateTableName = argv[2];

	ssize_t minPathCount = 2;
	ssize_t minMatePairCount = 2;
	int argi = 3;
	while (argi < argc) {
		if (strcmp(argv[argi], "-minPathCount") == 0) {
			minPathCount = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-minMatePairCount") == 0) {
			minMatePairCount = atoi(argv[++argi]);
		}
		else {
			PrintUsage();
			std::cout << "bad option: " << argv[argi] << std::endl;
		}
		++argi;
	}

	IntervalGraph graph;
	ReadMateList  mateTable;
	std::cout << "Reading interval graph." << std::endl;
	int vertexSize;
	ReadIntervalGraph(graphFileName, graph, vertexSize);
	std::cout << "Reading mate table." << std::endl;
	ReadMateTable(mateTableName, mateTable);

	PathIntervalList &paths       = graph.paths;
	PathLengthList   &pathLengths = graph.pathLengths;
	TEdgeList        &edges       = graph.edges;
	//UNUSED// TVertexList      &vertices    = graph.vertices;


	/*
	 * Use read-paths to determine which edges are unique or not.
	 */

	PathBranch pathTree, removedPathTree;

	CollectPathTree(paths, pathLengths, pathTree);

	
	PathTraceList pathTraces;
	PathTreeToPathList(pathTree, pathTraces);

	MarkEnclosedPaths(pathTraces);
	RemoveEnclosedPaths(pathTraces);

	TraceMapMatrix traceMaps;
	traceMaps.resize(edges.size());
	StoreTraceMaps(pathTraces, traceMaps);

	ssize_t e;
	for (e = 0; e < edges.size(); e++ ) {
		MateEdgeMap mateEdges;
		std::list<ssize_t> srcEdgePath;
		std::cout << "searching for mates from: " << e << std::endl;
		CollectMateEdges(graph, mateTable, e, mateEdges, -1);
		RemoveLowCountMateEdges(mateEdges, 3);
		if (FindMatePath(graph, e, mateEdges, srcEdgePath)) {
			std::cout << "Success!" << std::endl;
		}
	}

	/*	
	std::string bGraphOutName = graphOutName + ".bgraph";
	std::string intvOutName = graphOutName + ".intv";
	std::string gvzOutName  = graphOutName + ".dot";
	std::string pathOutName = graphOutName + ".path";
	std::string edgeOutName = graphOutName + ".edge";
	std::string euGraphOutName = graphOutName + ".graph";
	
	graph.PrintIntervalGraph(bGraphOutName, intvOutName);
  PrintGraph(graph.vertices,graph.edges, euGraphOutName);
  PrintEdges(graph.vertices, graph.edges, edgeOutName);
  GVZPrintBGraph(graph.vertices, graph.edges, gvzOutName);
	WriteReadPaths(pathOutName, graph.paths, graph.pathLengths);
	*/

}

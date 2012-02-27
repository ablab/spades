/***************************************************************************
 * Title:          TransformGraph.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  10/31/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "IntervalGraph.h"
#include "PathLib.h"
#include "MateLibrary.h"
#include "IntegralTupleStatic.h"

using namespace std;

void PrintUsage() {
	cout << "usage: transformGraph graphName graphOutName " << endl;
	cout << "  -notStrict   (FALSE) Resolve any particular instance of a repeat even if not all " << endl
						<< "               instances of that repeat are resolved." << endl << endl
						<< "  -minPathCount c(0) " << endl
						<< "               Remove paths if they occur less than c times." << endl << endl
						<< "  -erodeShortPaths E Erode the ends of paths if they are less than 'E'"<<endl<<endl
						<< "               into an edge." << endl << endl;
	
}

int main(int argc, char* argv[]) {
	string graphFileName, graphOutName;
	int vertexSize;
	if (argc <= 2) {
		PrintUsage();
		exit(0);
	}

	graphFileName = argv[1];
	graphOutName = argv[2];
	int argi= 3;
	ssize_t notStrict = 0;
	ssize_t minPathCount = 0;
	ssize_t erodeShortPaths = 0;
	while (argi < argc) {
		if (strcmp(argv[argi], "-notStrict") == 0) {
			notStrict = 1;
		}
		else if (strcmp(argv[argi], "-minPathCount") == 0) {
			minPathCount = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-erodeShortPaths") == 0 ){
			erodeShortPaths = atoi(argv[++argi]);
		}
		else {
			PrintUsage();
			exit(0);
		}
		++argi;
	}

	std::string reportFileName = FormReportName(graphFileName);
	ofstream report;
	openck(reportFileName, report, std::ios::app, report);
	BeginReport(argc, argv, report);


	IntervalGraph graph;
	ReadIntervalGraph(graphFileName, graph, vertexSize, report);
	cout << "got " << graph.edges.size() << " edges." << endl;

	PathIntervalList &paths       = graph.paths;
	PathLengthList   &pathLengths = graph.pathLengths;
	TEdgeList        &edges       = graph.edges;
	TVertexList      &vertices    = graph.vertices;
	//	ReadReadPaths(pathFileName, paths, pathLengths, report);


	ssize_t changeMade;

	do {
		PathBranch pathTree, removedPathTree;

	CollectPathTree(paths, pathLengths, pathTree);

	cout << "Done collecting paths.  Looking for redundant paths." << endl;
	// Make the list of paths more easy to search.
	list<ssize_t> curPath;


	TrimLowCoverageBranches(pathTree, minPathCount);

	// PrintPathTraceList(pathTraces);
	// Find paths that are subpaths of others.


	PathTraceList pathTraces;
	PathTreeToPathList(pathTree, pathTraces);
	PrintPathTraces(pathTraces, graph, cout);

	MarkEnclosedPaths(pathTraces);
	RemoveEnclosedPaths(pathTraces);
	//CheckPathTraceListBalance(graph.edges, pathTraces);
	cout << "removing low count paths." << endl;
	RemoveLowCountPaths(graph, paths, pathLengths, pathTree, minPathCount);
	//	PrintPathTraces(pathTraces, graph, cout);
	cout << "removing low count traces. " << endl;
	RemoveLowCountTraces(pathTraces, minPathCount, removedPathTree);
	//	#ifndef NDEBUG
		cout << "removed paths: " << endl;
		PrintPathTraces(pathTraces, graph, cout);
		cout << "done." << endl;
		//#endif
		//	assert(CheckPathTraceListBalance(edges, pathTraces));
		//		CheckPathBalance(pathTraces, tra

	//UNUSED+// ssize_t pi;
	ssize_t p ;

  assert(graph.CheckAllPathsBalance(1));

		changeMade = 0;
		TraceMapMatrix traceMaps;
		traceMaps.resize(edges.size());
		StoreTraceMaps(pathTraces, traceMaps);
	
		MarkResolvedPaths(pathTraces, traceMaps, notStrict);
		//		PrintCandidatePaths(pathTraces, traceMaps, graph.edges);

		// Make sure each path has it's own balance, but don't bother 
		// doing this in debug mode.
		assert(CheckPathBalance(edges, pathTraces, traceMaps));

	
		// Now check the balance of the resolved traces.
		ssize_t e;
		ssize_t resolvedTraceIndex;
		//UNUSED// ssize_t endEdgeIndex, beginEdgeIndex;
		ssize_t traceMapIndex;

		vector<ssize_t> balancedEdges;
		balancedEdges.resize(edges.size());
		for (e = 0; e < edges.size(); e++ ) {  
			balancedEdges[e] = edges[e].balancedEdge;
		}
		ssize_t numResolvedRepeats = 0;
		for (p = 0; p < pathTraces.size(); p++ ){
			ssize_t firstEdge = (*pathTraces[p].edges)[0];
			ssize_t traceLength = pathTraces[p].edges->size();
			ssize_t lastEdge = (*pathTraces[p].edges)[traceLength -1];
			ssize_t balancedEdge = balancedEdges[lastEdge];
			//UNUSED// ssize_t balancedLastEdge = balancedEdges[firstEdge];

			if (traceMaps[firstEdge].outResolved and traceMaps[lastEdge].inResolved) {

				ssize_t resolvedLength = (*pathTraces[p].edges).size();
				ssize_t lastEdge = (*pathTraces[p].edges)[resolvedLength-1];
				ssize_t resolvedTraceEndBal = balancedEdges[(*pathTraces[p].edges)[resolvedLength-1]];
				ssize_t resolvedTraceBeginBal = balancedEdges[(*pathTraces[p].edges)[0]];
				//			assert((*pathTraces[resolvedTraceIndex].edges)[0] == traceMaps[e].startEdge);
				ssize_t b = resolvedTraceEndBal;
				assert (traceMaps[b].outResolved);
				//					cout << "detaching bal  starting at: " << b << endl; 
				traceMapIndex       = traceMaps[b].GetFirstStartTraceIndex();
				resolvedTraceIndex  = traceMaps[b].traces[traceMapIndex].trace;

				
				if (resolvedTraceBeginBal != firstEdge and
						resolvedTraceBeginBal != lastEdge and
						resolvedTraceEndBal != firstEdge and 
						resolvedTraceEndBal != lastEdge)  {

				ssize_t pe;
				cout << "detaching path " << p << " starting at: " 
						 << firstEdge << " ["  << edges[firstEdge].index << "]" << " ";
				

				for (pe = 0; pe < pathTraces[p].edges->size(); pe++) {
					cout << (*pathTraces[p].edges)[pe] << " ";
				}
				cout << endl;

				cout << "detaching (bal) path " << resolvedTraceIndex << " starting at: " 
						 <<  balancedEdge << " [" << edges[balancedEdge].index << "] " << " ";
				for (pe = 0; pe < pathTraces[resolvedTraceIndex].edges->size(); pe++) {
					cout << (*pathTraces[resolvedTraceIndex].edges)[pe] << " ";
				}
				cout << endl;
				



				assert(pathTraces[p].edges->size() == pathTraces[resolvedTraceIndex].edges->size());

				for (pe = 0; pe < pathTraces[p].edges->size(); pe++) {
					assert(edges[(*pathTraces[p].edges)[pe]].balancedEdge == 
								 (*pathTraces[resolvedTraceIndex].edges)[resolvedLength - 1 - pe]);
				}

					//					assert(graph.CheckAllPathsBalance(1));
					changeMade = 1;
					BalancedDetachPaths(pathTraces, traceMaps, graph, vertices, edges,
															paths, pathLengths,
															pathTraces[p], p,
															pathTraces[resolvedTraceIndex], resolvedTraceIndex);
					
					traceMaps[firstEdge].outResolved = 0;
					traceMaps[lastEdge].inResolved = 0;
					traceMaps[b].outResolved = 0;
					traceMaps[resolvedTraceBeginBal].inResolved = 0;
					++numResolvedRepeats;

					//					assert(graph.CheckAllPathsBalance(1));
				}
			}
		}


		// Somehow some dirt ends up left over. Clean it here.
		for (e = 0; e < edges.size(); e++ ){ 
			if (edges[e].src == -1 and edges[e].dest == -1 and edges[e].intervals->size() > 0)
				edges[e].intervals->clear();
		}
			
		graph.RemoveMarkedIntervalsNoPaths();
		//		graph.RemoveEmptyEdges();
		graph.RemoveUnlinkedEdges();
		ssize_t nremoved = graph.RemoveEmptyVertices();
		cout << "resolved: " << numResolvedRepeats << " repeats and removed: " << nremoved << " vertices." << endl;
		graph.CondenseSimplePaths();

		if (!changeMade and erodeShortPaths) {
			graph.ErodeShortEndIntervals(erodeShortPaths);
			changeMade = 1;
			erodeShortPaths = 0;
		}
	}
	while (changeMade);
	//	while (0);

	string bGraphOutName = graphOutName + ".bgraph";
	string intvOutName = graphOutName + ".intv";
	string gvzOutName  = graphOutName + ".dot";
	string pathOutName = graphOutName + ".path";
	string edgeOutName = graphOutName + ".edge";
	string euGraphOutName = graphOutName + ".graph";
	CheckEdges(graph.vertices, graph.edges);

	graph.CondenseEdgeLists();
	graph.SetMultiplicities();
	graph.PrintIntervalGraph(bGraphOutName, intvOutName, report);


  PrintGraph(graph.vertices,graph.edges, euGraphOutName, report);
  PrintEdges(graph.vertices, graph.edges, edgeOutName, report);
  GVZPrintBGraph(graph.vertices, graph.edges, gvzOutName, report);
	WriteReadPaths(pathOutName, graph.paths, graph.pathLengths, report);

	EndReport(report);
	report.close();

	return 0;

}

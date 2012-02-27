/***************************************************************************
 * Title:          MateTransformGraph.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  11/27/2009
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
#include "Scaffold.h"


using namespace std;
void PrintUsage() {
	cout << "usage: mateTransformGraph graphName  "
			 << "mateTableName ruleFileName graphOutName  " << endl
			 << "      [-minMatePairCount m]  Remove mate pairs if less than m" << endl 
			 << "                             pairs confirm an edge-pair." 
			 << endl << endl
			 << "      [-notStrict]           Resolve repeats even when not all " << endl 
			 << "                             paths through a repeat are resolved." 
			 << endl << endl
			 << "      [-findPaths]           Find tough-to trace paths between paired edges." 
			 << endl << endl
			 << "      [-noIter]  Do not iteratively simplify paths."
			 << endl << endl 
			 << "      [-minScafCov N (40)]   Consider edges valid for scaffolding if" << endl 
			 << "                             at least N reads were used to build the edge." 
			 <<endl << endl
			 << "      [-scaffoldEdgeLength L (100)]" << endl
			 << "                             Only try and resolve repeats between edges of length L or greater."
			 << endl << endl
			 << "      [-onlyMatePaths]       Only resolve repeats using mate-paths, not read-paths." 
			 << endl << endl
			 << "      [-startPhase p]        Start on phase 1 (find paths), 2 (support), 3" <<endl
			 << "                             mate paths only, or 4, skip to scaffold" << endl
			 << endl << endl
			 << "      [-ruleType T ]         Only transform mates of clone type T." 
			 << endl << endl
			 << "      [-pathReads file]      Print paths as reads to 'file'." 
			 << endl << endl
			 << "      [-readLength l]        Use 'l' as an estimated read length (for printing paths)."
			 << endl << endl
			 << "      [-scaffold]            Use mate-pairs to scaffold disconnected contigs." << endl;
}

int main(int argc, char* argv[]) {
	string graphFileName, graphOutName;
	int vertexSize;
	if (argc < 5) {
		PrintUsage();
		exit(0);
	}
	string mateTableName, ruleFileName;
	string matePathSeqName;

	graphFileName   = argv[1];
	mateTableName   = argv[2];
	ruleFileName    = argv[3];
	graphOutName    = argv[4];

	matePathSeqName = graphOutName + ".pathseq";
	int argi = 5;
	ssize_t minMatePairCount = 3;
	ssize_t notStrict = 0;
	ssize_t verbose = 1; // default to verbose for now.
	ssize_t iterate = 1;
	ssize_t findPaths = 0;
	ssize_t uniquePathsOnly = 1;
	ssize_t minScaffoldCoverage = 40;
	//UNUSED// ssize_t onlyMatePaths = 0;
	ssize_t ignoreInconsistentPaths = 0;
	ssize_t scaffoldEdgeLength = 100;
	ssize_t exactPathsPhase      = 1;
	ssize_t scaffoldPathsPhase   = 2;
	ssize_t pairedPathsOnlyPhase = 3;
	ssize_t runScaffold  = 0;
	ssize_t currentPhase = 1;
	ssize_t endPhase     = 3;
	ssize_t ruleType     = 1;
	string pathReadsName = "";
	ssize_t readLength = 35;
	ssize_t startPhase = 1;
	while (argi < argc) {
		if (strcmp(argv[argi], "-minMatePairCount") == 0) {
			minMatePairCount = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-notStrict") == 0) {
			notStrict = 1;
		}
		else if (strcmp(argv[argi], "-verbose") == 0) {
			verbose = !verbose;
		}
		else if (strcmp(argv[argi], "-noIter") == 0) {
			iterate = 0;
		}
		else if (strcmp(argv[argi], "-findPaths") == 0) {
			findPaths = 1;
			uniquePathsOnly = 0;
			currentPhase = scaffoldPathsPhase;
		}
		else if (strcmp(argv[argi], "-minScafCov") == 0) {
			minScaffoldCoverage = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-onlyMatePaths") == 0) {
			//			startPhase = pairedPathsOnlyPhase;
			// onlyMatePaths = 1;
			endPhase = 4;
		}
		else if (strcmp(argv[argi], "-startPhase") == 0) {
			startPhase = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-ruleType") == 0) {
			ruleType = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-pathReads") == 0) {
			pathReadsName = argv[++argi];
		}
		else if (strcmp(argv[argi], "-readLength") == 0) {
			readLength = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-scaffoldEdgeLength") == 0) {
			scaffoldEdgeLength = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-scaffold") == 0) {
			runScaffold = 1;
		}
		else {
			PrintUsage();
			cout << "bad option: " << argv[argi] << endl;
			exit(1);
		}
		++argi;
	}
		
	string reportFileName = FormReportName(graphFileName);
	std::ofstream report;
	openck(reportFileName, report, std::ios::app, report);
	BeginReport(argc, argv, report);

	EdgePairMap pairedEdgeMap;

	RuleList rules;
	ParseRuleFile(ruleFileName, rules, report);
	
	IntervalGraph graph;
	ReadIntervalGraph(graphFileName, graph, vertexSize, report);
	vertexSize = graph.vertexSize;

	ReadMateList mateList;
	ReadMateTable(mateTableName, mateList, report);


	PathIntervalList &paths       = graph.paths;
	PathLengthList   &pathLengths = graph.pathLengths;
	TEdgeList        &edges       = graph.edges;
	TVertexList      &vertices    = graph.vertices;

	// Go from the mate paths to regular paths.

	//UNUSED+// ssize_t readIndex;
	ssize_t p ;
	PathBranch pathTree;
	
	// Compute some statistics on the number of mate-pairs spanning 
	// each pair of edges.
	//	int mateType, mateIndex, lastPathIntv, 
	//		lastEdge,	lastEdgeIntv,	lastEdgePos;

	ssize_t graphIsTransformed;
	ssize_t detachIter = 0;

	ssize_t numBadMates;
	ssize_t mateSep;
	double mateStddev;

	ssize_t curRule;

	for (curRule = 0; curRule < rules.size(); curRule++) {
		ComputeMatePairLengthDistribution(graph, mateList, curRule, mateSep, mateStddev); 
		currentPhase = startPhase;
		cout << "Computed mated pair length distro: " << mateSep << " " << mateStddev << endl;

		do {
		

			numBadMates = RemoveBadStartEndMatePairs(graph, mateList, rules, curRule);
			cout << "Removed " << numBadMates << " bad mate pairs." << endl;

			//UNUSED// ssize_t numEmptyEdges = 0;
			//			assert(graph.CheckAllPathsBalance(1));
			assert(graph.CheckBalance());
			// cout << "removing low coverage edges." << endl;
			// graph.RemoveLowCoverageEdges(5,3);
			ssize_t nCondensed = 0;
			//		nCondensed = graph.CondenseSimplePaths();
			cout << "after removing bad mate pairs, and the reads corresponding to them "
					 << nCondensed << " edges were condensed." << endl;

			assert(CheckEdges(vertices, edges));
		
			PathTraceList mateTraces;
		
			graphIsTransformed = 0;

			if (currentPhase == exactPathsPhase) {
				cout << "Exact paths phase." << endl;
				// Remove mate-pairs that are too low in frequency (perhaps remove the reads as well)
			
				// Count the mate pais on every edge, and store their average position.
				pairedEdgeMap.clear();
				StoreEdgePairMap(graph, mateList, pairedEdgeMap, curRule);
				/*
				EdgePairMap::iterator peIt, peEnd;
				for (peIt = pairedEdgeMap.begin(); peIt != pairedEdgeMap.end(); ++peIt) {
					cout << (*peIt).first.edge1 << " " << (*peIt).first.edge2 << " " << (*peIt).second.count << endl;
				}
				*/
				// If two edges are joined by mate pairs less than minMatePairCount
				// get rid of the mate pairings.
				RemoveLowFrequencyEdgePairs(graph, pairedEdgeMap, mateList, minMatePairCount, curRule);
			
				// Search for paths between mates, when a single path exists, store it.
				StoreUniqueMatePaths(graph, vertices, edges, rules, pairedEdgeMap, pathTree);

			}
			else if (currentPhase == scaffoldPathsPhase or 
							 currentPhase == pairedPathsOnlyPhase) {
				if (currentPhase == scaffoldPathsPhase) {
					cout << "Supported paths phase." << endl;
				}
				else if (currentPhase == pairedPathsOnlyPhase) {
					cout << "Supported paths phase. Ignoring read paths." << endl;
				}
				ignoreInconsistentPaths = 1;
				if (currentPhase == pairedPathsOnlyPhase) {
					graph.ErodeShortEndIntervals(5);
				}

				ssize_t e;
				for (e = 0; e < graph.edges.size(); e++ ){
					MateEdgeMap::iterator mit, mpeIt;
					MateEdgeMap mateEdges;
					CollectMateEdges(graph, mateList, e, mateEdges, 0); // collect mate-pairs from the src edge
					/*				cout << "mates for edge: " << e << " " << graph.edges[e].length << " ";
					for (mit = mateEdges.begin(); mit != mateEdges.end(); ++mit) {
												cout << (*mit).first  << "(" << graph.edges[(*mit).first].length << "," << (*mit).second.count << ") ";
					}
					cout << endl;
					*/
				}

				// Mark edges that have a minimal coverage.
				MarkScaffoldEdges(graph, minScaffoldCoverage, scaffoldEdgeLength);
			
				PathIntervalList pairedMatePathList;
				PathLengthList   pairedMatePathLengthList;
				FindScaffoldPaths(graph, mateList, 
													rules, curRule, scaffoldEdgeLength,
													pairedMatePathList, pairedMatePathLengthList);

				CollectPathTree(pairedMatePathList, pairedMatePathLengthList, pathTree);
				//				cout << "the mate traces are:" << endl;
				PathTraceList mateTraces;
				PathTreeToPathList(pathTree, mateTraces);
				TraceMapMatrix mateTraceMaps;
				mateTraceMaps.resize(edges.size());
				StoreTraceMaps(mateTraces, mateTraceMaps);
				//				PrintPathTraceResolution(edges, mateTraces, mateTraceMaps);
				//				PrintPathTraces(mateTraces, graph, cout);				
				/*
					if (pathReadsName != "") {
					PrintTracesAsReads(graph.vertices, graph.edges, mateTraces, readLength, pathReadsName);
					pathReadsName = "";
					}
				*/
			} // End checking paired paths only
			else { 
				break;
			}


			// Collect the paths due to mate-pairs.
			PathTreeToPathList(pathTree, mateTraces);
			if (!(currentPhase != pairedPathsOnlyPhase) ) {
				// Collect the paths due to reads so that mate transformations
				// are consistent with read paths.
				CollectPathTree(paths, pathLengths, pathTree);
				//	PrintPathTree(pathTree);
		
				cout << "Done collecting paths.  Looking for redundant paths." << endl;
				// Make the list of paths more easy to search.
			}
		
			PathTraceList pathTraces;
			PathTreeToPathList(pathTree, pathTraces);
			// Find paths that are subpaths of others.
			ssize_t pt;
			for (pt = 0; pt < pathTraces.size(); pt++ ){ 
				assert(pathTraces[pt].edges->size() >0);
			}
			MarkEnclosedPaths(pathTraces);
			RemoveEnclosedPaths(pathTraces);
			for (pt = 0; pt < pathTraces.size(); pt++ ){ 
				assert(pathTraces[pt].edges->size() >0);
			}


	
			TraceMapMatrix traceMaps;
			traceMaps.resize(edges.size());
			StoreTraceMaps(pathTraces, traceMaps);

			PrintCandidatePaths(pathTraces, traceMaps, edges);

			MarkResolvedPaths(pathTraces, traceMaps, notStrict);
			CountAdjacencies(pathTraces, traceMaps);
			PrintTracesAsReads(vertices, edges, pathTraces, readLength, matePathSeqName, report);

			// Print some useful information about the paths
			if (verbose) PrintPathTraceResolution(edges, pathTraces, traceMaps);

			// Make sure each path has it's own balance, but don't bother 
			// doing this in optimized mode.
			assert(CheckPathBalance(edges, pathTraces, traceMaps));

			// Store the balanced edge of all paths.  This is necessary 
			// in order to update balance while the graph is being modified.
			vector<ssize_t> balancedEdges;
			balancedEdges.resize(edges.size());
			ssize_t e;
			for (e = 0; e < edges.size(); e++) {
				balancedEdges[e] = edges[e].balancedEdge;
			}

			// Just for some bookkeeping and stats later on.
			set<ssize_t> internalDetachedEdges;
			// Now fix the balance of the paths.
			for (p = 0; p < pathTraces.size(); p++ ) {
				ssize_t firstEdge = (*pathTraces[p].edges)[0];
				ssize_t traceLength = pathTraces[p].edges->size();
				ssize_t lastEdge = (*pathTraces[p].edges)[traceLength -1];
				//UNUSED// ssize_t balancedEdge = balancedEdges[lastEdge];
				//UNUSED// ssize_t balancedLastEdge = balancedEdges[firstEdge];

				if (traceMaps[firstEdge].outResolved and traceMaps[lastEdge].inResolved) {


					ssize_t lastEdgeBal, firstEdgeBal;
					lastEdgeBal = balancedEdges[lastEdge];
					firstEdgeBal = balancedEdges[firstEdge];

					ssize_t balPathTraceIndex = traceMaps[lastEdgeBal].GetFirstStartTraceIndex();
					ssize_t balPathTrace = traceMaps[lastEdgeBal].traces[balPathTraceIndex].trace;

					//
					// Make sure the paths do not overlap
					//
					if (firstEdgeBal != firstEdge and 
							firstEdgeBal != lastEdge and 
							lastEdgeBal  != firstEdge and
							lastEdgeBal  != lastEdge) {
						ssize_t pi;
						/*
						cout << "detaching path: of length: " << pathTraces[p].edges->size() << ": ";
						for (pi = 0; pi < pathTraces[p].edges->size(); pi++) {
							cout << (*pathTraces[p].edges)[pi] << " ("
									 << edges[(*pathTraces[p].edges)[pi]].index << ") ";
						}
						cout << endl;

						cout << "detaching mate of length: " << pathTraces[balPathTrace].edges->size() << ": ";
						for (pi = 0; pi < pathTraces[balPathTrace].edges->size() ; pi++) {
							cout << (*pathTraces[balPathTrace].edges)[pi] << " (" 
									 << edges[(*pathTraces[balPathTrace].edges)[pi]].index << ") " ;
						}
						cout << endl;
						*/
						

						BalancedDetachPaths(pathTraces, traceMaps, graph, vertices, edges,
																paths, pathLengths,
																pathTraces[p], p,
																pathTraces[balPathTrace], balPathTrace,
																mateList, ignoreInconsistentPaths, 1);

						graphIsTransformed = 1;
						/*
						traceMaps[firstEdge].outResolved = 0;
						traceMaps[lastEdge].inResolved = 0;
						*/
						traceMaps[firstEdge].outResolved = traceMaps[lastEdge].outResolved;

						// The last edge should not be used since it is detached.
						traceMaps[lastEdge].inResolved   = 0;

						for (pi = 1; pi < pathTraces[p].edges->size() - 1; pi++) {
							internalDetachedEdges.insert((*pathTraces[p].edges)[pi]);
						}
			
						assert(traceMaps[lastEdgeBal].outResolved);
	
						for (pi = 1; pi < pathTraces[balPathTrace].edges->size() - 1; pi++) {
							internalDetachedEdges.insert((*pathTraces[balPathTrace].edges)[pi]);
						}
						/*
						traceMaps[lastEdgeBal].outResolved = 0;
						traceMaps[firstEdgeBal].inResolved = 0;
						*/
						traceMaps[lastEdgeBal].outResolved = traceMaps[firstEdgeBal].outResolved;
						traceMaps[firstEdgeBal].inResolved = 0;						

						//assert(graph.CheckAllPathsBalance(1));
						//assert(graph.CheckBalance());
					}
				}
			}
		
			// All of the internal edges should have been cleared.  Check out
			// how many weren't.
			set<ssize_t>::iterator setIt, setEnd;
			setEnd = internalDetachedEdges.end();
			for (setIt = internalDetachedEdges.begin(); setIt != setEnd; ++setIt) {
				//			cout << *setIt << " " << edges[*setIt].intervals->size() << " " << edges[*setIt].length << endl;
				ssize_t balEdge = balancedEdges[*setIt];
				if (internalDetachedEdges.find(balEdge) == internalDetachedEdges.end()) {
					cout << "edge: " << *setIt << " is detached, but not balance: " << balEdge << endl;
					exit(1);
				}
			}
			graph.RemoveErasedPaths();
			graph.RemoveMarkedIntervalsNoPaths();
			graph.RemoveUnlinkedEdges();
			//	assert(graph.CheckGraphStructureBalance());
			//graph.RemoveEmptyEdges();
			assert(graph.CheckGraphStructureBalance());
			ssize_t nremoved = graph.RemoveEmptyVertices();
			assert(nremoved >= 0);
			//			cout << "removed: " << nremoved << " vertices." << endl;
			graph.CondenseSimplePaths();
			assert(graph.CheckGraphStructureBalance());
			for (e = 0; e < edges.size(); e++ ){ 
				assert(edges[e].balancedEdge != -1);
			}
			// Free the allocated structures.
			for (p = 0; p < pathTraces.size(); p++) {
				if (pathTraces[p].edges != NULL) {
					pathTraces[p].edges->clear();
					delete pathTraces[p].edges;
				}
				pathTraces[p].edges = NULL;
			}
			pathTraces.clear();
			ssize_t m;
			for (m = 0; m < mateTraces.size(); m++ ){
				if (mateTraces[m].edges != NULL) {
					mateTraces[m].edges->clear();
					delete mateTraces[m].edges;
				}
				mateTraces[m].edges = NULL;
			}
			mateTraces.clear();

			traceMaps.clear();
			internalDetachedEdges.clear();
			balancedEdges.clear();
			pairedEdgeMap.clear();
			DeletePathTree(pathTree);
			++detachIter;
			//			assert(graph.CheckBalance());
			//			graph.CheckAllPathsBalance(1);


			if (!graphIsTransformed) {
				++currentPhase;
			}
			else {
				// Some extra cleaning.
				cout << "cleaning graph: " << curRule << " " << currentPhase << " " << detachIter << endl;
				//				graph.Erode(graph.vertexSize * 4);
				//				graph.RemoveLowCoverageEdges(5, -1);
			}
		}
		while (currentPhase < endPhase);
	}

	ssize_t numEmptyEdges = graph.RemoveEmptyEdges(); 
	cout << "Removed " << numEmptyEdges << " empty edges." << endl;

	graph.RemoveLowCoverageEdges(5, 5);
	graph.RemoveMarkedIntervalsNoPaths();
	if (runScaffold) {
		IntMatrix matchMat;
		InitScaffoldMatchMatrix(matchMat);
		for (curRule = 0; curRule < rules.size(); curRule++) {
			std::vector<ssize_t> vToRemove, eToRemove;
			MatePairScaffoldJoinEdges(graph, mateList, curRule, minMatePairCount, matchMat, vToRemove);
			graph.Prune(vToRemove, eToRemove);
			graph.CondenseSimplePaths();
		}
	}

	//UNUSED// ssize_t pi;

	string bGraphOutName = graphOutName + ".bgraph";
	string intvOutName = graphOutName + ".intv";
	string gvzOutName  = graphOutName + ".dot";
	string pathOutName = graphOutName + ".path";
	string edgeOutName = graphOutName + ".edge";
	string euGraphOutName = graphOutName + ".graph";
	CheckEdges(graph.vertices, graph.edges);
	graph.CondenseEdgeLists();
	graph.PrintIntervalGraph(bGraphOutName, intvOutName, report);


  PrintGraph(graph.vertices,graph.edges, euGraphOutName, report);
  PrintEdges(graph.vertices, graph.edges, edgeOutName, report);
  GVZPrintBGraph(graph.vertices, graph.edges, gvzOutName, report);
	WriteReadPaths(pathOutName, graph.paths, graph.pathLengths, report);


	EndReport(report);
	report.close();
	return 0;

}



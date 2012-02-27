/***************************************************************************
 * Title:          CleanGraphWithMates.cpp 
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
#include "PathLib.h"
#include <vector>
#include <iterator>
#include "IntegralTupleStatic.h"

using namespace std;


void PrintUsage() {
	std::cout << " usage: cleanGraphWithMates graphIn mateTable readDescription graphOut " << std::endl
						<< "  [-minPathCount c]  Remove paths with count less than 'c' " << std::endl
						<< "  [-minMatePairCount c] Remove paths with mate count less than 'c'" << std::endl;
}


ssize_t IsSetSourceContained(IntervalGraph& g, ssize_t source, std::set<ssize_t> &subset) {
	// Each vertex in 'subset' is either a dest of another vertex in subset, or
	// it is the source.

	std::set<ssize_t>::iterator subsetIt, subsetEnd;
	subsetEnd = subset.end();
	for (subsetIt = subset.begin(); subsetIt != subsetEnd; ++subsetIt) { 
		ssize_t inEdge, inEdgeIndex;
		if (*subsetIt == source) 
			continue;

		for (inEdgeIndex = g.vertices[*subsetIt].FirstIn();
				 inEdgeIndex != g.vertices[*subsetIt].EndIn();
				 inEdgeIndex = g.vertices[*subsetIt].NextIn(inEdgeIndex)) {
			inEdge = g.vertices[(*subsetIt)].in[inEdgeIndex];
			if (subset.find(g.edges[inEdge].src) == subset.end()){
				return 0;
			}
		}
		ssize_t outEdge, outEdgeIndex;
		for (outEdgeIndex = g.vertices[*subsetIt].FirstOut();
				 outEdgeIndex != g.vertices[*subsetIt].EndOut();
				 outEdgeIndex = g.vertices[*subsetIt].NextOut(outEdgeIndex)) {
			outEdge = g.vertices[(*subsetIt)].in[outEdgeIndex];
			if (subset.find(g.edges[outEdge].dest) == subset.end()){
				return 0;
			}
		}
	}
	return 1;
}


ssize_t FindSpanningVertexSet(IntervalGraph &g, ssize_t srcVertex, ssize_t destVertex, 
													ssize_t curVertex,
													ssize_t maxDepth,
													std::set<ssize_t> &vertices,
													std::set<ssize_t> &traversed) {
	// Didn't find the dest vertex, quit.
	if (maxDepth < 0) {
		return 0;
	}

	if (maxDepth >= 0 and curVertex == destVertex) {
		return 1;
	}
	ssize_t outEdge, outEdgeIndex;
	ssize_t spanFound = 0;
	ssize_t dest;
	for (outEdgeIndex = g.vertices[curVertex].FirstOut();
			 outEdgeIndex != g.vertices[curVertex].EndOut();
			 outEdgeIndex = g.vertices[curVertex].NextOut(outEdgeIndex)) {
		outEdge = g.vertices[curVertex].in[outEdgeIndex];
		dest = g.edges[outEdge].dest;
		if (traversed.find(dest) != traversed.end())
			continue;

		traversed.insert(dest);

		if (!FindSpanningVertexSet(g, srcVertex, destVertex, dest,
															 maxDepth - 1, vertices, traversed)) {
			// Not all the dest vertices from this vertex may reach the dest, quit.
			return 0;
		}
		spanFound = 1;
	}
	
	// There was no path to the dest
	if (spanFound == 0) {
		return 0;
	}

	// Now it is guaranteed that dest vertices from this vertex may reach destVertex.
	// Add them to the set vertices.
	for (outEdgeIndex = g.vertices[curVertex].FirstOut();
			 outEdgeIndex != g.vertices[curVertex].EndOut();
			 outEdgeIndex = g.vertices[curVertex].NextOut(outEdgeIndex)) {
		outEdge = g.vertices[curVertex].in[outEdgeIndex];
		vertices.insert(g.edges[outEdge].dest);
	}
	
	return 1;
}


int main(int argc, char* argv[]) {
	
	std::string graphFileName, mateTableName, graphOutName;
	std::string ruleFileName;
	if (argc < 5) {
		PrintUsage();
		exit(0);
	}
	graphFileName   = argv[1];
	mateTableName   = argv[2];
	ruleFileName    = argv[3];
	graphOutName    = argv[4];
	
	ssize_t minPathCount = 2;
	ssize_t minMatePairCount = 2;
	int argi = 5;
	ssize_t mateType = -1;
	while (argi < argc) {
		if (strcmp(argv[argi], "-minPathCount") == 0) {
			minPathCount = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-minMatePairCount") == 0) {
			minMatePairCount = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-mateType") == 0) {
			mateType = atoi(argv[++argi]);
		}
		else {
			PrintUsage();
			std::cout << "bad option: " << argv[argi] << std::endl;
			exit(0);
		}
		++argi;
	}

	IntervalGraph graph;
	ReadMateList  mateTable;
	RuleList rules;
	ParseRuleFile(ruleFileName, rules);

	std::cout << "Reading interval graph." << std::endl;
	int vertexSize;
	ReadIntervalGraph(graphFileName, graph, vertexSize);
	std::cout << "Reading mate table." << std::endl;
	ReadMateTable(mateTableName, mateTable);

	PathIntervalList &paths       = graph.paths;
	PathLengthList   &pathLengths = graph.pathLengths;
	TEdgeList        &edges       = graph.edges;
	TVertexList      &vertices    = graph.vertices;


	//UNUSED// ssize_t numRemoved = 0;
	/*
	numRemoved = RemoveBadStartEndMatePairs(graph, mateTable, rules);
	// Remove edges that have had the reads mapping to them removed.
	if (numRemoved > 0) {
		cout << "Removed " << numRemoved << " mate pairs." << endl;
		graph.RemoveMarkedIntervalsNoPaths();
		graph.RemoveEmptyEdges();
		graph.RemoveLowCoverageEdges(5,3);
		graph.CondenseSimplePaths();
	}
	*/

	// Now try and remove disjoint edges.
	ssize_t numCut = 0;
	ssize_t numSaved = 0;
	TVertexList newVertices;
	ssize_t appVertexIndex = vertices.size();
	ssize_t newVertexIndex = 0;
	TVertex tmpV;
	ssize_t e;
	ssize_t minPathSpan = 2;
	vector<ssize_t> cutIn, cutOut;
	for (e = 0; e < graph.edges.size(); e++) { 
		if (graph.vertices[edges[e].src].InDegree() > 0) {
			if (graph.IsEdgeInDisjoint(e, minPathSpan)) {

				// Now try and see if there is extra information that links 
				// e with edges before it.
				std::set<ssize_t> dsEdgeSet, traversedVertices;
				dsEdgeSet.insert(e);
				StoreDownstreamEdges(graph, edges[e].dest, 200 - 
														 (edges[e].length - vertices[edges[e].dest].vertexSize), 
														 dsEdgeSet);
				if (QueryUpstreamEdgesForDownstreamConnection(graph, mateTable,
																											edges[e].src, 200,
																											dsEdgeSet, 2, traversedVertices)) {
					std::cout << "Would have cut the edge, except mate pairs cross it." 
										<< std::endl;
					++numSaved;
				}
				else {
					// cut this edge from the source vertex
					cutIn.push_back(e);
				}
			}
		}
		if (graph.vertices[graph.edges[e].dest].OutDegree() > 0) {
			if (graph.IsEdgeOutDisjoint(e, minPathSpan)) {
				std::set<ssize_t> dsEdgeSet, traversedVertices;
				dsEdgeSet.insert(e);
				StoreDownstreamEdges(graph, edges[e].dest, 200, dsEdgeSet);
				if (QueryUpstreamEdgesForDownstreamConnection(graph, mateTable,
																											edges[e].src, 200,
																											dsEdgeSet, 2, traversedVertices)) {
					std::cout << "Would have cut the edge, except mate pairs cross it." 
										<< std::endl;
					++numSaved;
				}
				else {
				//				StoreDownstreamEdges(graph, edges[curEdge].dest, 200, dsEdgeSet);
					cutOut.push_back(e);
				}
			}
		}
	}

	ssize_t c;
	for (c = 0; c < cutIn.size(); c++) {
		e = cutIn[c];
		graph.vertices[graph.edges[e].src].EraseOutEdge(e);
		graph.edges[e].src = appVertexIndex;
		newVertices.push_back(tmpV);
		newVertices[newVertexIndex].AddOutEdge(e);
		cout << "cutting out edge " << e << " of length: " << graph.edges[e].length << endl;
		++appVertexIndex;
		++newVertexIndex;
		++numCut;
	}
	for (c = 0; c < cutOut.size(); c++ ){
		e = cutOut[c];
		graph.vertices[graph.edges[e].dest].EraseInEdge(e);
		graph.edges[e].dest = appVertexIndex;
		newVertices.push_back(tmpV);
		newVertices[newVertexIndex].AddInEdge(e);
		cout << "cutting in edge " << e << " of length " << graph.edges[e].length << endl;
		++appVertexIndex;
		++newVertexIndex;
		++numCut;
	}

	cout << "cut: " << numCut << " and saved: " << numSaved << endl;
	exit(0);
	/*
	 * Use read-paths to determine which edges are unique or not.
	 */

	PathBranch pathTree, removedPathTree;

	CollectPathTree(paths, pathLengths, pathTree);

	std::cout << "before trimming: " << std::endl;
	
	// Now check the graph for disjoint edges.
	std::set<ssize_t> origDisjointEdgeSet;
	for (e = 0; e < edges.size(); e++) {
		ssize_t src = edges[e].src;
		if (vertices[src].InDegree() != 0) {
			// Check src to see if it was made disjoint (no 
			// edges pass through the previous edge to here.
			if (graph.IsEdgeInDisjoint(e,2)) {
				std::cout << "EDGE " << e << " IS DISJOINT " << std::endl;
				origDisjointEdgeSet.insert(e);
			}
			
		}
	}

	TrimLowCoverageBranches(pathTree, minPathCount);
	RemoveLowCountPaths(graph, paths, pathLengths, pathTree, minPathCount);
	
	// Now check the graph for disjoint edges.
	std::set<ssize_t> lce;
	for (e = 0; e < edges.size(); e++) {
		ssize_t src = edges[e].src;
		if (vertices[src].InDegree() != 0) {
			// Check src to see if it was made disjoint (no 
			// edges pass through the previous edge to here.
			if (graph.IsEdgeInDisjoint(e, 2)) {
				std::cout << "EDGE " << e << " IS DISJOINT " << std::endl;
				lce.insert(e);
			}
		}
	}

std::set<ssize_t> newDisjointEdges1;
	std::insert_iterator<std::set<ssize_t> > setDiffInsert1(newDisjointEdges1, newDisjointEdges1.begin());
	std::set_difference(lce.begin(), lce.end(),
											origDisjointEdgeSet.begin(), origDisjointEdgeSet.end(),
											setDiffInsert1);

	std::cout << "The NEW disjoint edges are: " << std::endl;
	std::set<ssize_t>::iterator newDisIt;
	for (newDisIt = newDisjointEdges1.begin(); newDisIt != newDisjointEdges1.end(); ++newDisIt) {
		std::cout << *newDisIt << " ";
	}
	std::cout << std::endl;
	
	
	PathTraceList pathTraces;
	PathTreeToPathList(pathTree, pathTraces);

	MarkEnclosedPaths(pathTraces);
	RemoveEnclosedPaths(pathTraces);

	TraceMapMatrix traceMaps;
	traceMaps.resize(edges.size());
	StoreTraceMaps(pathTraces, traceMaps);
	MarkResolvedPaths(pathTraces, traceMaps);

	/*
	 * Use mate-pairs to try and find some paths that are bad.
	 */

	for (e = 0; e < edges.size(); e++ ){ 
		ssize_t intv;
		ssize_t numIntv = (*edges[e].intervals).size();
		for (intv = 0; intv < numIntv; intv++) {
			if ((*edges[e].intervals)[intv].markedForDeletion)
				continue;
			ssize_t path, pathPos;
			path = (*edges[e].intervals)[intv].read;
			pathPos = (*edges[e].intervals)[intv].pathPos;
			if (pathPos == 1) {
				// This path didn't start in this edge.  Maybe it's a bad path?
				// Do a test to see if it has a mate-pair that can't be reached from
				// the end of this path, but can be reached from the beginning.
				ssize_t firstEdge = graph.paths[path][0].edge;
				ssize_t firstEdgeIntv = graph.paths[path][0].index;
				ssize_t firstEdgePos = (*graph.edges[firstEdge].intervals)[firstEdgeIntv].edgePos;
				ssize_t firstEdgeReadPos = (*graph.edges[firstEdge].intervals)[firstEdgeIntv].readPos;

				ssize_t pathLength = graph.pathLengths[path];
				ssize_t lastEdge  = graph.paths[path][pathLength - 1].edge;
				ssize_t lastEdgeIntv = graph.paths[path][pathLength - 1].edge;
				ssize_t lastEdgePos =  (*graph.edges[lastEdge].intervals)[lastEdgeIntv].edgePos;
				//UNUSED// int lastEdgeReadPos = (*graph.edges[lastEdge].intervals)[lastEdgeIntv].readPos;

				if (firstEdge == lastEdge)
					// This path isn't branching
					continue;

				// Try and find the mate of these.
				if (path % 2 == 0) {
					ssize_t mateIndex = mateTable[path / 2].mateIndex;
					//UNUSED// ssize_t mateType  = mateTable[path / 2].mateType;

					// for now just use the 200 base inserts.
					if (mateIndex == -1)
						continue;

					ssize_t mateRead = mateIndex * 2 + 1;
					if (pathLengths[mateRead] == 0)
						continue;
					
					ssize_t mateEdge = paths[mateRead][0].edge;
					ssize_t mateEdgeIndex = paths[mateRead][0].index;
					ssize_t mateEdgePos = (*graph.edges[mateEdge].intervals)[mateEdgeIndex].edgePos;
					
					ssize_t readLength = graph.CalculateReadLength(path);

					//UNUSED// ssize_t endValid = 0;
					ssize_t lengthOffset = readLength - firstEdgeReadPos;
					std::list<ssize_t> lastPath, firstPath;
					ssize_t numLastPaths, numFirstPaths;
					ssize_t storePath = 1;
					MatePathList  lastMatePath, firstMatePath;
					ssize_t lastMatePathLength, firstMatePathLength;
					lastMatePathLength = firstMatePathLength = 0;
					map<ssize_t,ssize_t> visited;
					numLastPaths = CountValidMatePaths(graph, 
																						 lastEdge, lastEdgePos,
																						 mateEdge, mateEdgePos,
																						 0,
																						 200 - 100,
																						 200 + 100,
																						 10, 4, numLastPaths, 
																						 storePath, lastMatePath, 0, lastMatePathLength, visited);
					if (numLastPaths == 0) {
						visited.clear();
						numFirstPaths = CountValidMatePaths(graph, 
																								firstEdge, firstEdgePos,
																								mateEdge, mateEdgePos,
																								0,
																								200 - 100 + lengthOffset, 
																								200 + 100 + lengthOffset,
																								10, 4, numFirstPaths, 
																								storePath, firstMatePath, 0, firstMatePathLength, visited);
						if (numFirstPaths == 1) {
							std::cout << firstEdge << " " << lastEdge
												<< " (e: " << e << " fi " << edges[firstEdge].index 
												<< " li " << edges[lastEdge].index <<  ") num last paths was 0 and first: " 
												<< numFirstPaths << std::endl;
							// This path is likely bad. Remove it.
							graph.MarkPathForRemoval(path);
							ssize_t pathRC = path + 1;
							if (path % 2 == 1) {
								pathRC = path - 1;
							}
							
							graph.MarkPathForRemoval(pathRC);
							graph.pathLengths[pathRC] = 0;
							delete[] graph.paths[pathRC];
							graph.paths[pathRC] = NULL;
						}
					}
				}
			}
		}
	}
	graph.RemoveMarkedIntervalsNoPaths();
	// Now check the graph for disjoint edges.
	std::set<ssize_t> cleanedDisjointEdgeSet;
	for (e = 0; e < edges.size(); e++) {
		ssize_t src = edges[e].src;
		if (vertices[src].InDegree() != 0) {
			// Check src to see if it was made disjoint (no 
			// edges pass through the previous edge to here.
			if (graph.IsEdgeInDisjoint(e,2)) {
				std::cout << "edge: " << e << " of length: " << edges[e].length << " is disj with in " 
									<< vertices[src].InDegree() << std::endl;
				ssize_t in;
				std::cout << "in edge lengths: ";
				for (in = 0; in < 4; in++) {
					if (vertices[src].in[in] != -1) {
						std::cout << edges[vertices[src].in[in]].length << " ";
					}
				}
				std::cout << std::endl;
				cleanedDisjointEdgeSet.insert(e);
			}
		}
	}

	std::set<ssize_t> newDisjointEdges;
	std::insert_iterator<std::set<ssize_t> > setDiffInsert(newDisjointEdges, newDisjointEdges.begin());
	std::set_difference(cleanedDisjointEdgeSet.begin(), cleanedDisjointEdgeSet.end(),
											origDisjointEdgeSet.begin(), origDisjointEdgeSet.end(),
											setDiffInsert);

	std::cout << "The NEW disjoint edges are: " << std::endl;
	std::set<ssize_t>::iterator newDisIt1;

	ssize_t maxEdgeLength = 0;

	std::set<ssize_t> rcDisjointEdges;

	for (newDisIt1 = newDisjointEdges.begin(); newDisIt1 != newDisjointEdges.end(); ++newDisIt1) {
		std::cout << *newDisIt1 << " ";

		if (edges[*newDisIt1].length > maxEdgeLength) 
			maxEdgeLength = edges[*newDisIt1].length;

		ssize_t rcEdge = edges[*newDisIt1].balancedEdge;
		rcDisjointEdges.insert(rcEdge);
	}
	std::cout << std::endl;


	graph.RemoveMarkedIntervalsNoPaths();
	graph.RemoveEmptyEdges();
	graph.CondenseSimplePaths();


	//	exit(0);
											

	for (e = 0; e < edges.size(); e++ ) {
		// This edge is a branching edge.  That means that
		// it is possible that a branch is incorrect.  

		MateEdgeMap srcMateEdges;
		CollectMateEdges(graph, mateTable, e, srcMateEdges, mateType);
		MateEdgeMap::iterator mateIt;
		RemoveLowCountMateEdges(srcMateEdges, 3);

		MateEdgeMap::iterator mateEdgeIt, mateEdgeEnd;
		mateEdgeEnd = srcMateEdges.end();
		//UNUSED// ssize_t outEdge, outEdgeIndex;
		ssize_t dest = edges[e].dest;

		for (mateEdgeIt = srcMateEdges.begin(); mateEdgeIt != mateEdgeEnd; ++mateEdgeIt) {
			//UNUSED// ssize_t mateIsAdjacent = 0;
				
			std::set<ssize_t> containedVertices;
			std::set<ssize_t> traversedVertices;
			ssize_t mateEdgeSrc = edges[(*mateEdgeIt).first].src;
			if (mateEdgeSrc == dest) {
				// We can't find bulges between adjacent edges, so continue.
				continue;
			}
			// Consider a few rules: 
			//   This can't remove bulges in known repeat areas
			//   This can't remove bulges if there is a gap from the source edge.

			if (vertices[dest].InDegree() > 1)
				continue;
			if (vertices[dest].OutDegree() == 0)
				continue;
			if (vertices[mateEdgeSrc].OutDegree() > 1)
				continue;

			if (FindSpanningVertexSet(graph, dest, mateEdgeSrc, dest,
																4, containedVertices, traversedVertices)) {
				// Now look to see if the vertices in-between are not part of any
				// other paths.
				if (IsSetSourceContained(graph, dest, containedVertices)) {
					std::cout << "THE set: ";
					std::set<ssize_t>::iterator cit;
					for (cit = containedVertices.begin(); cit != containedVertices.end(); ++cit){ 
						std::cout << *cit << " ";
					}
					std::cout << " is wholly contained between two paired edges." << std::endl;
				}
			}
		}
	}
	graph.RemoveLowCoverageEdges(5,3);
	//	exit(0);

	std::string bGraphOutName = graphOutName + ".bgraph";
	std::string intvOutName = graphOutName + ".intv";
	std::string gvzOutName  = graphOutName + ".dot";
	std::string pathOutName = graphOutName + ".path";
	std::string edgeOutName = graphOutName + ".edge";
	std::string euGraphOutName = graphOutName + ".graph";
	CheckEdges(graph.vertices, graph.edges);

	graph.PrintIntervalGraph(bGraphOutName, intvOutName);


  PrintGraph(graph.vertices,graph.edges, euGraphOutName);
  PrintEdges(graph.vertices, graph.edges, edgeOutName);
  GVZPrintBGraph(graph.vertices, graph.edges, gvzOutName);
	WriteReadPaths(pathOutName, graph.paths, graph.pathLengths);
	
}

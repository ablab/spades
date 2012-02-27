/***************************************************************************
 * Title:          PathLib.h 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  11/01/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef PATH_LIB_H_
#define PATH_LIB_H_
#include <list>
#include "ReadPaths.h"
#include "PathBranch.h"
#include "mctypes.h"
#include "IntervalGraph.h"
#include "MateTable.h"
#include "Trace.h"


typedef std::pair<ssize_t, ssize_t> InOutEdgePair;

void CollectPathTree(PathIntervalList &paths,
										 PathLengthList &pathLengths,
										 PathBranch &pathTree);

void CollectPathTreeOnPath(PathInterval* path, ssize_t pathLength,
													 PathBranch &pathTree);

void PrintPathTree(PathBranch &pathTree);
void PrintPathTree(PathBranch &pathTree, std::list<ssize_t> &curPath);
void PrintEdgeTraces(PathTraceList &pathTraces, TraceMapMatrix &traceMaps, ssize_t edge);

ssize_t FindPathCoverage(PathInterval *path, ssize_t pathLength,
										 PathBranch &pathTree);

void DeletePathTraceList(PathTraceList &list);

ssize_t TrimLowCoverageBranches(PathBranch &pathTree, ssize_t minCount);

void RemoveLowCountPaths(IntervalGraph    &g,
												 PathIntervalList &paths,
												 PathLengthList   &pahtLengths,
												 PathBranch       &removedPathTree,
												 ssize_t minPathCount);

void RemoveLowCountTraces(PathTraceList &pathTraces, ssize_t minCount, PathBranch &removedPathTree);

void PathTreeToPathList(PathBranch &pathTree, PathTraceList &pathTraces);

void PathTreeToPathList(PathBranch &pathTree, 
												std::list<ssize_t> &curPath,
												PathTraceList &pathTraces);

void RemoveEnclosedPaths(PathTraceList &pathTraces);
ssize_t LocatePathStart(PathTraceList &pathTraces, ssize_t startEdge);

void MarkEnclosedPaths(PathTraceList &pathTraces);
void RemoveEnclosedPaths(PathTraceList &pathTraces);

void MarkExInternal(PathTraceList &pathTraces, TraceMapMatrix &traceMap);

void PrintPathTraces(PathTraceList &pathTraces, IntervalGraph &graph,
										 std::ostream &out);


ssize_t CheckForwardConsistency(PathTrace& path1, ssize_t pos1, PathTrace &path2, ssize_t pos2,
														ssize_t minMatchLength = 1);

ssize_t CheckReverseConsistency(PathTrace &path1, ssize_t pos1, PathTrace &path2, ssize_t pos2,
														ssize_t minMatchLength = 1);


// Store a map from an edge to the paths that go through it.

void StoreTraceMaps(PathTraceList &pathTraces, 
										TraceMapMatrix & traceMap  );

// Count how many edges are adjacent to this one in a postman traversal.
void CountAdjacencies(PathTraceList &pathTraces, TraceMapMatrix &traceMap);

// Checks to see if an edge is contained as the center
// of other paths.

ssize_t IsEdgePathContained(PathTraceList &pathTraces,
												TraceMapMatrix &traceMap,
												ssize_t edgeIndex);


// Looks to see if there is a 1-1 correspondence of entrance and exit edges
// for a repeat.
ssize_t AreEntranceAndExitEdgesPaired(PathTraceList &pathTraces,
																	TraceMapMatrix &traceMap,
																	TVertexList &vertices,
																	TEdgeList &edges,
																	ssize_t centerEdgeIndex, std::set<InOutEdgePair> &edgePairs);

void CollectPathSequence(TVertexList &vertices, 
												 TEdgeList &edges,
												 PathTrace &path, 
												 SimpleSequence &seq);


ssize_t  ArePathAndTraceReverseConsistent(PathInterval* path, ssize_t pathPos, ssize_t pathLength,
																			PathTrace &trace, ssize_t tracePos, ssize_t &pathEnd);
ssize_t  ArePathAndTraceForwardConsistent(PathInterval* path, ssize_t pathPos, ssize_t pathLength,
																			PathTrace &trace, ssize_t tracePos, ssize_t &pathBegin);

ssize_t FindLongestPathTraceReverseConsistency(PathInterval *path,
																					 ssize_t pathPos, ssize_t pathLength,
																					 PathTrace &trace, ssize_t tracePos,
																					 ssize_t &pathBegin);

ssize_t  FindLongestPathTraceForwardConsistency(PathInterval *path, 
																						ssize_t pathPos, ssize_t pathLength,
																						PathTrace &trace, ssize_t tracePos, 
																						ssize_t &pathEnd);

ssize_t FindLongestPathTraceOverlap(PathInterval *path, ssize_t pathLength, ssize_t &pathBegin, ssize_t &pathEnd,
																PathTrace &trace, 
																ssize_t &traceBegin, ssize_t &traceEnd, ssize_t findFirst);

ssize_t TraceContainsDuplications(PathTrace &trace );


ssize_t CountForwardConsistent(PathTrace &refPath, ssize_t refPos, ssize_t matchLength,
													 PathTraceList &pathTraces,
													 TraceMapMatrix &traceMap,
													 ssize_t edgeIndex, ssize_t &forwardConsistentTrace);

ssize_t CountReverseConsistent(PathTrace &refPath, ssize_t refPos, ssize_t matchLength,
													 PathTraceList &pathTraces,
													 TraceMapMatrix &traceMap,
													 ssize_t edgeIndex, ssize_t &reverseConsistentTrace);

ssize_t IsPathEndResolved(PathTrace &trace, TraceMapMatrix &traceMap);

ssize_t AreInternalPathEdgesResolved(PathTraceList &pathTraces,
																 TraceMapMatrix &traceMap,
																 ssize_t pathIndex);

ssize_t IsTangleEdgeResolved(PathTraceList &pathTraces,
												 TraceMapMatrix &traceMap,
												 ssize_t edgeIndex);

ssize_t MarkResolvedPaths(PathTraceList &pathTraces, TraceMapMatrix &traceMaps, ssize_t notStrict = 0 );

void PrintPathTraceList(PathTraceList &pathTraces);

ssize_t CheckPathTraceListBalance(TEdgeList &edges, PathTraceList &pathTraces);

ssize_t CheckPathBalance(TEdgeList &edges, PathTraceList &pathTraces, TraceMapMatrix &traceMaps);

ssize_t AreTracesForwardConsistent(PathTrace &traceA, ssize_t traceAPos, 
															 PathTrace &traceB, ssize_t traceBPos);
ssize_t AreTracesReverseConsistent(PathTrace &traceA, ssize_t traceAPos, 
															 PathTrace &traceB, ssize_t traceBPos);
ssize_t AreTracesConsistent(PathTrace &traceA, ssize_t traceAPos, 
												PathTrace &traceB, ssize_t traceBPos);

ssize_t FindCompatibleTrace(PathTraceList &pathTraces,
												TraceMapMatrix &traceMap,
												PathTrace &trace, ssize_t traceIndex, ssize_t tracePos,
												ssize_t edgeIndex,
												ssize_t &compatibleTrace, ssize_t &compatibleTraceIndex);

ssize_t ExtendPath(PathTraceList &pathTraces,
							 TraceMapMatrix &traceMap,
							 PathTrace &trace,
							 ssize_t traceIndex, ssize_t tracePos,
							 PathTrace &newTrace);

ssize_t PrintCandidatePaths(PathTraceList &pathTraces,
												TraceMapMatrix &traceMaps, TEdgeList &edges);



// Provide easy common interface to detach path, one with mate pairs,
// one without.

void DetachPath(PathTraceList &pathTraces,
								TraceMapMatrix &traceMap,
								IntervalGraph &graph,
								TVertexList &vertices, TEdgeList &edges,
								PathIntervalList &paths, PathLengthList &pathLengths,
								PathTrace &trace, ssize_t traceIndex, ssize_t removeInconsistentPaths, ssize_t doPrint);

void DetachPath(PathTraceList &pathTraces,
								TraceMapMatrix &traceMap,
								IntervalGraph &graph,
								TVertexList &vertices, TEdgeList &edges,
								PathIntervalList &paths, PathLengthList &pathLengths,
								PathTrace &trace, ssize_t traceIndex, 
								ReadMateList &mateList, ssize_t removeInconsistentPaths, ssize_t doPrint);


void DetachPath(PathTraceList &pathTraces,
								TraceMapMatrix &traceMap,
								IntervalGraph &graph,
								TVertexList &vertices, TEdgeList &edges,
								PathIntervalList &paths, PathLengthList &pathLengths,
								PathTrace &trace, ssize_t traceIndex, 
								ssize_t isMatePath, ReadMateList &mateList, ssize_t removeInconsistentPaths, ssize_t doPrint);

ssize_t PathContainsEdge(PathInterval *path, ssize_t pathLength, ssize_t edgeIndex, ssize_t &pathIndex);

void PrintTracesAsReads(TVertexList &vertices, TEdgeList &edges, 
												PathTraceList &traces, ssize_t endEdgeLength,
												std::string pathSeqName,
												std::ostream &report = std::cout);

void PrintPathTraceResolution(TEdgeList &edges,
															PathTraceList &pathTraces, 
															TraceMapMatrix &traceMaps);

ssize_t TraceContainsDuplications(PathTrace &trace );

void BalancedDetachPaths(PathTraceList &pathTraces,
												 TraceMapMatrix &traceMap,
												 IntervalGraph &graph, 
												 TVertexList &vertices, TEdgeList &edges,
												 PathIntervalList &paths, PathLengthList &pathLengths,
												 PathTrace &trace, ssize_t traceIndex,
												 PathTrace &balTrace, ssize_t balTraceIndex);


void BalancedDetachPaths(PathTraceList &pathTraces,
												 TraceMapMatrix &traceMap,
												 IntervalGraph &graph, 
												 TVertexList &vertices, TEdgeList &edges,
												 PathIntervalList &paths, PathLengthList &pathLengths,
												 PathTrace &trace, ssize_t traceIndex,
												 PathTrace &balancedTrace, ssize_t balancedTraceIndex,
												 ReadMateList &mateList, 
												 ssize_t removeInconsistentPathIntervals,
												 ssize_t isMatePath);


/*
class PathBranch {
 public:
	ssize_t edge;
	std::vector<ssize_t> branches;
	PathBranch(ssize_t e) {
		edge = e;
	}
};


typedef std::vector<PathBranch> PathTree;
*/
void AddBranch(ssize_t trunk, ssize_t branch);
#endif

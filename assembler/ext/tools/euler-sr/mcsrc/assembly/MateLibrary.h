/***************************************************************************
 * Title:          MateLibrary.h 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  11/26/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef MATE_LIBRARY_H_
#define MATE_LIBRARY_H_

#include "MateTable.h"
#include "RuleList.h"
#include "IntervalGraph.h"
#include "Trace.h"
#include "PathBranch.h"
#include "PathLib.h"
#include <list>
#include <set>
using namespace std;
class MatePathStatistics {
public:
	ssize_t numPaths;
	ssize_t meanPathLength;
	ssize_t pathCount;
	MatePathStatistics(ssize_t np, ssize_t mpl, ssize_t pc) {
		numPaths = np;
		meanPathLength = mpl;
		pathCount = pc;
	}
};

class EdgePair {
public:
	ssize_t edge1, edge2;
	ssize_t mateType;
	bool operator==(const EdgePair &p) const {
		return (p.edge1 == edge1 and p.edge2 == edge2);
	}
	bool operator<(const EdgePair &p) const {
		if (p.edge1 == edge1) {
			if (p.edge2 == edge2) {
				return mateType < p.mateType;
			}
			else {
				return edge2 < p.edge2;
			}
		}
		else {
			return edge1 < p.edge1;
		}
	}
	EdgePair &operator=(const EdgePair &p) {
		if (this != &p) {
			edge1 = p.edge1;
			edge2 = p.edge2;
			mateType = p.mateType;
		}
		return *this;
	}
};


class MatePairMap {
public:
	ssize_t count;
	ssize_t edge1, edge2;
	ssize_t meanEdge1End, meanEdge2Start;
	ssize_t numPaths;
	ssize_t averageNumEdges;
	ssize_t read1Length, read2Length;
};

typedef std::map<EdgePair, MatePairMap> EdgePairMap;

class MatePathInterval {
 public:
	ssize_t edge;
	MatePathInterval(ssize_t e) {
		edge = e;
	}
};

typedef std::list<MatePathInterval> MatePathList;

void ReadMateTable(std::string &MateTableName,
									 ReadMateList &mateList,
									 std::ostream &report = std::cout
									 );

ssize_t CountValidMatePaths(IntervalGraph &g,
												ssize_t curEdge, ssize_t curEdgePos,
												ssize_t endEdge, ssize_t endEdgePos,
												ssize_t curLength,
												ssize_t minLength, ssize_t maxLength,
												ssize_t maxDepth, ssize_t maxValid, ssize_t &numValid,
												ssize_t &storePath, MatePathList &matePath, 
												ssize_t print, ssize_t &totalPathLength,
												map<ssize_t,ssize_t> &visited);

void MarkEdgesOnValidMatePaths(IntervalGraph &g,
															 ssize_t curEdge, ssize_t curEdgePos,
															 ssize_t endEdge, ssize_t endEdgePos,
															 ssize_t curLength, ssize_t minLength, ssize_t maxLength,
															 ssize_t maxDepth,	std::vector<ssize_t> &path, ssize_t curDepth, 
															 std::vector<ssize_t> &edgeOnPathCount, 
															 ssize_t &numPaths, ssize_t &totalLength, ssize_t &totalEdges,
															 std::vector<ssize_t> &edgeTraversalCount);

ssize_t AreEdgesPaired(IntervalGraph &g,
									 ssize_t curEdge, ssize_t pairedEdge,
									 ReadMateList &mateList,
									 ssize_t curInterval);

void RemoveLowFrequencyEdgePairs(IntervalGraph &graph, 
																 EdgePairMap &matePairs,
																 ReadMateList &mateList, 
																 ssize_t minMatePairCount,
																 ssize_t ruleType=-1);

ssize_t StoreEdgePairMap(IntervalGraph &graph, 
										 ReadMateList &matePairs,
										 EdgePairMap &edgePairs,
										 ssize_t ruleType = -1);

void AssignMateOrder(ssize_t p, ssize_t mateIndex, ssize_t &mp, ssize_t&fisrtMate, ssize_t&secondMate);


void StoreMeanEdgePosition(EdgePairMap &matePairs);

void StoreUniqueMatePaths(IntervalGraph &graph, TVertexList &vertices, TEdgeList &edges,
													RuleList &rules,
													EdgePairMap &matePairs, PathBranch &pathTree);

class MatePairLocation {
 public:
	ssize_t count;     // count including read intervals
	ssize_t avgEndPos, avgStartPos;
	ssize_t cloneCount; // count of clones starting on corresponding edge
	MatePairLocation() {
		count = 0;
		avgEndPos = avgStartPos = 0;
		cloneCount = 0;
	}
	MatePairLocation(ssize_t c, ssize_t s, ssize_t a, ssize_t cc) {
		count = c;
		avgStartPos = s;
		avgEndPos = a;
		cloneCount = cc;
	}
};

typedef std::map<ssize_t, MatePairLocation> MateEdgeMap;

ssize_t StoreTreeDepth(IntervalGraph &g, ssize_t srcEdge, MateEdgeMap &mateEdges,
									 std::vector<ssize_t>  &distToEnd);

ssize_t SearchForMateEdge(IntervalGraph &g, ssize_t rootEdge, ssize_t srcEdge, 
											ssize_t maxSearchLength, ssize_t maxSearchDepth, 
											ssize_t mateEdge, std::list<ssize_t> &path);

ssize_t SearchForMateEdge(IntervalGraph &g, ssize_t rootEdge, ssize_t srcEdge, 
											ssize_t maxSearchLength, ssize_t maxSearchDepth, 
											MateEdgeMap &mateEdges, std::list<ssize_t> &path, ssize_t &altDestEdge);

void FindTreeDepth(IntervalGraph &g, ssize_t srcEdge, 
									 MateEdgeMap &mateEdges, std::set<ssize_t> &extraEdges,
									 std::vector<ssize_t> &distToSrc);

void CollectMateEdges(IntervalGraph &g, ReadMateList &readMates, ssize_t edge, 
											MateEdgeMap &mateEdges, ssize_t mateType, ssize_t strand=0);

ssize_t GetAverageMateStartPos(IntervalGraph &g, ReadMateList &readMates, 
													 ssize_t srcEdge, ssize_t destEdge, ssize_t &avgSrcEnd, ssize_t &avgDestEnd);

/*int MarkPathTree(IntervalGraph &g, int srcEdge,
								 MateEdgeMap &srcMateEdges, 
								 std::set<ssize_t> &additionalEdges );
*/

void ClearDistances(MateEdgeMap &readMates, std::set<ssize_t> &extraEdges,
										std::vector<ssize_t> distToSrc);

void UnmarkMateEdges(IntervalGraph &g, MateEdgeMap &readMates,
										 std::set<ssize_t> &extraEdges);

void UntraverseMateEdges(IntervalGraph &g, MateEdgeMap &readMates,
												 std::set<ssize_t> &extraEdges);

void RemoveLowCountMateEdges(MateEdgeMap &mateEdges, ssize_t minCount);

ssize_t ExtractMin(std::map<ssize_t, ssize_t> &distMap, std::set<ssize_t> &removed);

ssize_t FindClosestVertex(IntervalGraph &g, ssize_t startVertex, std::set<ssize_t> &destVertices);

ssize_t FindMatePath(IntervalGraph &g, ssize_t srcEdge,
								 MateEdgeMap &srcMateEdges,
								 std::list<ssize_t> &srcEdgePath);

void FindScaffoldPaths(IntervalGraph &g, 
											 ReadMateList &mateTable,
											 RuleList &mateRules, ssize_t mateType, ssize_t minScaffoldLength,
											 PathIntervalList &paths,
											 PathLengthList &pathLengths);

ssize_t FindPairedOutEdge(IntervalGraph &g, MateEdgeMap &pairedEdges, ssize_t curVertex, ssize_t &pairedOutEdge);


ssize_t CollectPairedOutEdges(IntervalGraph &g, MateEdgeMap &pairedEdges, ssize_t curVertex, 
													set<ssize_t> &pairedOutEdges);


ssize_t AdvancePathAlongPairedEdges(IntervalGraph &g, ssize_t curEdge, MateEdgeMap &pairedEdges, 
																ssize_t stopEdge, 
																std::list<ssize_t> &pairedPath, ssize_t &pathSeqLength,
																map<ssize_t,ssize_t> &edgeTarversalCount,
																ssize_t &lastEdge, ssize_t &numPairedOutEdges);

ssize_t FindPairedPath(IntervalGraph &g, MateEdgeMap &pairedEdges, 
									 ssize_t curEdge, ssize_t curEdgePos, ssize_t destEdge, ssize_t destEdgePos,
									 ssize_t curLength, ssize_t minLength, ssize_t maxLength,
									 std::list<ssize_t> &curPath, ssize_t &numCurPathPairedEdges,
									 std::list<ssize_t> &optimalPath, ssize_t &numOptPathPairedEdges, ssize_t &numOptPaths);

void GetLastStartPosition(MateEdgeMap &mateEdgeMap, ssize_t &startPos, ssize_t &pairedEdge);


ssize_t IsLengthValid(ssize_t length, ssize_t min, ssize_t max);

void MarkScaffoldEdges(IntervalGraph &g, ssize_t minCoverage, ssize_t minLength);

ssize_t FindPairedScaffoldEdges(IntervalGraph &g, ReadMateList &readMates, ssize_t edge, 
														std::set<ssize_t> &pairedScaffoldEdges);

ssize_t FindMaximallySupportedPath(IntervalGraph &g, ReadMateList &mateTable, 
															 RuleList &mateRules, ssize_t mateType,
															 ssize_t srcEdge, ssize_t destEdge,
															 std::list<ssize_t> &optimalPath) ;

ssize_t FindMaximallySupportedPath(IntervalGraph &g, MateEdgeMap &pairedEdges, 
															 ssize_t curEdge, ssize_t curEdgePos, ssize_t destEdge, ssize_t destEdgePos,
															 ssize_t curLength, ssize_t minPathLength, ssize_t maxPathLength,
															 ssize_t maxSearchDepth,
															 std::list<ssize_t> &curSupportedPath, ssize_t curSupportedPathScore,
															 std::list<ssize_t> &maxSupportedPath, ssize_t &maxSupportedPathScore,
															 map<ssize_t,ssize_t> &edgeTraversalCount,
															 ssize_t &numOptPaths);

ssize_t QueryUpstreamEdgesForDownstreamConnection(IntervalGraph &g, ReadMateList &mateTable,
																							ssize_t curVertex, ssize_t radius, 
																							std::set<ssize_t> &dsEdgeSet, ssize_t minPairedEdges,
																							std::set<ssize_t> &traversedVertices);


void StoreDownstreamEdges(IntervalGraph &g, ssize_t curVertex, ssize_t radius, 
													std::set<ssize_t> &dsEdgeSet);

ssize_t RemoveBadStartEndMatePairs(IntervalGraph &g, ReadMateList &mateTable, RuleList &rules, ssize_t ruleType);

ssize_t ComputeMatePairLengthDistribution(IntervalGraph g, ReadMateList &readMates,
																			ssize_t mateType, ssize_t &meanSep, double &stddevSep);


ssize_t  PathsMayOverlap(Path &forPath, ssize_t forPathLength,
										 Path &revPath, ssize_t revPathLength,
										 ssize_t minCloneSize, ssize_t maxCloneSize );

ssize_t  FindOverlappingMatePaths(Path &forPath, ssize_t forPathLength,
															Path &revPath, ssize_t revPathLength,
															Path &overlappingPath, ssize_t &overlappingPathLength);

#endif

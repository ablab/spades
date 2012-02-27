/***************************************************************************
 * Title:          IntervalGraph.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/19/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef BRANCHING_GRAPH_H_
#define BRANCHING_GRAPH_H_

#include "BVertex.h"
#include "BEdge.h"
#include "DeBruijnGraph.h"
#include "ReadIntervals.h"
#include "ReadPaths.h"
#include "SimpleSequence.h"
#include "graph/GraphAlgo.h"
#include "graph/MSTAlgo.h"
#include "ThreadPath.h"
#include "BinomialHeap.h"
#include "AlternativeEdge.h"
#include "compatibility.h"
#include <fstream>
#include <map>
#include <list>


using namespace std;

// TODO: Fix sparc compiler warning
// "warning: `class TVertex' has virtual functions but non-virtual destructor"
class TVertex : public BVertex {
 public:
	// HAS BITFIELD
  bool toDelete : 1;
  bool traversed: 1;

  // Ok, this is probably not the best coding style, but I
  // need another code for flagging, and I can probably have up to 32 flags
  // without changing the space requirements.

  TVertex() : BVertex() { 
    marked  = GraphVertex::NotMarked;
    flagged = GraphVertex::NotMarked;
  }
  char IsMarked() { return marked;}
  void Mark() {marked = GraphVertex::Marked;}
  void Unmark() {
    marked  = GraphVertex::NotMarked;
    flagged = GraphVertex::NotMarked;
    traversed = GraphVertex::NotMarked;
  }
	int vertexSize;
	TVertex& operator=(const TVertex &rhs) {
		if (this != &rhs) {
			*((BVertex*)this) = (BVertex) rhs;
			toDelete = rhs.toDelete;
			traversed = rhs.traversed;
			vertexSize = rhs.vertexSize;
		}
		return *this;
	}
};

class TEdge : public IntervalEdge {
 public:

  TEdge() : IntervalEdge() {
    marked       = GraphEdge::NotMarked;
    mst          = GraphAlgo::MSTOut;
    flagged      = GraphEdge::NotMarked;
    balPreferred = GraphEdge::NotMarked;
		guarded      = GraphEdge::NotMarked;
    // This is used to record if an edge is likely to be erroneous.
    suspect      = GraphEdge::NotMarked;
    // Sometimes it is desirable to mark an edge as preferred
    // so that if there are operations that should be done in order
    // to preserve the balancedness of a graph, the proper order may be maintained.
    // Normally this isn't a problem, but if there is a tie for preference in 
    // processing order, then use balPreferred to break the tie.
    // What that means that when the operation is first performed to the graph,
    // the balPreferred should be set (so that later when the balaned edge is processed
    // preference will be set.
    balPreferred = GraphEdge::NotMarked;
  }
	ssize_t IsShort(ssize_t minLength) {
		return (guarded == GraphEdge::NotMarked and
						length < minLength);
	}

	ssize_t IsLowCoverage(ssize_t coverage) {
		return (guarded == GraphEdge::NotMarked and
						(ssize_t) intervals->size() < coverage);
	}
	ssize_t CountReadsStartingInEdge() {
		size_t i;
		ssize_t numStarts = 0;
		for (i = 0; i < intervals->size(); i++) {
			if ((*intervals)[i].readPos == 0)
				numStarts++;
		}
		return numStarts;
	}
	
	ssize_t CountReadsCoveringEdge() {
		size_t i;
		ssize_t numCover = 0;
		for (i = 0; i < intervals->size(); i++) {
			if ((*intervals)[i].length == length) 
				numCover++;
		}
		return numCover;
	}

	void MarkIntervalForRemoval(ssize_t intervalIndex) {
		assert(intervalIndex < (ssize_t) intervals->size());
		if ((*intervals)[intervalIndex].markedForDeletion == 0) {
			multiplicity--;
			(*intervals)[intervalIndex].markedForDeletion = 1;
		}
	}

	ssize_t IsUnexpectedlyLowCoverage(double readsPerNucleotide, 
																double maxStdDev, 
																ssize_t absoluteCutoff,
																ssize_t acceptedLength) {
		if ((ssize_t) intervals->size() < absoluteCutoff)
			return 1;

		if (length >= acceptedLength) 
			return 0;
		
		if (intervals->size() > 6)
					return 0;
		double expIntervals;
		expIntervals = readsPerNucleotide * length;
		// Assume intervals are distributed as the sum of a Poisson 
		double stddevIntervals = sqrt(expIntervals * readsPerNucleotide);

		double normStdDev;
		ssize_t numCoveringReads; 
    numCoveringReads = intervals->size();

		normStdDev = (expIntervals - numCoveringReads) / stddevIntervals;
		if (normStdDev > maxStdDev) {
			return 1;
		}
		return 0;
	}

  char IsMarked() {return marked;}
  void Mark() {marked = GraphEdge::Marked;}
  void Unmark() {
    marked       = GraphEdge::NotMarked;
    flagged      = GraphEdge::NotMarked;
    traversed    = GraphEdge::NotMarked;
    balPreferred = GraphEdge::NotMarked;
		guarded      = GraphEdge::NotMarked;
  }
  // This plus one more field inherited from GraphEdge
  // Will need to upgrade to an integer field soon
	// HAS BITFIELD
  bool toDelete     : 1;
  bool traversed    : 1;
  bool balPreferred : 1;
  bool suspect      : 1;
	bool guarded      : 1;
  TEdge &operator=(const TEdge &edge) {
		if (this != &edge) {
			IntervalEdge::operator=(edge);
			this->marked       = edge.marked;
			this->flagged      = edge.flagged;
			this->traversed    = edge.traversed;
			this->balPreferred = edge.balPreferred;
			this->suspect      = edge.suspect;
			this->guarded      = edge.guarded;
		}
		return *this;
  }
	void Clear() {
		((IntervalEdge*)this)->Clear();
	}
	AlternativeEdgeVect altEdges;
};

typedef std::vector<TVertex> TVertexList;
typedef std::vector<TEdge> TEdgeList;



class IntervalGraph {
 public:
  TVertexList vertices;
  TEdgeList edges;
  _INT_ vertexSize;
	_INT_ containsIntervals;
	_INT_ storeAlternativeEdges;
  std::vector<SimpleSequence> *edgeSeqPtr;
	std::vector<AlternativeEdge> altEdges;
  IntervalGraph() {
    edgeSeqPtr = NULL;
    vertexSize = -1;
    isBalanced = 1;
		containsIntervals = 1;
		storeAlternativeEdges = 0;
  }
	ssize_t GetBalancedPathIndex(ssize_t p) {
		if (p % 2 == 0) 
			return p+1;
		else 
			return p - 1;
	}
	void ReadAlternativeEdges(const char *altEdgeInName, std::ostream &report = std::cout);
	//	void ReadAlternativeEdges(string &altEdgeInName);
	void WriteAlternativeEdges(string &altEdgeOutName, std::ostream &report = std::cout);
  void CalcDMST(); 
	void PathToThread(ssize_t p, ThreadPath &path);
  ssize_t FindAlternatePaths();
	void PrintAlternativeEdges(string &altEdgeOutName);
	
  void RemoveAllPathsThroughEdge(ssize_t edge);

  void RemovePath(std::vector<ssize_t> &pathEdges,
									std::vector<ssize_t> &pathIndices);

  void TracePath(ssize_t curEdge, ssize_t curIndex, 
								 std::vector<ssize_t> &pathEdges,
								 std::vector<ssize_t> &pathIndices);

  void MergeSimplePath(ssize_t vertex, ssize_t toEdge, ssize_t fromEdge);

  void ReadIntervalGraph(std::string &bGraphFileName, 
												 std::string &intervalFileName,
												 std::string pathName = "",
												 _INT_ skipIntervals = 0,
												 std::ostream &report = std::cout
												 );

  void ReadIntervalGraph(std::string &bGraphFileName, 
												 std::string &intervalFileName,
												 std::string pathName,
												 std::ostream &report
												 );
  
  void PrintIntervalGraph(std::string &bGraphName,
													std::string &intervalName,
													std::ostream &report = std::cout
													);

  void RemoveAllSimpleBulges(ssize_t minBulgeSize);
  ssize_t RemoveSmallComponents(ssize_t size, 
														std::string &readsFile,
														std::string &componentReadsFile,
														std::ostream &report = std::cout
														);

  void Erode(ssize_t minEdgeLength);
  void RemoveAllButMST();
	ssize_t CondenseSimplePaths();

  ssize_t TracePathReverse(ssize_t curEdge, ssize_t curIndex, ssize_t &prevEdge, ssize_t& prevIndex);
  ssize_t TracePathForwards(ssize_t curEdge, ssize_t curIndex, ssize_t &nextEdge, ssize_t& nextIndex);
	ssize_t RemoveEmptyEdges();
	ssize_t RemoveEmptyVertices();
  ssize_t RemoveLowCoverageEdges(double lowCoverageStddev, ssize_t absoluteCutoff);
  ssize_t FindShortContainedDAG(ssize_t sourceVertex, ssize_t maxPathLength, 
														std::set<ssize_t> &dagVertices,
														std::set<ssize_t> &dagEdges,
														std::vector<ssize_t> &shortestPath);

  void Unmark();
  void RemoveBalPreference();
  void Unflag();
  void Untraverse();
  void RemoveBalPreferred();
  void Unsuspect();
  void InitializeFlags();
  void IncrementOptimalPathCount(ssize_t source, 
																 std::vector<ssize_t> &sinks, 
																 std::vector<ssize_t> &edgeTraversalCount);

  ssize_t FindHighestMultiplicityPath(ssize_t begin, ssize_t end, std::vector<ssize_t> &path);
	void CalcEdgeMultiplicityStats(double &avgMultip);
  void MarkSuspectEdges(ssize_t polyNucLength);
  ssize_t RemoveSuspectBulges();
  void FindHighestScoringPathSourceToSink(ssize_t source, ssize_t sink, std::vector<ssize_t> &edgeTraversalCount);
  void TraceSourceToSinkPaths(std::vector<ssize_t> &edgeTraversalCount);
  void FindSourcesAndSinks(std::vector<ssize_t> &sources,
													 std::vector<ssize_t> &sinks);
	void RemoveBulges(ssize_t bulgeLength, ssize_t useDMST);
	void RemoveWhirls(ssize_t whirlLength);
	ssize_t StorePathForwards(ssize_t curEdge, ssize_t curIndex,
												std::vector<ssize_t> &pathEdges,
												std::vector<ssize_t> &pathIndices);

	ssize_t StorePathReverse(ssize_t curEdge, ssize_t curIndex,
											 std::vector<ssize_t> &pathEdges,
											 std::vector<ssize_t> &pathIndices);

	void DiscardGappedPaths();
	void RemoveMarkedPathIntervals();
	void RemoveTruncatedPathIntervals();
	void GrowIntervalsOnSimplePath(ssize_t edgeIndex);
	ssize_t CollectSimplePathIntervals(ssize_t startEdge, ssize_t numNewIntervals);
	void MarkIntervalForRemoval(ssize_t edgeIndex, ssize_t intervalIndex) {
		assert(edgeIndex >= 0);
		assert(intervalIndex >=0);
		if (IsIntervalMarkedForRemoval(edgeIndex, intervalIndex))
			return;
		ssize_t readIndex =	(*edges[edgeIndex].intervals)[intervalIndex].read;
		ssize_t pathPos   =	(*edges[edgeIndex].intervals)[intervalIndex].pathPos;
		// Clear the path
		assert(pathPos < pathLengths[readIndex]);
		paths[readIndex][pathPos].edge = -1;
		paths[readIndex][pathPos].index = -1;
		
		// Clear the edge interval
		(*edges[edgeIndex].intervals)[intervalIndex].markedForDeletion = 1;
		edges[edgeIndex].multiplicity--;
	}
	PathLengthList pathLengths;
	PathIntervalList paths;

	size_t IsIntervalMarkedForRemoval(ssize_t edgeIndex, ssize_t intervalIndex) {
		if (edgeIndex >= 0 and intervalIndex >= 0)
			return (*edges[edgeIndex].intervals)[intervalIndex].markedForDeletion;
		else
			return 1;
	}
	ssize_t RemoveMarkedIntervals();
	ssize_t RemoveMarkedIntervalsNoPaths();
	void SetMaxReadIndex(ssize_t mri) {maxReadIndex=  mri;}
	ssize_t IsForwardRead(ssize_t read) {
		return (read % 2 == 0);
	}
	ssize_t GetCompReadIndex(ssize_t read) {
		if (read % 2 == 0) 
			return read + 1;
		else 
			return read - 1;
	}

	ssize_t RemovingEdgeCutsGraph(ssize_t e) {
		// The edge cuts the 
		/*
	std::cout << vertices[edges[e].src].OutDegree() << " "
						<< vertices[edges[e].src].InDegree() << " "
						<< vertices[edges[e].dest].OutDegree() << " "
						<< vertices[edges[e].dest].InDegree() << std::endl;
		*/
		return ((vertices[edges[e].src].OutDegree() == 1 and 
						 vertices[edges[e].src].InDegree() != 0) or
						(vertices[edges[e].dest].InDegree() == 1 and
						 vertices[edges[e].dest].OutDegree() != 0));
	}
	ssize_t CheckPathContinuity();
	ssize_t CheckAllPathsBalance(ssize_t fatal=0);
	ssize_t CheckAllPathsContinuity(ssize_t fatal=0);
	ssize_t CheckBalance();
	ssize_t CheckGraphStructureBalance();				
	ssize_t CheckPathContinuity(ssize_t p);
	ssize_t CheckPathBalance(ssize_t p, ssize_t balp);

	void ProtectEdges(std::string &protectedEdgeFileName,
										_INT_ edgeTupleLength,
										std::ostream &report = std::cout
										);

  ssize_t SearchForDirectedCycle(ssize_t sourceVertex, ssize_t curVertex, std::set<ssize_t> &visitedEdges,
														 ssize_t curPathLength, ssize_t maxCycleLength);

  ssize_t SearchForUndirectedCycle2(ssize_t sourceVertex, ssize_t curVertex, std::set<ssize_t> &curPath, 
																ssize_t curPathLength, ssize_t maxCycleLength, ssize_t isUndirected,
																std::set<ssize_t> &visitedEdges, std::set<ssize_t> &visitedVertices);
	

  void SkipOutEdge(ssize_t vertex, ssize_t inEdge, ssize_t outEdge);
  void Prune(std::vector<ssize_t> &verticesToRemove,
						 std::vector<ssize_t> &edgesToRemove);

	ssize_t CountReadsContainedInEdge(ssize_t edge);
	ssize_t CountReadsPassingThroughEdge(ssize_t edge);
	ssize_t CountReadsExtendingIntoEdge(ssize_t edge, ssize_t limit);

	void RemoveLowPathEdges(ssize_t minPaths, ssize_t minExtend);

	ssize_t ReplacePathRangeForward(ssize_t readIndex,
															ssize_t intvEdge, ssize_t intvIndex,
															ssize_t pathStart, ssize_t pathEnd,
															std::vector<ssize_t> &newEdges);

	ssize_t ReplacePathRangeReverse(ssize_t readIndex,
															ssize_t intvEdge, ssize_t intvIndex,
															ssize_t pathStart, ssize_t pathEnd,
															std::vector<ssize_t> &newEdges);

	ssize_t ReplacePathRange(ssize_t readIndex, 
											 ssize_t intvEdge, ssize_t intvIndex,
											 ssize_t pathStart, ssize_t pathEnd,
											 std::vector<ssize_t> &newEdges);
	void PrintPath(ssize_t p, std::ostream &pathOut);
	void PrintPathReverse(ssize_t p, std::ostream &pathOut);
	void UpdateAllPathIndices();
	void UpdatePathIndices(ssize_t edge);
	ssize_t CalculateReadLength(ssize_t read);

	void RemoveEdgeAndMarkIntervalsForRemoval(ssize_t edge,
																						std::vector<ssize_t> &removedVertices);

	void MarkIntervalsOnPathForRemoval(ssize_t path);
  void RemoveEdge(ssize_t edge, std::vector<ssize_t> &removedVertices);
	void MarkPathForRemoval(ssize_t path);
	ssize_t IsEdgeInDisjoint(ssize_t edge, ssize_t minSpanningPaths);
	ssize_t IsEdgeOutDisjoint(ssize_t edge, ssize_t minSpanningPaths);
	ssize_t CutDisjointEdges(ssize_t minSpanningPaths);
  void DisconnectEdgesAtSource(std::vector<ssize_t> &edgeList);
	void DisconnectEdgesAtDest(std::vector<ssize_t> &edgeList);
	ssize_t isBalanced;
	ssize_t PathToSequence(vector<ssize_t> &path, SimpleSequence &seq);
	ssize_t ComputePathLength(vector<ssize_t> &path);
	ssize_t GetPathLength(ssize_t path);
	void ComputeIntervalPositions(list<ssize_t> &path, vector<ssize_t> &positions);
	ssize_t AlignmentBoundariesToPath(ssize_t *alignment, ssize_t alignmentLength,
																ssize_t start, ssize_t end,
																vector<ssize_t> &path, vector<ssize_t> &edgeStarts, 
																vector<ssize_t> &pathEdges, 
																vector<ssize_t> &pathEdgeStarts, vector<ssize_t> &pathEdgeLengths);
	ssize_t ListToVector(list<ssize_t> &l, vector<ssize_t> &v);
	ssize_t LookupAlignedPosition(ssize_t edge, ssize_t edgePos, ssize_t *alignment, ssize_t seqLength, 
														vector<ssize_t> &path,
														vector<ssize_t> &edgeStarts, ssize_t traversal);

	ssize_t LocateEdgeOnPath(ssize_t edge, list<ssize_t> &path);


	void MergeOutEdges(ssize_t vertex,
										 ssize_t toEdge, ssize_t fromEdge);

	void AppendIntervals(ssize_t toEdge, ssize_t fromEdge);
	ssize_t CountIntervalsOnSimplePath(ssize_t edge);
	void MoveIntervals(ssize_t toEdge, ssize_t fromEdge, ssize_t toEdgeIntvStartIndex, ssize_t lengthOffset);
	void MoveIntervals(ssize_t toEdge, ssize_t fromEdge);

	void MergeInEdges(ssize_t vertex,
										ssize_t toEdge, ssize_t fromEdge);

	void RemovePath(ssize_t path);
	void MarkPathIntervalsForRemoval(ssize_t path);
	void RemoveUnlinkedEdges();
	void ErodeShortEndIntervals(ssize_t minIntvLength);
	void RemoveErasedPaths();
	void CondenseEdgeLists();
	void SortAllEdgeIntervalsByReadPos();
	void Free();
	ssize_t TraceBackEdges(ssize_t startVertex, ssize_t endVertex, ssize_t backEdges[], 
										 std::vector<ssize_t> &path);
	// The maximum read index among all read intervals
	ssize_t ExtractMin(set<ssize_t> &keyQueue, map<ssize_t,ssize_t> &values, ssize_t &minKey, ssize_t &minValue);
	ssize_t maxReadIndex;
	ssize_t VertexMapToEdgePath(map<ssize_t,ssize_t> &vertexPaths, ssize_t startVertex,
													list<ssize_t> &edgePath, ssize_t dir);

	void AssignVertexSizes();

	ssize_t SearchTwoDirectedCycle(ssize_t redVertex, ssize_t blackVertex, ssize_t maxLength,
														 ssize_t maxDegree, ssize_t &hasAboveMax,
														 BinomialHeapNode<ssize_t,ssize_t>* redNodeList[],
														 BinomialHeapNode<ssize_t,ssize_t>* blackNodeList[],
														 ssize_t redDistList[],
														 ssize_t blackDistList[],
														 ssize_t redInv[], ssize_t blackInv[], ssize_t invocation,
														 ssize_t &minVertex, ssize_t &minCycleLength,
														 ssize_t redPath[], ssize_t blackPath[]);

	void AddOutEdgeLengths(ssize_t curVertex, ssize_t lengthToCurVertex, ssize_t maxLength,
												 BinomialHeap<ssize_t,ssize_t> &lengthPQueue,
												 BinomialHeapNode<ssize_t,ssize_t> *nodeRef[],
												 ssize_t invocations[], ssize_t inv,
												 ssize_t useFlagged, ssize_t path[]);

	ssize_t MarkPathRangeForRemoval(ssize_t pathIndex, ssize_t start, ssize_t end);

	ssize_t FindMinimalMatchedArrivalLength(map<ssize_t,ssize_t> &setA, map<ssize_t,ssize_t> &setB, 
																			ssize_t &minElement, ssize_t&minLength);
	ssize_t  DoesPathRepeat(ssize_t readIndex, ssize_t pathPos);
	void ProtectEdges(ReadPositions &protectedPositions, 
										SimpleSequenceList &protectedEdges, 
										_INT_ edgeTupleLength);

	void FormComplimentPath(std::vector<ssize_t> &origPath,
													std::vector<ssize_t> &compPath);

						
	void PrintPaths(std::string pathFile, std::ostream &report = std::cout);
	ssize_t NumForwardReads() {
		return ((maxReadIndex+1)/2);
	}
	void SetMultiplicities();

	void RemoveEdgeList(std::vector<ssize_t> &edgesToRemove);
	void RemoveVertexList(std::vector<ssize_t> &verticesToRemove);

	void MarkPathForRemoval(std::vector<ssize_t> &pathEdges,
													std::vector<ssize_t> &pathIndices);

	void RemoveEdgeAndMarkPathsForRemoval(ssize_t edge,
																				std::vector<ssize_t> &removedVertices);

	void RemoveEdgeAndMarkIntervalsForRemoval(std::vector<ssize_t> &edgeList,
																						std::vector<ssize_t> &removedVertices);


	ssize_t MakeRoomForEdges(ssize_t readIndex, ssize_t intvEdge, ssize_t intvIndex,
											 ssize_t pathStart, ssize_t pathEnd,
											 std::vector<ssize_t> &newEdges);

	void CopyNewEdges(ssize_t readIndex, ssize_t intvEdge, ssize_t intvIndex,
										ssize_t pathStart, std::vector<ssize_t> &newEdges);
	void PrintImbalancedPaths(ssize_t p, ssize_t balp);
	void MarkPathsThroughEdgeForRemoval(ssize_t edgeIndex);
	void MarkIntervalsInEdgeForRemoval(std::vector<ssize_t> &edgeIndices);
	void MarkIntervalsInEdgeForRemoval(ssize_t edgeIndex);
	ssize_t RouteRemovedIntervals(ssize_t maxSearchLength);
	void SplicePathRange(ssize_t readIndex,
											 ssize_t spliceStart, ssize_t spliceEnd);
	ssize_t SearchAlternatePath(ssize_t curPathEdge, ssize_t nextPathEdge, 
													std::list<ssize_t> &altPathEdges,
													ssize_t maxSearchLength,
													ssize_t curPathLength = 0);

	ssize_t StoreAlternatePath(ssize_t curPathEdge, ssize_t nextPathEdge, 
												 std::vector<ssize_t> &altPathEdges,
												 ssize_t maxSearchLength);
	void DeleteReadInterval(ssize_t readIndex, ssize_t pathPos);
	void UntraverseReadIntervals();
	void AssignPathOrderToEdges();

	void AssignIntervalPathOrder();
	void AssignEdgesToIntervals();
	void SortAllEdgeIntervalsByEdgePos();
  void RemoveEdges(std::vector<ssize_t> &edgesToRemove, std::vector<ssize_t> &orphanedVertices);
  void MarkBalancedEdges();
	ssize_t RemoveBulgingEdges(ssize_t bulgeLength, ssize_t useDMST, ssize_t iter, ssize_t maxDegree, ssize_t &hasAboveMax);
	ssize_t SearchForUndirectedCycle(ssize_t sourceVertex, ssize_t cycleEndVertex, ssize_t maxCycleLength);
	ssize_t LookupBalVertex(ssize_t vertex);

  template<typename Container>
		void FlagEdgeSet(Container &edgeSet) {
		std::set<ssize_t>::iterator edgeIt;
    for (edgeIt = edgeSet.begin(); edgeIt != edgeSet.end(); ++edgeIt) {
      edges[*edgeIt].flagged = GraphEdge::Marked;
    }
  }


  template<typename Container>
		void TraverseEdgeSet(Container &edgeSet) {
		std::set<ssize_t>::iterator edgeIt;
    for (edgeIt = edgeSet.begin(); edgeIt != edgeSet.end(); ++edgeIt) {
      edges[*edgeIt].traversed = GraphEdge::Marked;
      vertices[edges[*edgeIt].src].traversed = GraphVertex::Marked;
      vertices[edges[*edgeIt].dest].traversed = GraphVertex::Marked;
    }
  }

  template<typename Container>
		void UntraverseEdgeSet(Container &edgeSet) {
		std::set<ssize_t>::iterator edgeIt;
    for (edgeIt = edgeSet.begin(); edgeIt != edgeSet.end(); ++edgeIt) {
      edges[*edgeIt].flagged = GraphEdge::NotMarked;
      vertices[edges[*edgeIt].src].traversed = GraphVertex::NotMarked;
      vertices[edges[*edgeIt].dest].traversed = GraphVertex::NotMarked;
    }
  }
  ssize_t IsSetCyclic(std::set<ssize_t> &vertexSet, ssize_t curVertex, std::list<ssize_t> &path);
  ssize_t IsDAGSuspect(std::set<ssize_t> &dagEdges);
  ssize_t PruneDAG(std::set<ssize_t> &dagVertices, 
							 std::set<ssize_t> &dagEdges,
							 std::vector<ssize_t> &path,
							 std::vector<ssize_t> &verticesToDelete,
							 std::vector<ssize_t> &edgesToDelete);

  ssize_t IsEdgeSuspect(ssize_t e, ssize_t polyNucLength);
  void RerouteCompoundEdgeIntervals(ssize_t vertex, ssize_t inEdge, ssize_t outEdge);
  void ConcatenateEdgeSequence(ssize_t vertex, ssize_t inEdge, ssize_t outEdge);
  void RerouteSimplePathIntervals(ssize_t vertex, ssize_t inEdge, ssize_t outEdge);
  ssize_t  RemoveSimpleBulges(ssize_t minBulgeSize);
  ssize_t  GetOutEdgeIndex(ssize_t vertex, ssize_t edge);
  ssize_t  GetInEdgeIndex(ssize_t vertex, ssize_t edge);
  ssize_t  RemoveSimpleBulge(ssize_t vertex, ssize_t e1, ssize_t e2, ssize_t &edgeIndex);
  ssize_t  RemoveSimpleBulgeLowerMult(ssize_t vertex, ssize_t e1, ssize_t e2, ssize_t &edgeIndex);
  void MergeEdgeIntervals(ssize_t source, ssize_t dest);

  void ClearMST();

  ssize_t FindRouteToMaxPath(ssize_t curVertex, 
												 std::vector<ssize_t> &optPathTraversals, 
												 std::vector<ssize_t> &vertexPath);

  ssize_t FindVertexOnMaxPath(ssize_t curVertex, 
													std::list<ssize_t> &bestPath,
													std::list<ssize_t> &curPath, 
													std::vector<ssize_t> &optPathTraversals);

  ssize_t VertexOnMaxPath(ssize_t curVertex, std::vector<ssize_t> &optPathTraversals);
  ssize_t EdgeOnMaxPath(ssize_t edge, std::vector<ssize_t> &optPathTraversals);
  ssize_t VertexOnAlternatePath(ssize_t curVertex, std::vector<ssize_t> &optPathTraversals);
  ssize_t FindAlternatePaths(ssize_t sourceVertex, 
												 std::vector<ssize_t> &optPathTraversals, 
												 std::vector<ssize_t> &vertexPath,
												 std::vector<ssize_t> &verticesToDelete,
												 std::vector<ssize_t> &edgesToDelete);
	void DeleteEdgeReadIntervals(ssize_t edge);
	void DeleteEdgeListReadIntervals(std::vector<ssize_t> &edgeList);
	void UpdatePathEdges(std::vector<ssize_t> &edgesToRemove);
  ssize_t VertexSetContained(std::set<ssize_t> &dag, ssize_t source);
  void CollectVertices(ssize_t curVertex, ssize_t destVertex, std::set<ssize_t> &dagVertices);
  void CollectEdges(ssize_t lastVertex, std::set<ssize_t> &vertexSet, std::set<ssize_t> &edgeSet);
	ssize_t SearchForCycle(ssize_t sourceVertex, ssize_t prevVertex, ssize_t curVertex, ssize_t curPathLength, ssize_t maxPathLength,
										 std::string padding = "");
};

class RemoveAllButMSTFunctor {
 public:
  TVertexList *vertices;
  TEdgeList   *edges;
  std::vector<ssize_t> erodedVertexList;
  std::vector<ssize_t> erodedEdgeList;
  void operator()(ssize_t vertexIndex) {
    std::cout << "THIS ISN't IMPLEMENTED" << std::endl;
    exit(0);
    ssize_t edgeIndex, balancedIndex;
    ssize_t outEdgeIndex;
    for (outEdgeIndex = (*vertices)[vertexIndex].FirstOut();
				 outEdgeIndex < (*vertices)[vertexIndex].EndOut();
				 outEdgeIndex = (*vertices)[vertexIndex].NextOut(outEdgeIndex)) {
      edgeIndex = (*vertices)[vertexIndex].out[outEdgeIndex];
      balancedIndex = (*edges)[edgeIndex].balancedEdge;
      if ((*edges)[edgeIndex].mst != GraphAlgo::MSTIn) {
				assert((*edges)[balancedIndex].mst != GraphAlgo::MSTIn);
				erodedEdgeList.push_back(edgeIndex);
      }
    }
  }
};


class ErodeLeavesFunctor {
 public:
  TVertexList *vertices;
  TEdgeList   *edges;
	IntervalGraph *graph;
  std::vector<ssize_t> erodedVertexList;
  std::vector<ssize_t> erodedEdgeList;
  ssize_t minLength;
  void operator()(ssize_t vertexIndex);
  void Clear() {
    erodedVertexList.clear();
    erodedEdgeList.clear();
  }
};

void FormGraphFileNames(std::string &base, 
												std::string &bgraph, std::string &graph,
												std::string &intv, 
												std::string &path, std::string &edge,
												std::string &altEdge);


void ReadIntervalGraph(std::string &base, 
											 IntervalGraph &graph, _INT_ &vertexSize,
											 _INT_ skipIntervals = 0,
											 std::ostream &report = std::cout
											 );

void ReadIntervalGraph(std::string &base, 
											 IntervalGraph &graph, _INT_ &vertexSize,
											 std::ostream &report
											 );


void WriteIntervalGraph(std::string &base,
												IntervalGraph &graph, int vertexSize = 1,
												std::ostream &report = std::cout
												);



#endif

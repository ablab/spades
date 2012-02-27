/***************************************************************************
 * Title:          IntervalGraph.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/14/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "HashedSpectrum.h"
#include "IntervalGraph.h"
#include "DeBruijnGraph.h"
#include "graph/GraphAlgo.h"
#include "graph/MSTAlgo.h"
#include "Spectrum.h"
#include "align/alignutils.h"
#include "AlignmentPrinter.h"
#include <iomanip>
#include <algorithm>
#include <iterator>
#include <functional>
#include <iostream>
#include <queue>
#include "compatibility.h"

//#define VERBOSE
#undef VERBOSE

using namespace std;

void FormGraphFileNames(std::string &base, 
												std::string &bgraph, std::string &graph,
												std::string &intv, 
												std::string &path, std::string &edge,
												std::string &altEdge) {
	bgraph = base + ".bgraph";
	graph  = base + ".graph";
	intv   = base + ".intv";
	path   = base + ".path";
	edge   = base + ".edge";
	altEdge= base + ".altEdges";
}


void ReadIntervalGraph(std::string &baseIn, 
											 IntervalGraph &graph, _INT_ &vertexSize,
											 std::ostream &report) {
	ReadIntervalGraph(baseIn, graph, vertexSize, 0, report);
}


void ReadIntervalGraph(std::string &baseIn, 
											 IntervalGraph &graph, _INT_ &vertexSize, _INT_ skipIntervals,
											 std::ostream &report) {

	std::string bGraphFileName, graphFileName,
		intervalFileName, pathFileName, edgeFileName, altEdgeFileName;

	FormGraphFileNames(baseIn, bGraphFileName, graphFileName,
										 intervalFileName, pathFileName, edgeFileName, 
										 altEdgeFileName);

	//	graph.vertexSize    = vertexSize;

	graph.ReadIntervalGraph(bGraphFileName, intervalFileName, pathFileName, skipIntervals, report);

	vertexSize = graph.vertexSize;

	ReadSequences(edgeFileName, graph.edges, report);
	
	// 
	// If alternative edges exist, read them in.
	//
	
	graph.ReadAlternativeEdges(altEdgeFileName.c_str(), report);

}




void WriteIntervalGraph(std::string &baseOut,
												IntervalGraph &graph, int vertexSize,
												std::ostream &report) {

	std::string bGraphFileName, graphFileName,
		intervalFileName, pathFileName, edgeFileName, altEdgeFileName;

	graph.SetMultiplicities(); // TODO: Determine correct place to put this

	FormGraphFileNames(baseOut, bGraphFileName, graphFileName,
										 intervalFileName, pathFileName, edgeFileName,
										 altEdgeFileName);

	graph.PrintIntervalGraph(bGraphFileName, intervalFileName, report);

  PrintGraph(graph.vertices, graph.edges, graphFileName, report);
  PrintEdges(graph.vertices, graph.edges, edgeFileName, report);
	WriteReadPaths(pathFileName, graph.paths, graph.pathLengths, report);

	//
	// Read in the alternative edges (but don't try if the 
	// alt edge file does not exist.)
	//
	graph.WriteAlternativeEdges(altEdgeFileName, report);
	
}


void IntervalGraph::CalcDMST() {
  GraphAlgo::CalcDirectedMST(vertices, edges);
}


void IntervalGraph::ErodeShortEndIntervals(ssize_t minIntvLength) {
	ssize_t e;
	ssize_t numEroded = 0;
	ssize_t changeMade = 0;
	ssize_t prevNumEroded = -1;
	do {
		changeMade = 0;
		for (e = 0; e < edges.size(); e++ ){ 
			ReadIntervalList *intervals = edges[e].intervals;
			ssize_t numIntervals = intervals->size();
			ssize_t i;
			ssize_t pathIndex;
			ssize_t src = edges[e].src;
			ssize_t srcLen = vertices[src].vertexSize;
			for (i = 0; i < numIntervals; i++ ) {
				// Is this interval the last on the path?
				pathIndex = (*intervals)[i].read;
				if ((*intervals)[i].markedForDeletion){
					continue;
				}
				// This path is harmless if it corresponds to 1 or fewer edges.
				if (pathLengths[pathIndex] <= 1) {
					continue;
				}
				if (((*intervals)[i].pathPos == pathLengths[pathIndex]-1) and
						(*intervals)[i].length - srcLen < minIntvLength) {
					// This interval is too short, remove it.

					edges[e].MarkIntervalForRemoval(i);
					pathLengths[pathIndex]--;

					// Erode the reverse complement.
					ssize_t rcPathIndex;
					if (pathIndex % 2 == 0) {
						rcPathIndex = pathIndex + 1;
					}
					else {
						rcPathIndex = pathIndex - 1;
					}
					ssize_t rcIntvEdge  = paths[rcPathIndex][0].edge;
					ssize_t rcIntvIndex = paths[rcPathIndex][0].index;
					edges[rcIntvEdge].MarkIntervalForRemoval(rcIntvIndex);
					ssize_t p;
					for (p = 0; p < pathLengths[rcPathIndex] - 1; p++) {
						paths[rcPathIndex][p] = paths[rcPathIndex][p+1];
						(*edges[paths[rcPathIndex][p].edge].intervals)[paths[rcPathIndex][p].index].pathPos = p;
					}
					pathLengths[rcPathIndex]--;
				
					++numEroded;
					changeMade = 1;
					// Bail out early if only a small change was made.
					if (prevNumEroded != -1 and 
							(prevNumEroded < numEroded * .01))
							//							(1 - double(numEroded - prevNumEroded) / numEroded) < .01)
						changeMade = 0;
				}
			}
		}
		cout << CurTimeString() << ": eroded: " << numEroded << " paths. " << endl;
		RemoveMarkedIntervals();
		prevNumEroded = numEroded;
	}
	while (changeMade);
}

void IntervalGraph::PrintIntervalGraph(std::string &bGraphName,
																			 std::string &intervalName,
																			 std::ostream &report) {
	std::ofstream bGraphOut;
	openck(bGraphName, bGraphOut, std::ios::out, report);
  PrintBGraph(vertices, edges, vertexSize, bGraphOut);
	bGraphOut << vertexSize;
	//	SortAllEdgeIntervalsByEdgePos();
	PrintReadIntervals(edges, intervalName, report);
	//	SortAllEdgeIntervalsByReadPos();
} 

ssize_t IntervalGraph::CalculateReadLength(ssize_t read) {

	if (pathLengths[read] <= 0) {
		return 0;
	}

  ssize_t pathPos;
  pathPos = 0;

  ssize_t readLength = 0;
  ssize_t edgeIndex;
  ssize_t intvIndex;
  ssize_t destVertex = -1;
  for (pathPos = 0; pathPos < pathLengths[read]; pathPos++) {
    edgeIndex = paths[read][pathPos].edge;
    intvIndex = paths[read][pathPos].index;
    destVertex = edges[edgeIndex].dest;
    readLength+= (*edges[edgeIndex].intervals)[intvIndex].length - vertices[destVertex].vertexSize;
  }
  /* 
     The last interval goes up to and including the last vertex, so the loop
     over-subtracts by one vertex.  Add that back now.
  */
  readLength += vertices[destVertex].vertexSize;
  return readLength;
}

ssize_t IntervalGraph::CheckGraphStructureBalance() {
  //UNUSED// ssize_t  v;
  ssize_t e;
  ssize_t destV, balV, balE;
  for (e = 0; e < edges.size(); e++) {
    destV = edges[e].dest;
    balE  = edges[e].balancedEdge;
    balV  = edges[balE].src;
    if (edges[balE].balancedEdge != e) {
      std::cout << "EDGE BALANCE ERROR: " << e << " != " << edges[balE].balancedEdge << std::endl;
    }
    else {
      if (destV >= 0 and 
					(vertices[destV].InDegree() != vertices[balV].OutDegree() or
					 vertices[destV].OutDegree() != vertices[balV].InDegree())) {
				std::cout << "Graph is not complementary. " << std::endl;
				std::cout << "vertex " << destV << " " << vertices[destV].InDegree() 
									<< " " << vertices[destV].OutDegree() << " versus: " << balV <<"  "
									<< vertices[balV].InDegree() << " "
									<< vertices[balV].OutDegree() << std::endl;
				return 0;
      }
    }
  }
	return 1;
}

void IntervalGraph::AssignPathOrderToEdges() {
  ssize_t p, e;
  for (p = 0; p < paths.size(); p++ ) {
    for (e = 0; e < pathLengths[p]; e++ ) {
			assert(paths[p][e].edge < edges.size());
			assert(paths[p][e].index < edges[paths[p][e].edge].intervals->size());
      (*edges[paths[p][e].edge].intervals)[paths[p][e].index].pathPos = e;
    }
  }
}
	
ssize_t IntervalGraph::CheckAllPathsBalance(ssize_t fatal) {
  //UNUSED// ssize_t pos;
  ssize_t p ;
  //UNUSED// ssize_t edge, intv, balEdge, balIntv;

  for (p = 0; p < paths.size(); p+=2 ) {
    ssize_t balp = p + 1;
    if (CheckPathBalance(p, balp) == 0) {
      PrintImbalancedPaths(p, balp);
      assert(!fatal);
    }
  }
	return 1;
}

void IntervalGraph::PrintImbalancedPaths(ssize_t p, ssize_t balp) {
  std::cout << "path: " << p << " is not a mirror image of " << balp << std::endl;
  std::cout << "     ";
  PrintPath(p, std::cout);
  std::cout << "bal: ";
  PrintPathReverse(balp, std::cout);
}


ssize_t IntervalGraph::CheckPathContinuity(ssize_t p) {
  ssize_t pos;
  ssize_t edge, index;
  ssize_t nextEdge, nextIndex;
  for (pos = 0; pos < pathLengths[p]-1; pos++ ) {
    edge = paths[p][pos].edge;
    index = paths[p][pos].index;
    nextEdge = paths[p][pos+1].edge;
    nextIndex = paths[p][pos+1].index;
    if (edge == -1 or nextEdge == -1)
      continue;
    if ((*edges[edge].intervals)[index].readPos +
				(*edges[edge].intervals)[index].length -
				vertices[edges[edge].dest].vertexSize != 
				(*edges[nextEdge].intervals)[nextIndex].readPos) {
			/*      std::cout << "path: " << p << " is off by " <<
				((*edges[edge].intervals)[index].readPos +
				 (*edges[edge].intervals)[index].length -
				 vertices[edges[edge].dest].vertexSize) -
				(*edges[nextEdge].intervals)[nextIndex].readPos << std::endl;
			*/
      return 0;
    }
		ssize_t dest = edges[edge].dest;
		if (nextEdge != edge and vertices[dest].LookupOutIndex(nextEdge) == -1) {
			cout << "path: " << p << " is not continuous" << endl;
			return 0;
		}
  }
  return 1;
}

ssize_t IntervalGraph::CheckAllPathsContinuity(ssize_t fatal) {
  ssize_t p;
  for (p = 0; p < paths.size(); p++) {
    if (!CheckPathContinuity(p)) 
			return 0;
  }
	return 1;
}

ssize_t IntervalGraph::CheckPathBalance(ssize_t p, ssize_t balp) {
  ssize_t pos;
  ssize_t edge, intv, balEdge, balIntv;
  if (pathLengths[p] > 0) {
    if (pathLengths[balp] != pathLengths[p]) {
      return 0;
    }
    else {
      for (pos = 0; pos < pathLengths[p]; pos++) {
				edge = paths[p][pos].edge;
				intv = paths[p][pos].index;
				balEdge = paths[balp][pathLengths[balp] - pos - 1].edge;
				balIntv = paths[balp][pathLengths[balp] - pos - 1].index;
				assert(edge < 0 or intv < 0 or  (*edges[edge].intervals)[intv].length >= 0);
				assert(balEdge < 0 or intv < 0 or (*edges[balEdge].intervals)[balIntv].length >= 0);
				if ((edge < 0 and balEdge >= 0) or
						(edge >= 0 and balEdge < 0)) {
					return 0;
				}
				if (edge >= 0 and balEdge >= 0 and 
						(*edges[edge].intervals)[intv].length != 
						(*edges[balEdge].intervals)[balIntv].length) {
					return 0;
				}
				if (edge >= 0 and balEdge >= 0) {
					if (edge != edges[balEdge].balancedEdge) {
						//						cout << "edge mismatch: " << p << " " << pos << " " << edge << " " << balEdge << " " << edges[balEdge].balancedEdge << endl;
					}
				}
      }
    }
  }
  // No imbalances detected, return 0K
  return 1;
}
	

void IntervalGraph::PrintPath(ssize_t p, std::ostream &pathOut) {
  ssize_t pos, edge, intv;
  pathOut << std::setw(8) <<  p << " len: " << pathLengths[p] << endl;
  for (pos = 0; pos < pathLengths[p]; pos++) {
    edge = paths[p][pos].edge;
    intv = paths[p][pos].index;
    if (edge != -1)
      pathOut << " (" << p << ", " << pos << "), e " << edge << ", i " << edges[edge].index << ", I " 
							<< intv 
							<< ", R " << setw(8) << (*edges[edge].intervals)[intv].readPos 
							<< ", E " << setw(8) << (*edges[edge].intervals)[intv].edgePos 
							<< ", L " << setw(8) << (*edges[edge].intervals)[intv].length << ", v "
							<< vertices[edges[edge].dest].vertexSize << ") " << endl;
    else
      pathOut << " -1 (0) ";
  }
  pathOut << std::endl;
}


void IntervalGraph::PrintPathReverse(ssize_t p, std::ostream &pathOut) {
  ssize_t pos, edge, intv;
  pathOut << std::setw(8) << p << " len: " << pathLengths[p] << endl;
  ssize_t len = pathLengths[p];
  for (pos = 0; pos < pathLengths[p]; pos++) {
    edge = paths[p][len - pos - 1].edge;
    intv = paths[p][len - pos - 1].index;
    if (edge != -1)
      pathOut << " (" << p << ", " << pos << "), e " << edge << ", i " << edges[edge].index << ", I " 
							<< intv 
							<< ", R "  << setw(8) << (*edges[edge].intervals)[intv].readPos 
							<< ", E "  << setw(8) << (*edges[edge].intervals)[intv].edgePos 
							<< ", L "  << setw(8) << (*edges[edge].intervals)[intv].length << ", v "
							<< vertices[edges[edge].dest].vertexSize << ") " << endl;

    else
      pathOut << " -1 (0) ";
  }
  pathOut << std::endl;
}

void IntervalGraph::PrintPaths(std::string pathName,
															 std::ostream &report) {
  std::ofstream pathOut;
  openck(pathName, pathOut, std::ios::out, report);
  ssize_t p;
  //UNUSED// ssize_t pos;
  //UNUSED// ssize_t edge, intv;
  for (p = 0; p < paths.size(); p+=2 ) {
    if (pathLengths[p] > 0) {
      pathOut << "     ";
      PrintPath(p, pathOut);
    }
    ssize_t balp = p+1;
    if (pathLengths[balp] > 0) {
      pathOut << "bal: ";
      PrintPath(balp, pathOut);
    }
  }
  pathOut.close();
}

void IntervalGraph::AssignVertexSizes() {
	ssize_t v;
	for (v = 0; v < vertices.size(); v++) {
		vertices[v].vertexSize = vertexSize;
	}
}

void IntervalGraph::ReadIntervalGraph(std::string &bGraphName,
																			std::string &intervalName, 
																			std::string pathName, 
																			std::ostream &report) {
	ReadIntervalGraph(bGraphName, intervalName, pathName, 0, report);
}

void IntervalGraph::ReadIntervalGraph(std::string &bGraphName,
																			std::string &intervalName, 
																			std::string pathName, 
																			_INT_ skipIntervals,
																			std::ostream &report) {
  // The interval graph is a combination of graph and intervals
  std::ifstream in;
  openck(bGraphName, in, std::ios::in, report);
  ReadVertexList(in, vertices);
  ReadEdgeList(in, edges);
	int graphVertexSize;
	if ((in >> graphVertexSize)) {
		vertexSize = graphVertexSize;
	}
	//  ReadBGraph(bGraphName, vertices, edges, report);

	AssignVertexSizes();
  CheckVertices(vertices, edges);
  CheckEdges(vertices, edges);

	if (skipIntervals) {
		CheckBalance();
		InitializeFlags();
		return;
	}

  maxReadIndex = ReadReadIntervals(intervalName, edges, report);

  // All the methods here work with intervals sorted by read pos
  // but other code uses them sorted by edge pos
  // Restore the order when writing the graph
	// SortAllEdgeIntervalsByReadPos();
	
  // Assign each read interval its order in the path it is part of
  if (pathName != "") {
		std::cout << CurTimeString() << ": Reading graph from " << pathName << std::endl;
    ReadReadPaths(pathName, paths, pathLengths, report);
		std::cout << CurTimeString() << ": Done reading graph." << std::endl;
    AssignPathOrderToEdges();
  }
  else {
    AssignIntervalPathOrder();
  }

  //	AssignEdgesToIntervals();
  CheckBalance();
  assert(CheckAllPathsBalance(1));
  CheckGraphStructureBalance();
  InitializeFlags();
}


void IntervalGraph::MergeEdgeIntervals(ssize_t source, ssize_t dest) {

  //UNUSED// ssize_t sourceLength = edges[source].length;
  //UNUSED// ssize_t destLength   = edges[dest].length;

  ssize_t i;
  ssize_t numDestIntervals = edges[dest].intervals->size();

  edges[dest].intervals->resize(edges[dest].intervals->size() + 
																edges[source].intervals->size());

  ssize_t destPos;

  // Copy over the read intervals
  for (i = 0; i < edges[source].intervals->size(); i++ ) {
    (*edges[dest].intervals)[i + numDestIntervals] = (*edges[source].intervals)[i];
    destPos = i + numDestIntervals;
  }
	//  SortReadIntervalsByReadPos(*edges[dest].intervals);
  ssize_t read, pathPos;
  for (i = 0; i < (*edges[dest].intervals).size(); i++) {
    /*
      The read and pathPos do not change when merging
      the edge intervals, but the original edge, and position
      on the edge do.  Reset all of them since the order
      changes when sorting.
    */
    read = 		(*edges[dest].intervals)[i].read;
    pathPos = (*edges[dest].intervals)[i].pathPos;
    paths[read][pathPos].edge  = dest;
    paths[read][pathPos].index = i;
  }
}

ssize_t IntervalGraph::GetOutEdgeIndex(ssize_t vertex, ssize_t edge) {
  ssize_t e;
  for (e = 0; e < 4; e++ ) {
    if (vertices[vertex].out[e] == edge)
      return e;
  }
  return 4;
}

ssize_t IntervalGraph::GetInEdgeIndex(ssize_t vertex, ssize_t edge) {
  ssize_t e;
  for (e = 0; e < 4; e++ ) {
    if (vertices[vertex].in[e] == edge)
      return e;
  }
  return 4;
}

void IntervalGraph::RemoveAllSimpleBulges(ssize_t minBulgeSize) {
  ssize_t nRemoved = 0;
  ssize_t removeBulgeIter = 0;
  do {
    nRemoved = RemoveSimpleBulges(minBulgeSize);
    ++removeBulgeIter;
    std::cout << CurTimeString() << ": remove simple bulges iteration " << removeBulgeIter << " removed " << nRemoved << " edges " << std::endl;
  }
  while (nRemoved > 0);

}

ssize_t IntervalGraph::RemoveSimpleBulges(ssize_t minBulgeSize) {
  
  // Remove simple bulges, directed cycles of length 2.  This was 
  // mainly written as practice for merging intervals, to see if it 
  // would work, and to eek out some better results for component size.
  ssize_t v, a, b;
  std::vector<ssize_t> removedEdges;
  ssize_t removedEdge;
  Unflag();
  //UNUSED// ssize_t edgeA, edgeB;
  ssize_t balancedEdge;
  ssize_t edgeIndex;
  
  for (v = 0; v < vertices.size(); v++ ){
		// Find which edges go into the same position
		//    if (vertices[v].OutDegree() == 2) {
		for (a = 0; a < 3 ; a++ ){
			if (vertices[v].out[a] == -1
					or
					edges[vertices[v].out[a]].length >= minBulgeSize)
				continue;
			for (b = a + 1; b < 4 ; b++ ) {
				if (
						//vertices[v].out[a] != -1 and
						vertices[v].out[b] != -1 and
						edges[vertices[v].out[a]].dest == 
						edges[vertices[v].out[b]].dest) {
					if (
							// edges[vertices[v].out[a]].length < minBulgeSize and
							edges[vertices[v].out[b]].length < minBulgeSize) {
						// There is a simple bulge.  Need to check if there is a preference on 
						// which one to remove.  There is a preference when the balanced
						// edge of one of the edges has already been removed.
						if (edges[vertices[v].out[a]].flagged == GraphEdge::Marked) {
							removedEdge = RemoveSimpleBulge(v, vertices[v].out[b], vertices[v].out[a], edgeIndex);
							removedEdges.push_back(removedEdge);
						}
						else if (edges[vertices[v].out[b]].flagged == GraphEdge::Marked) {
							removedEdge = RemoveSimpleBulge(v, vertices[v].out[a], vertices[v].out[b], edgeIndex);
							removedEdges.push_back(removedEdge);
						}
						else {
							// There is no preference on which edge to merge into the other.
							// Keep the edge of higher multiplicity
							if (edges[vertices[v].out[a]].balancedEdge == vertices[v].out[b]) {
								// The bulge involves two balanced edges.  Since we'll be deleting
								// the balance, we mark the balance of the originals to themselves 
								// so that the bookkeeping is correct (no edges remain without balanced
								// edges).
								edges[vertices[v].out[a]].balancedEdge = vertices[v].out[a];
								edges[vertices[v].out[b]].balancedEdge = vertices[v].out[b];
							}
							// Remove the edge of lower multiplicity
							removedEdge = RemoveSimpleBulgeLowerMult(v, vertices[v].out[a], vertices[v].out[b], edgeIndex);
							removedEdges.push_back(removedEdge);
							balancedEdge = edges[removedEdge].balancedEdge;
							// std::cout << "removing bulged edge: " << removedEdge << " and giving preference to " << balancedEdge << std::endl;
							edges[balancedEdge].flagged = GraphEdge::Marked;
							//  std::cout << "removing edge " << removedEdge << " and giving preference to " << balancedEdge << std::endl;
						}
					}
				}
			}
		}
		//    }
	}
	// Now do a sanity check on the removed edges.  For every edge 
	// in the removed set, it's balanced edge should have been removed.

	//UNUSED// ssize_t e2;
	ssize_t e ;
	std::sort(removedEdges.begin(), removedEdges.end());
	for (e = 0; e < removedEdges.size(); e++) {
		removedEdge = removedEdges[e];
		balancedEdge = edges[removedEdge].balancedEdge;
		if (std::binary_search(removedEdges.begin(), removedEdges.end(),
													 balancedEdge) == 0) {
			std::cout << "error, removed " << removedEdge
								<< "( " << edges[removedEdge].index << ") but not "  
								<< edges[removedEdge].balancedEdge << std::endl;
		}
	}

	std::cout << CurTimeString() << ": Simple bulge removal removed: " << removedEdges.size() << " edges " << std::endl;
	std::vector<ssize_t> nullVertices;
	/* 
		 Remove the bulges from the graph.
		 At this point in time the paths are marked for deletion, but they are not 
		 removed so that they may be re-routed.
	*/

	Prune(nullVertices, removedEdges);

	/* 
		 For now re-route before condensing simple paths.  If it were done afterwards, there 
		 would be gapped intervals being joined together. That's not so bad, but just not 
		 yet figured out.
	*/
	//	RouteRemovedIntervals();
	SortAllEdgeIntervalsByReadPos();
	UpdateAllPathIndices();
	CondenseSimplePaths();
	assert(CheckAllPathsBalance(1));
	return removedEdges.size();
}



ssize_t IntervalGraph::IsEdgeSuspect(ssize_t e, ssize_t polyNucLength) {
	// edge e 
	ssize_t p;
	if (polyNucLength > edges[e].length)
		return 0;

	ssize_t polyNucPos = 0;
	ssize_t lastPolyNuc = -1;
	ssize_t polyNucFound = 0;
	for (p = 1; p < edges[e].length; p++ ) {
		if (edges[e].seq.seq[polyNucPos] != edges[e].seq.seq[p])
			polyNucPos = p;
		else if (p - polyNucPos >= polyNucLength) {
			polyNucFound = 1;
			lastPolyNuc = p;
		}
		if (p > vertexSize ) {
			// There exists a k-mer here that does not have 
			// a homopolymer in it.  There fore this edge is trusted
			if (lastPolyNuc < p - vertexSize)
				return 0;
		}
	}
	return 1;
}

ssize_t IntervalGraph::LookupAlignedPosition(ssize_t edge, ssize_t edgePos, ssize_t *alignment, ssize_t seqLength, 
																				 vector<ssize_t> &path,
																				 vector<ssize_t> &edgeStarts, ssize_t traversal) {
	ssize_t nTraversal = 0;
	//UNUSED// ssize_t pathLength = edgeStarts.size();
	vector<ssize_t>::iterator pathIt, pathEnd;
	pathEnd = path.end();
	ssize_t pathPos;
	for(pathPos = 0, pathIt = path.begin();
			pathIt != pathEnd; ++pathIt, ++pathPos) {
		if (*pathIt == edge) {
			++nTraversal;
			if (nTraversal == traversal) {
				// Found the proper path index
				break;
			}
		}
	}
	assert(pathPos < edgeStarts.size());

	ssize_t alignPos = edgeStarts[pathPos];
	alignPos += edgePos;
	assert(alignPos < seqLength);
	return alignPos;
}

ssize_t IntervalGraph::AlignmentBoundariesToPath(ssize_t *alignment, ssize_t alignmentLength,
																						 ssize_t start, ssize_t end,
																						 vector<ssize_t> &path, vector<ssize_t> &edgeStarts, 
																						 vector<ssize_t> &pathEdges, 
																						 vector<ssize_t> &pathEdgeStarts, vector<ssize_t> &pathEdgeLengths) {
	// Find the edge containing the start
	ssize_t edgeIndex;
	ssize_t startEdgeIndex = 0;
	ssize_t endEdgeIndex = 0;
	ssize_t startEdgePos =  0;
	//UNUSED// ssize_t startEdgePathLength = 0;
	startEdgeIndex = -1;

	if (start >= edgeStarts[edgeStarts.size() - 1]) {
		startEdgeIndex = edgeStarts.size() - 1;
	}
	else {
		for (edgeIndex = 0; edgeIndex < edgeStarts.size(); edgeIndex++) {
			if (start < edgeStarts[edgeIndex]) {
				startEdgeIndex = edgeIndex - 1;
				break;
			}
		}
	}

	startEdgePos = start - edgeStarts[startEdgeIndex];
	endEdgeIndex = startEdgeIndex;
	/*	if (end > edgeStarts[edgeStarts.size() - 1] + edges[path[path.size()-1]].length) {
		endEdgeIndex = edgeStarts.size() - 1;
		}
		else {
	*/
	for (edgeIndex = startEdgeIndex; edgeIndex < edgeStarts.size(); edgeIndex++) {
		if (end < edgeStarts[edgeIndex] + edges[path[edgeIndex]].length) {
			endEdgeIndex = edgeIndex;
			break;
		}
	}
	//	}

	ssize_t pathLength = endEdgeIndex - startEdgeIndex + 1;
	assert(pathLength > 0);
	pathEdges.resize(pathLength);
	pathEdgeStarts.resize(pathLength);
	pathEdgeLengths.resize(pathLength);

	// Store the edges corresponding to this path.
	ssize_t i;
	for (i = 0; i < pathLength; i++) {
		pathEdges[i] = path[startEdgeIndex + i];
	}

	// Now store the boundaries corresponding to this path.
	pathEdgeStarts[0] = start - edgeStarts[startEdgeIndex];
	if (pathLength == 1) {
		pathEdgeLengths[0] = end - start + 1;
	}
	else {
		pathEdgeLengths[0] = edges[path[startEdgeIndex]].length  - pathEdgeStarts[0];

		for (i = 1; i < pathLength; i++) {
			// This is an internal interval, so it starts at the beginning of the edge.
			pathEdgeStarts[i] = 0;
			pathEdgeLengths[i] = edges[path[startEdgeIndex + i]].length;
		}
		// Fix the end length if the interval ends before the end of the edge.
		pathEdgeLengths[pathLength-1] = (end - edgeStarts[endEdgeIndex] + 1);
	}
	return pathLength;
}


void IntervalGraph::MarkSuspectEdges(ssize_t polyNucLength) {
	ssize_t e;
	Unsuspect();
	if (vertexSize <= 0) {
		std::cout << "ERROR, the vertex size needs to be set to mark " 
							<< "suspect edges! " << std::endl;
		exit(0);
	}

	//  std::cout << "marking suspect edges!! " << std::endl;
	for (e = 0; e < edges.size(); e++ ) {
		if (IsEdgeSuspect(e, polyNucLength)) {
			edges[e].suspect = GraphEdge::Marked;
			/*
				std::cout << " edge : " << e << " is suspect " << std::endl;
			*/
		}
		else
			edges[e].suspect = GraphEdge::NotMarked;
	}
}

ssize_t IntervalGraph::RemoveSimpleBulgeLowerMult(ssize_t vertex, ssize_t e1, ssize_t e2, ssize_t &edgeIndex) {

	// Merge the edge of lower multiplicity into the edge of higher.
	ssize_t sourceEdge, destEdge;
	ssize_t destVertex;
	destVertex = edges[e1].dest;

	sourceEdge = e1;
	destEdge   = e2;

	// check the lower mult. edge
	if (edges[e2].multiplicity > edges[e1].multiplicity) {
		sourceEdge = e2;
		destEdge = e1;
	}

	RemoveSimpleBulge(vertex, destEdge, sourceEdge, edgeIndex);

	return sourceEdge;
}

ssize_t IntervalGraph::RemoveSimpleBulge(ssize_t vertex, ssize_t destEdge, ssize_t sourceEdge, ssize_t &edgeIndex) {

	// Merge the edge of lower multiplicity into the edge of higher.
	ssize_t destVertex;
	destVertex = edges[destEdge].dest;

	MergeEdgeIntervals(sourceEdge, destEdge);

	//  RemoveEdge(sourceEdge);
	// Unlink the out edge
	edgeIndex = GetOutEdgeIndex(vertex, sourceEdge);
	vertices[vertex].out[edgeIndex] = -1;

	// Unlink the in edge
	edgeIndex  = GetInEdgeIndex(destVertex, sourceEdge);
	vertices[destVertex].in[edgeIndex] = -1;

	/* 
		 Remove this for now since we'll just merge the edge intervals since 
		 the source and dest edges do not change.
		 MarkIntervalsInEdgeForRemoval(sourceEdge);
	*/
	return sourceEdge;
}

class CountAndStoreVertices {
public:
	ssize_t size;
	ssize_t length;
	ssize_t maxSize;
	std::vector<TEdge> *edges;
	std::vector<TVertex> *vertices;
	std::vector<ssize_t> componentIndices, edgeIndices;

	void operator()(ssize_t vertexIndex) {
		ssize_t e;
		size++;
		if (maxSize == -1 and componentIndices.size() > maxSize) {
			return;
		}
		componentIndices.push_back(vertexIndex);
		ssize_t edgeIndex;
		ssize_t srcIndex;
		//UNUSED// ssize_t destIndex;
		for (e = (*vertices)[vertexIndex].FirstIn(); 
				 e < (*vertices)[vertexIndex].EndIn(); 
				 e = (*vertices)[vertexIndex].NextIn(e) ) {
			assert((*vertices)[vertexIndex].in[e] != -1);
			edgeIndex = (*vertices)[vertexIndex].in[e];
			srcIndex = (*edges)[edgeIndex].src;
			if (edgeIndex >= 0 and
					(*edges)[edgeIndex].marked == GraphEdge::NotMarked ) {
				edgeIndices.push_back(edgeIndex);
				length += (*edges)[edgeIndex].length;
				(*edges)[edgeIndex].marked = GraphEdge::Marked;
			}
		}

		for (e = (*vertices)[vertexIndex].FirstOut(); 
				 e < (*vertices)[vertexIndex].EndOut(); 
				 e = (*vertices)[vertexIndex].NextOut(e) ) {
			assert((*vertices)[vertexIndex].out[e] != -1);
			edgeIndex = (*vertices)[vertexIndex].out[e];
			if (edgeIndex >= 0 and
					(*edges)[edgeIndex].marked == GraphEdge::NotMarked) {
				edgeIndices.push_back(edgeIndex);
				length += (*edges)[edgeIndex].length;
				(*edges)[edgeIndex].marked = GraphEdge::Marked;
			}
		}
	}
};

void IntervalGraph::RemoveBalPreferred() {
	ssize_t e;
	for (e = 0; e < edges.size(); e++) 
		edges[e].balPreferred = GraphEdge::NotMarked;
}

void IntervalGraph::Unsuspect(){
	ssize_t e;
	for (e = 0; e < edges.size(); e++) 
		edges[e].suspect = GraphEdge::NotMarked;
}

void IntervalGraph::Unflag() {
	ssize_t v, e;
	for (v = 0; v < vertices.size(); v++ )
		vertices[v].flagged = GraphVertex::NotMarked;
	for (e = 0; e < edges.size(); e++)
		edges[e].flagged = GraphEdge::NotMarked;
}

void IntervalGraph::Untraverse() {
	ssize_t v, e;
	for (v = 0; v < vertices.size(); v++ )
		vertices[v].traversed = GraphVertex::NotMarked;
	for (e = 0; e < edges.size(); e++)
		edges[e].traversed = GraphEdge::NotMarked;
}

void IntervalGraph::Unmark() {
	ssize_t v, e;
	for (v = 0; v < vertices.size(); v++ )
		vertices[v].marked = GraphVertex::NotMarked;
	for (e = 0; e < edges.size(); e++)
		edges[e].marked = GraphEdge::NotMarked;
}

void IntervalGraph::InitializeFlags() {
	ssize_t v, e;
	for (v = 0; v < vertices.size(); v++ )
		vertices[v].Unmark();
	for (e = 0; e < edges.size(); e++)
		edges[e].Unmark();
}

ssize_t IntervalGraph::RemoveSmallComponents(ssize_t length, 
																				 std::string &readsFile,
																				 std::string &componentReadsFile,
																				 std::ostream &report) {

	ssize_t v, cv, ce;
	CountAndStoreVertices countSize;
	countSize.vertices = &vertices;
	countSize.edges    = &edges;
	countSize.maxSize  = length;
	Unmark();
	std::vector<ssize_t> verticesToRemove, edgesToRemove;
	ssize_t component = 0;
	for (v = 0; v < vertices.size(); v++ ) {
		if (vertices[v].marked == GraphVertex::NotMarked ) {
			countSize.size = 0;
			countSize.length = 0;
			countSize.componentIndices.clear();
			countSize.edgeIndices.clear();
			TraverseDFS(vertices, edges, v, countSize);
			if (countSize.length < length) {
				// Need to push back everything in the component
				for (cv = 0; cv < countSize.componentIndices.size(); cv++) {
					verticesToRemove.push_back(countSize.componentIndices[cv]);
				}
				for (ce = 0; ce < countSize.edgeIndices.size(); ce++) {
					edgesToRemove.push_back(countSize.edgeIndices[ce]);
				}
			} 
			++component;
		}
	}
	std::cout << CurTimeString()
						<< ": component removal removed: " << verticesToRemove.size() << " vertices  / "
						<< vertices.size() << " and : " 
						<< edgesToRemove.size() << " edges / " 
						<< edges.size() << std::endl;

	ssize_t e, i;
	MarkIntervalsInEdgeForRemoval(edgesToRemove);
	SimpleSequenceList reads;
	std::ofstream compReadsOut;
	if (readsFile != "") {
		std::cout << CurTimeString() << ": Printing small component reads to: " << componentReadsFile << std::endl;
		openck(componentReadsFile, compReadsOut, std::ios::out, report);
		ReadSimpleSequences(readsFile, reads, report);
		AppendReverseComplements(reads);
		ssize_t er;
		ssize_t read;
		std::stringstream titlestrm;
		for (er = 0; er < edgesToRemove.size(); er++ ) {
			e = edgesToRemove[er];
			for (i = 0; i < (*edges[e].intervals).size(); i++ ) {
				read = (*edges[e].intervals)[i].read;
				titlestrm.str("");
				titlestrm << read;
				if (pathLengths[read] > 0) {
					reads[read].PrintSeq(compReadsOut, titlestrm.str());
					pathLengths[read] = 0;
				}
			}
		}
		compReadsOut.close();
	}

	Prune(verticesToRemove, edgesToRemove);
	ssize_t numRemoved = RemoveMarkedIntervals();
	std::cout << CurTimeString() << ": after removing components, removed " << numRemoved << " paths " << std::endl;


	for (e = 0; e < edges.size();e++ ){ 
		for (i = 0; i < (*edges[e].intervals).size(); i++) { 
			assert((*edges[e].intervals)[i].pathPos >= 0);
		}
	}
	assert(CheckAllPathsBalance(1));
	return verticesToRemove.size();
}


void IntervalGraph::RemoveEdges(std::vector<ssize_t> &edgesToRemove, std::vector<ssize_t> &orphanedVertices) {
	ssize_t e;
	for (e = 0; e < edgesToRemove.size(); e++ ){ 
		RemoveEdge(edgesToRemove[e], orphanedVertices);
	}
}

void IntervalGraph::ConcatenateEdgeSequence(ssize_t vertex, ssize_t inEdge, ssize_t outEdge) {

	// Copy the edge from outEdge into in edge
	// and update the edge lengths
	ssize_t oldEdgeLength = edges[inEdge].length;
	edges[inEdge].seq.length += edges[outEdge].seq.length - vertexSize;

	unsigned char* newSeq;
	// keep the old sequence
	newSeq = new unsigned char[edges[inEdge].seq.length];
	memcpy(newSeq, edges[inEdge].seq.seq, oldEdgeLength);

	// Make sure we are not about to write over memory 
	// starting before the end of the edge
	assert(edges[inEdge].seq.length > vertexSize);

	memcpy((unsigned char*) (newSeq + oldEdgeLength - vertexSize),
				 edges[outEdge].seq.seq, edges[outEdge].seq.length);

	delete[] edges[inEdge].seq.seq;
	edges[inEdge].seq.seq = newSeq;
	edges[inEdge].length = edges[inEdge].seq.length;
}

void IntervalGraph::MergeSimplePath(ssize_t vertex, ssize_t inEdge, ssize_t outEdge) {

	// do a couple of sanity checks.
	assert(vertices[vertex].InDegree() == 1);

	// Add outedge.seq to inedge.seq
	ConcatenateEdgeSequence(vertex, inEdge, outEdge);


	// Update the length of the edge, and the sequence

	// Modify the adjacencies of vertex so that inEdge.dest is outEdge.dest
	SkipOutEdge(vertex, inEdge, outEdge);

	// Mark this edge as not used
	edges[outEdge].dest = -1;
	edges[outEdge].src  = -1;

	// Mark this vertex as not used.
	ssize_t inIndex, outIndex;
	inIndex  = vertices[vertex].LookupInIndex(inEdge);
	outIndex = vertices[vertex].LookupOutIndex(outEdge);
	assert(inIndex >= 0);
	//	assert(outIndex >= 0);
	vertices[vertex].in[inIndex]   = -1;
	if (outIndex >= 0)
		vertices[vertex].out[outIndex] = -1;

	//  SortReadIntervalsByReadPos(*edges[inEdge].intervals);
	//  UpdatePathIndices(inEdge);

	// Free up some memory.
	edges[outEdge].Clear();
	//	delete edges[outEdge].intervals;

	// Propagate some flags.
	if (edges[outEdge].guarded == GraphEdge::Marked)
		edges[inEdge].guarded = GraphEdge::Marked;
}

void IntervalGraph::UpdateAllPathIndices() {
	ssize_t e;
	for (e = 0; e < edges.size(); e++) 
		UpdatePathIndices(e);
}

void IntervalGraph::UpdatePathIndices(ssize_t edge) {
	/*
	 * Make path intervals point towards the correct index in an edge.
	 */
	ssize_t i;
	ssize_t read, pathPos;
	for (i = 0; i < edges[edge].intervals->size(); i++) { 
		if (!IsIntervalMarkedForRemoval(edge, i)) {
			read    = (*edges[edge].intervals)[i].read;
			pathPos = (*edges[edge].intervals)[i].pathPos;
			assert(read >= 0);
			assert(read < paths.size());
			assert(pathPos >= 0);
			assert(pathPos < pathLengths[read]);
			paths[read][pathPos].index = i;
		}
	}
}

void IntervalGraph::SortAllEdgeIntervalsByReadPos() {
	ssize_t e;
	for (e = 0; e < edges.size(); e++) {
		SortReadIntervalsByReadPos(*edges[e].intervals);
	}
}


void IntervalGraph::SortAllEdgeIntervalsByEdgePos() {
	ssize_t e;
	for (e = 0; e < edges.size(); e++) {
		SortReadIntervals(*edges[e].intervals);
	}
}

void IntervalGraph::SkipOutEdge(ssize_t vertex, ssize_t inEdge, ssize_t outEdge) {

	// Transform part of the graph of 
	//     the form  src--inEdge-->vertex--outedge-->dest
	//     to   --inEdge->dest
	//     
	// In edge now becomes the balance of outEdge->balancedEdge.
	// 
	// Unlink outEdge from graph
	ssize_t destVertex;
	ssize_t srcVertex;

	srcVertex  = edges[inEdge].src;
	destVertex = edges[outEdge].dest;

	// Fix the destination of the in edge so that 
	// it skips 'vertex'
	edges[inEdge].dest = destVertex;

	// Now fix the in edge from the dest vertex
	// have   vertex -- outEdge --> destVertex
	// we want the index of this.
	ssize_t inEdgeIndex = vertices[destVertex].LookupInIndex(outEdge);
	assert(inEdgeIndex != vertices[destVertex].EndIn());

	// point back past the skipped edge.
	vertices[destVertex].in[inEdgeIndex] = inEdge;

	// Fix the balanced edge of the skipped edge to reference the 
	// in edge.
	// We hope that only balanced operations will be performed on the graph

	ssize_t balOutEdge;
	balOutEdge = edges[outEdge].balancedEdge;

	edges[balOutEdge].balancedEdge = inEdge;

	ssize_t balInEdge;
	balInEdge = edges[inEdge].balancedEdge;
	if (edges[balInEdge].balancedEdge == -1) {
		// The in-edge references an edge that has been removed.
		// This can only happen if balanced edges have already been
		// skipped.  So, it's safe to make inedge point to the balanced
		// edge that the out edge uses, since it's already been merged with that.
		// I should draw an ascii diagram for this, but I don't want to now.
		edges[inEdge].balancedEdge = edges[outEdge].balancedEdge;
	}
	edges[outEdge].balancedEdge = -1;
}


ssize_t IntervalGraph::CondenseSimplePaths() {
	std::cout << CurTimeString() << ": CondenseSimplePaths()" << std::endl;

	std::vector<ssize_t> removedVertices;
	std::vector<ssize_t> removedEdges;

	ssize_t vertex;
	ssize_t destVertex;
	ssize_t outEdge, outEdgeIndex;
	ssize_t followingVertex, followingEdge, followingEdgeIndex;
	ssize_t prevLength;

	// Phase 1, modify read intervals that are condensed.
	vertex = 0;
	for (vertex = 0; vertex < vertices.size(); vertex++) { 
		// Found a branching vertex, attempt to simplify paths from it.
		if (vertices[vertex].OutDegree() != 1 or vertices[vertex].InDegree() != 1) {
			for (outEdgeIndex = vertices[vertex].FirstOut();
					 outEdgeIndex < vertices[vertex].EndOut();
					 outEdgeIndex = vertices[vertex].NextOut(outEdgeIndex)) {
				outEdge    = vertices[vertex].out[outEdgeIndex];
				//UNUSED// ssize_t destVertex0;
				destVertex = edges[outEdge].dest;
				prevLength = edges[outEdge].length;
				
				ssize_t beginSrc = edges[outEdge].src;
				followingEdge = outEdge;
				//UNUSED// ssize_t lengthOffset0,curIntv0;
				ssize_t lengthOffset,  curIntv ;
				lengthOffset = edges[outEdge].length - vertices[destVertex].vertexSize;
				curIntv = edges[outEdge].intervals->size();

				//
				// Count how many new intervals will be added.
				//
				ssize_t numSimplePathIntervals = CountIntervalsOnSimplePath(outEdge);
				if (numSimplePathIntervals > edges[outEdge].intervals->size()) {

					//
					// Only bother adding all the new intervals if more exist.
					//
					edges[outEdge].intervals->resize(numSimplePathIntervals);
					while (destVertex != beginSrc and 
								 vertices[destVertex].OutDegree() == 1 and
								 vertices[destVertex].InDegree() == 1) {

						// followingEdge is the next edge on the simple path
						followingEdgeIndex = vertices[destVertex].FirstOut();
						followingEdge      = vertices[destVertex].out[followingEdgeIndex];
						destVertex         = edges[followingEdge].dest;
						ssize_t numFollowingEdgeIntervals = edges[followingEdge].intervals->size();
						// Make the intervals on the following edge reference out edge
						MoveIntervals(outEdge, followingEdge, curIntv, lengthOffset);

						curIntv      += numFollowingEdgeIntervals;
						lengthOffset += edges[followingEdge].length - vertices[destVertex].vertexSize;
					}
					//					GrowIntervalsOnSimplePath(followingEdge);
				}
			}
		}
	}	

	// Phase 2. Remove all intervals that were marked as continuations
	// in the previous step.
	//	RemoveMarkedIntervals();
	//	RemoveMarkedPathIntervals();

	// Phase 3. Combine intervals along the path into the original edge

	// Phase 4. Merge the edges along the simple path
	// This could be sped up as well.


	// Phase 5. Merge adjacent edges on the simple path.

	for (vertex = 0; vertex < vertices.size(); vertex++) {
		// If this vertex is the beginning of a path, try to simplify it.
		if (vertices[vertex].OutDegree() != 1 or vertices[vertex].InDegree() != 1) {
			for (outEdgeIndex = vertices[vertex].FirstOut();
					 outEdgeIndex < vertices[vertex].EndOut();
					 outEdgeIndex = vertices[vertex].NextOut(outEdgeIndex)) {
				outEdge = vertices[vertex].out[outEdgeIndex];
				destVertex = edges[outEdge].dest;
				prevLength = edges[outEdge].length;

				while (vertices[destVertex].OutDegree() == 1 and
							 vertices[destVertex].InDegree() == 1) {

					// Found an vertex --> destVertex --> some other vertex
					// Find out where this vertex goes.
					followingEdgeIndex = vertices[destVertex].FirstOut();
					followingEdge      = vertices[destVertex].out[followingEdgeIndex];
					followingVertex    = edges[followingEdge].dest;

					ssize_t curOutEdgeLength = edges[outEdge].length;
					MergeSimplePath(destVertex, outEdge, followingEdge);

					// Remove this vertex and the outEdge
					removedVertices.push_back(destVertex);
					removedEdges.push_back(followingEdge);

					// Re-assign the destination edge
					destVertex = followingVertex;
					AppendAlternativeEdges(edges[followingEdge].altEdges,
																 curOutEdgeLength - 
																 vertices[edges[outEdge].dest].vertexSize,
																 edges[outEdge].altEdges);
					edges[followingEdge].altEdges.clear();
																 
				}
			}
		}
	}

	Prune(removedVertices, removedEdges);

	std::cout << CurTimeString() << ": Exit CondenseSimplePaths()" << std::endl;

	return removedEdges.size();
}

ssize_t IntervalGraph::CollectSimplePathIntervals(ssize_t startEdge, ssize_t numNewIntervals) {

	ssize_t destVertex = edges[startEdge].dest;
	ssize_t curIntv = edges[startEdge].intervals->size();
	ssize_t storeIntervals = 0;
	ssize_t prevNumIntervals = edges[startEdge].intervals->size();
	if (numNewIntervals > 0) {
		storeIntervals = 1;
		numNewIntervals += prevNumIntervals;
		edges[startEdge].intervals->resize(numNewIntervals );
	}
	ssize_t nextEdge, nextEdgeIndex;
	ssize_t firstSource = edges[startEdge].src;
	ssize_t edgeOffset  = edges[startEdge].length - vertices[destVertex].vertexSize;
	while (destVertex != firstSource and
				 vertices[destVertex].OutDegree() == 1 and
				 vertices[destVertex].InDegree() == 1) {
		nextEdgeIndex = vertices[destVertex].FirstOut();
		nextEdge      = vertices[destVertex].out[nextEdgeIndex];
		ssize_t nextIntv;
		ssize_t nextIntvPath, nextIntvPathPos;
		for (nextIntv = 0; nextIntv < edges[nextEdge].intervals->size(); nextIntv++ ) { 
			if (!IsIntervalMarkedForRemoval(nextEdge, nextIntv)) {
				if (storeIntervals) {
					// We should only be adding edges that start in the next edge.c
					//			assert((*edges[nextEdge].intervals)[nextIntv].pathPos == 0);
					assert(curIntv < numNewIntervals);
					(*edges[startEdge].intervals)[curIntv] = (*edges[nextEdge].intervals)[nextIntv];

					// Fix offsets
					(*edges[startEdge].intervals)[curIntv].edgePos += edgeOffset;
					nextIntvPath = (*edges[nextEdge].intervals)[nextIntv].read;
					nextIntvPathPos = (*edges[nextEdge].intervals)[nextIntv].pathPos;
					paths[nextIntvPath][nextIntvPathPos].edge = startEdge;
					++curIntv;
				}
				else {
					// Simply counting them this time.  Will store on the next pass
					// so that the array may be resized only once.
					++numNewIntervals;
				}
			}
		}
		destVertex = edges[nextEdge].dest;
		edgeOffset += edges[nextEdge].length - vertices[destVertex].vertexSize;
	}

	if (storeIntervals) {
		assert(curIntv == numNewIntervals);
		SortReadIntervalsByReadPos(*edges[startEdge].intervals);
		UpdatePathIndices(startEdge);
	} 

	return numNewIntervals;
}

void IntervalGraph::GrowIntervalsOnSimplePath(ssize_t edgeIndex) {
	ssize_t intv;
	ssize_t numIntervals = edges[edgeIndex].intervals->size();
	ssize_t pathIndex, pathPos;
	for (intv = 0; intv < numIntervals; intv++ ){ 
		if (!IsIntervalMarkedForRemoval(edgeIndex, intv)) {
			// search ahead for this path to be extended.
			pathIndex = (*edges[edgeIndex].intervals)[intv].read;
			pathPos   = (*edges[edgeIndex].intervals)[intv].pathPos;
			ssize_t p;
			ssize_t isSimplePath = 1;
			ssize_t destVertex = edges[edgeIndex].dest;
			ssize_t nextIntvEdge, nextIntvPos;

			for (p = pathPos + 1; p < pathLengths[pathIndex] and isSimplePath; p++ ) {
				nextIntvEdge = paths[pathIndex][p].edge;
				nextIntvPos  = paths[pathIndex][p].index;
				if (nextIntvEdge < 0 or nextIntvPos < 0)
					continue;
				if ((*edges[nextIntvEdge].intervals)[nextIntvPos].markedForDeletion)
					continue;

				assert(nextIntvEdge >= 0);
				assert(nextIntvPos >= 0);
				(*edges[edgeIndex].intervals)[intv].length += 
					(*edges[nextIntvEdge].intervals)[nextIntvPos].length - 
					vertices[destVertex].vertexSize;

				destVertex = edges[nextIntvEdge].dest;

				MarkIntervalForRemoval(nextIntvEdge, nextIntvPos);

				// Look to see of more of this paht should be extended
				if (vertices[destVertex].InDegree() == 1 and 
						vertices[destVertex].OutDegree() == 1) {
					isSimplePath = 1;
				}
				else {
					isSimplePath = 0;
				}
			}
		}
	}
}

void IntervalGraph::MarkPathsThroughEdgeForRemoval(ssize_t edgeIndex) {
	std::vector<ssize_t> pathEdges, pathIndices;
	ssize_t i;
	// Note: I check to see if an interval is traversed so that
	// I don't re-mark paths for removal.  There is nothing terribly
	// wrong with doing that, and so to save space, the traversed
	// flag may be removed.
	for (i = 0; i < edges[edgeIndex].intervals->size(); i++) {
		if ((*edges[edgeIndex].intervals)[i].read >= 0 and
				(*edges[edgeIndex].intervals)[i].traversed == 0) {
			// This edge has not yet been flagged for removal, do it now
			ssize_t readIndex = (*edges[edgeIndex].intervals)[i].read;
			MarkIntervalsOnPathForRemoval(readIndex);
			delete[] paths[readIndex];
			paths[readIndex] = NULL;
			pathLengths[readIndex] = 0;
			//      TracePath(edgeIndex, i, pathEdges, pathIndices);
			//      MarkPathForRemoval(pathEdges, pathIndices);
		}
	}
}

void IntervalGraph::MarkIntervalsInEdgeForRemoval(std::vector<ssize_t> &edgeIndices) {
	ssize_t e;
	for (e = 0; e < edgeIndices.size(); e++) {
		MarkIntervalsInEdgeForRemoval(edgeIndices[e]);
	}
}

void IntervalGraph::MarkIntervalsInEdgeForRemoval(ssize_t edgeIndex) {
	ssize_t i;
	//UNUSED// ssize_t read;
	//UNUSED// ssize_t pathPos;
	/* 
		 Given an edge, mark all intervals in it for removal.
	*/ 
	//	std::cout << "marking intervals in edge " << edgeIndex << std::endl;
	for (i = 0; i < edges[edgeIndex].intervals->size(); i++ ) { 
		MarkIntervalForRemoval(edgeIndex, i);
	}
}

void IntervalGraph::DiscardGappedPaths() {

	/* 
		 Discard paths that have gaps in them (a path is a series of edges and indexes
		 on the edge ( (e1,i1), (e2,i2), ... , (en,in) ).
		 The path is gapped if one of the pairs is cleared:
		 ((e1,i1), (-1,-1), ... (en,in)).

	*/
	ssize_t p, pi;
	ssize_t pathIsGapped;
	for (p = 0; p < paths.size(); p++) {
		pathIsGapped = 0;
		for (pi = 0; pi < pathLengths[p]; pi++) {
			/* Either both are not -1, or both are */
			if (paths[p][pi].edge == -1 and 
					paths[p][pi].index == -1) {
				pathIsGapped = 1;
				break;
			}
			else {
				assert(paths[p][pi].edge != -1 and paths[p][pi].index != -1);
			}
		}
		if (pathIsGapped) {
			/*
				std::cout << "path: " << p  << " of len: " << pathLengths[p] << " is gapped " << std::endl;
				for (pi = 0; pi < pathLengths[p]; pi++) {
				std::cout << paths[p][pi].edge << " " << paths[p][pi].index << ", ";
				}
				std::cout << std::endl;
			*/
			/*
				There exists a gap in this path.  There are a few options.
				- Remove this path from existence since it won't help resolving
				any equivalent transformations.
				- Keep the longest stretch of this path.  In some small
				cases this may help, but I don't anticipate that.

				for now,remove the entre path.
			*/

			std::vector<ssize_t> pathEdges, pathIndices;

			for (pi = 0; pi < pathLengths[p]; pi++) {
				if (paths[p][pi].edge != -1) {
					pathEdges.push_back(paths[p][pi].edge);
					pathIndices.push_back(paths[p][pi].index);
				}
				else {
					if (pathEdges.size() > 0) {
						/*
							std::cout << "removing subpath: " << std::endl;
							ssize_t pei;
							for (pei = 0; pei < pathEdges.size(); pei++) {
							std::cout << " " << pathEdges[pei] << " " << pathIndices[pei] << " ";
							}
							std::cout << std::endl;
						*/
						MarkPathForRemoval(pathEdges, pathIndices);
						pathEdges.clear();
						pathIndices.clear();
					}
				}
				paths[p][pi].edge = -1;
				paths[p][pi].index = -1;
			}
			if (pathEdges.size() > 0) {
				// process the last set of edges, if they exist
				/*
					std::cout << "removing subpath: " << std::endl;
					ssize_t pei;
					for (pei = 0; pei < pathEdges.size(); pei++) {
					std::cout << " " << pathEdges[pei] << " " << pathIndices[pei] << " ";
					}
					std::cout << std::endl;
				*/
				MarkPathForRemoval(pathEdges, pathIndices);
				pathEdges.clear();
				pathIndices.clear();
			}
			pathLengths[p] = 0;
		}
	}
	RemoveMarkedIntervals();
	//	RemoveLowCoverageEdges(2);
}

ssize_t IntervalGraph::DoesPathRepeat(ssize_t readIndex, ssize_t pathPos) {
	ssize_t p;
	assert(pathPos < pathLengths[readIndex]);
	for (p = 0; p < pathPos-1; p++) {
		if (paths[readIndex][p].edge == paths[readIndex][pathPos].edge)
			return 1;
	}
	return 0;
}


ssize_t IntervalGraph::RouteRemovedIntervals(ssize_t maxSearchLength) {
	//UNUSED// ssize_t p;
	ssize_t e, i ;
	ssize_t readIndex;
	ssize_t curPathPos;
	ssize_t numRouted = 0;
	//	std::cout << "Rerouting removed intervals"  << std::endl;
	UntraverseReadIntervals();
	Untraverse();
	std::vector<ssize_t> pathEdges, pathIndices;
	ssize_t compReadIndex;
	for (e = 0; e < edges.size(); e++ ) {
		for (i = 0; i < edges[e].intervals->size(); i++) {
			if (! IsIntervalMarkedForRemoval(e, i)) {
				// and
				//					edges[e].intervals[i].traversed == 0) {
				readIndex = (*edges[e].intervals)[i].read;
				/*
					std::cout << "searching for alternate paths for " << e << " " 
					<< edges[e].index << " " << readIndex << std::endl;
					PrintPath(readIndex, std::cout);
				*/
				ssize_t compEdge, compIndex;
				ssize_t compPathPos, compPathNext;

				// If this graph is balanced, we process the complimentary path
				// at the same time as the forward path.  If this path is 
				// complimentary, wait until the forward path is reached.
				if (isBalanced and !IsForwardRead(readIndex))
					continue;

				compReadIndex = GetCompReadIndex(readIndex);
				// Try and trace this path
				curPathPos = (*edges[e].intervals)[i].pathPos;
				if (curPathPos < pathLengths[readIndex]-1) {
					// This is not the last edge on the path, so it may have been gapped.
					ssize_t curPathIndex, curPathEdge;
					assert(readIndex >= 0);
					assert(curPathPos >= 0);
					curPathIndex = paths[readIndex][curPathPos].index;
					curPathEdge =  paths[readIndex][curPathPos].edge;

					assert(curPathEdge != -1);
					assert(curPathIndex != -1);
					if (paths[readIndex][curPathPos+1].edge == -1) {
						// The path is missing some edges.  If it's possible 
						// to re-join this path back to the graph, do so.
						//						std::cout << e << " " << i << " has been deleted! " << std::endl;
						ssize_t nextPathPos, nextPathIndex, nextPathEdge;
						for (nextPathPos = curPathPos+2; 
								 nextPathPos < pathLengths[readIndex]; 
								 nextPathPos++ ) {
							if (paths[readIndex][nextPathPos].edge != -1)
								/* 
									 Found the continuation of the path that has
									 not been deleted.  Break now and try and link
									 the path at the previous undeleted pos to next. 
								*/
								break;
						}
						// Only route the paths where the curPathPos and nextPathPos are 
						// within the range 0 .. pathLengths-1, so that part of the 
						// path is known.
						if (nextPathPos < pathLengths[readIndex]) {
							nextPathIndex = paths[readIndex][nextPathPos].index;
							nextPathEdge  = paths[readIndex][nextPathPos].edge;
							// std::cout << "the path resumes at " << nextPathEdge << " " << nextPathIndex << std::endl;
							std::vector<ssize_t> altPathEdges;


							if (curPathEdge != nextPathEdge) {
								//std::cout << "searching path for " << readIndex << std::endl;
								ssize_t foundAlternatePath = 0;
								ssize_t advance = 0;
								ssize_t gappedLength;
								ssize_t pathRepeats = 0;
								for (advance = 0; advance < 3 
											 and !foundAlternatePath
											 and nextPathPos < pathLengths[readIndex]
											 and paths[readIndex][nextPathPos].edge != -1; nextPathPos++, advance++ ){ 
									if (DoesPathRepeat(readIndex, nextPathPos)) {
										pathRepeats = 1;
										continue;
									}
									nextPathIndex = paths[readIndex][nextPathPos].index;
									nextPathEdge  = paths[readIndex][nextPathPos].edge;
									gappedLength = (*edges[nextPathEdge].intervals)[nextPathIndex].readPos - 
										((*edges[curPathEdge].intervals)[curPathIndex].readPos + 
										 (*edges[curPathEdge].intervals)[curPathIndex].length -
										 vertices[edges[curPathEdge].dest].vertexSize);
									// Search for an alternate path that is of similar length to the original
									// path.  
									if (StoreAlternatePath(curPathEdge, nextPathEdge, 
																				 altPathEdges, (ssize_t)( gappedLength * 1.20) )) {
										if (pathRepeats) {
											//std::cout << "found alternate path for cyclic path: " << readIndex << std::endl;
										}
										foundAlternatePath = 1;
										break;
									}
								}
								if (foundAlternatePath) {
									/*std::cout << "found alternate path for " << readIndex << std::endl;*/
									if (isBalanced and IsForwardRead(readIndex)) {
										compReadIndex = GetCompReadIndex(readIndex);
										compPathPos   = pathLengths[compReadIndex] - nextPathPos - 1;
										compPathNext  = pathLengths[compReadIndex] - curPathPos - 1;
										compEdge      = paths[compReadIndex][compPathPos].edge;
										compIndex     = paths[compReadIndex][compPathPos].index;
									}						
									else {
										// reset compliment indices since the graph is not complementary
										// or the complement has already been processed
										compPathPos = compPathNext = compEdge = compIndex = -1;
									}

									ssize_t toRemove;
									ssize_t toRemoveEdge, toRemoveInterval;
									for (toRemove = curPathPos + 1; toRemove < nextPathPos; toRemove++ ) {
										toRemoveEdge     = paths[readIndex][toRemove].edge;
										toRemoveInterval = paths[readIndex][toRemove].index;
										if (toRemoveEdge >= 0) {
											MarkIntervalForRemoval(toRemoveEdge, toRemoveInterval);
										}
									}
									for (toRemove = compPathPos + 1; toRemove < compPathNext; toRemove++) {
										toRemoveEdge = paths[compReadIndex][toRemove].edge;
										toRemoveInterval = paths[compReadIndex][toRemove].index;
										if (toRemoveEdge >= 0) {
											MarkIntervalForRemoval(toRemoveEdge, toRemoveInterval);
										}
									}

									if (curPathPos+1 == nextPathPos) {
										std::cout << "ERROR: the positions in the read path should not be " << std::endl
															<< " adjacent when replacing a gapped path " << std::endl;
										assert(0);
									}
									if (altPathEdges.size() > 0) {
										numRouted++;
										ssize_t newLength, ne;
										newLength = 0;
										for (ne = 0; ne < altPathEdges.size(); ne++) {
											newLength += edges[altPathEdges[ne]].length - 
												vertices[edges[altPathEdges[ne]].dest].vertexSize;
										}

										ReplacePathRangeForward(readIndex,
																						curPathEdge, curPathIndex,
																						curPathPos+1, nextPathPos, altPathEdges);
										assert(CheckPathContinuity(readIndex));

										assert(isBalanced and compReadIndex >= 0);
										std::vector<ssize_t> compPathEdges;
										FormComplimentPath(altPathEdges, compPathEdges);
										ReplacePathRangeReverse(compReadIndex, compEdge, compIndex,
																						compPathPos + 1, compPathNext, compPathEdges);
										assert(CheckPathContinuity(compReadIndex));
										if (CheckPathBalance(readIndex, compReadIndex)==0) {
											std::cout << "Tried to replace path range " << compPathPos +1 
																<< " " << compPathNext - 1 << " " << compPathEdges.size() <<" ";
											ssize_t c;
											for (c = 0; c < compPathEdges.size(); c++) { std::cout <<" " << compPathEdges[c];}
											std::cout << std::endl;
											PrintImbalancedPaths(readIndex, compReadIndex);
										}

										ssize_t ap;
										for (ap = 0; ap < altPathEdges.size(); ap++) {
											assert(edges[altPathEdges[ap]].intervals->size() ==
														 edges[compPathEdges[compPathEdges.size() - ap - 1]].intervals->size());
										}
									}
									else {
										SplicePathRange(readIndex, curPathPos+1, nextPathPos);
										assert(isBalanced and compReadIndex >= 0);
										SplicePathRange(compReadIndex, compPathPos+1, compPathNext);
									}
								}
							}
							else {
								// The deleted section loops back to itself, so just delete the intervening 
								// section
								SplicePathRange(readIndex, curPathPos+1, nextPathPos);

								compReadIndex = GetCompReadIndex(readIndex);
								compPathPos   = pathLengths[compReadIndex] - nextPathPos - 1;
								compPathNext  = pathLengths[compReadIndex] - curPathPos - 1;
								compEdge      = paths[compReadIndex][compPathPos].edge;
								compIndex     = paths[compReadIndex][compPathPos].index;

								assert(isBalanced and compReadIndex >= 0);
								SplicePathRange(compReadIndex, compPathPos, compPathNext);
							}
						} // end if nextPathPos < end
					} // end if curPathPos+1 is gapped 
				} // end if cur path pos is not at end
			} // end checking if edge is deleted
		}
	}
	//SortAllEdgeIntervalsByReadPos();
	//UpdateAllPathIndices();
	return numRouted;
}



void IntervalGraph::SplicePathRange(ssize_t readIndex, ssize_t spliceStart, ssize_t spliceEnd) {

	//
	// Remove spliceStart ... spliceEnd entirely from a path. 
	// This will effectively remove a read from the graph.  Unfortunately, for now
	// the sequence will be a bit messed up, since I doin't actually splice the read
	// I just shift the readPos offsets back.  In the future, I should
	// shift the readPos offsets back, and actually modify the corresponding read.
	assert (spliceEnd > spliceStart);
	//UNUSED// PathInterval *newPath;
	ssize_t splicedLength = spliceEnd - spliceStart;
	ssize_t newPathLength = pathLengths[readIndex] - splicedLength;
	ssize_t p;
	ssize_t edge, intv;
	// Calculate the length of the read that is deleted.
	// The read position is moved back for all subsequent reads.
	//UNUSED// ssize_t deletedLength = 0;
	if (spliceEnd == spliceStart) {
		// can't really splice anything here.
		return;
	}
	/*
		if (spliceStart == 0) {
		if (spliceEnd < pathLengths[readIndex]) {
			ssize_t edge, index;
			edge = paths[readIndex][spliceEnd].edge;
			index = paths[readIndex][spliceEnd].index;
			if ( edge >= 0 and index >= 0) {
				deletedLength = (*edges[edge].intervals)[index].readPos;
			}
		}
	}
	else {
		if (spliceEnd < pathLengths[readIndex]) {
			assert(spliceStart>0);
			ssize_t startEdge  = paths[readIndex][spliceStart-1].edge;
			ssize_t startIntv  = paths[readIndex][spliceStart-1].index;
			ssize_t destVertex = edges[startEdge].dest;
			ssize_t resumeEdge = paths[readIndex][spliceEnd].edge;
			ssize_t resumeIntv = paths[readIndex][spliceEnd].index;

			if (startEdge == -1 or resumeEdge == -1) {
				// Something is wrong with this path, remove the entire thing.
				for (p = 0; p < pathLengths[readIndex]; p++) {
					edge = paths[readIndex][p].edge;
					intv = paths[readIndex][p].index;
					if (edge >= 0 && index >= 0 and
							!IsIntervalMarkedForRemoval(edge, intv)) {
						MarkIntervalForRemoval(edge, intv);
					}
				}
				return;
			}

			assert(resumeEdge >= 0);
			assert(resumeIntv >= 0);
			assert(startEdge >= 0);
			assert(startIntv >= 0);
			assert(spliceEnd < pathLengths[readIndex]);
			deletedLength = (*edges[resumeEdge].intervals)[resumeIntv].readPos -
				((*edges[startEdge].intervals)[startIntv].readPos + 
				 (*edges[startEdge].intervals)[startIntv].length - 
				 vertices[destVertex].vertexSize);
		}
		// else, if spliceEnd == readLengths[readIndex]-1, the entire end
		// of the path is deleted, so no readPos offsets are necessary
	}
	*/
	// First remove all intervals corresponding to spliced region.
	for (p = spliceStart; p < spliceEnd; p++) {
		edge = paths[readIndex][p].edge;
		intv = paths[readIndex][p].index;
		if (edge != -1 and intv != -1 and
				!IsIntervalMarkedForRemoval(edge, intv)) {
			MarkIntervalForRemoval(edge, intv);
		}
	}

	for (p = spliceEnd ; p < pathLengths[readIndex]; p++) {
		assert(newPathLength > 0);
		paths[readIndex][p - splicedLength] = paths[readIndex][p];
		edge = paths[readIndex][p].edge;
		intv = paths[readIndex][p].index;
		if (edge >= 0 and intv >= 0 and !IsIntervalMarkedForRemoval(edge, intv)) {
			(*edges[edge].intervals)[intv].pathPos = p - splicedLength;
			//		(*edges[edge].intervals)[intv].readPos -= deletedLength;
		}
	}

	if (newPathLength > 0) {
		//    paths[readIndex] = newPath;
	}
	else {
		paths[readIndex] = NULL;
	}
	pathLengths[readIndex] = newPathLength;
}

void IntervalGraph::FormComplimentPath(std::vector<ssize_t> &origPath,
																			 std::vector<ssize_t> &compPath) {
	// Given orig path, find the path in the graph that
	// traversed all the complimentary edges
	compPath.resize(origPath.size());
	ssize_t e;
	for (e = 0; e < origPath.size(); e++) {
		compPath[origPath.size() - e - 1] = edges[origPath[e]].balancedEdge;
	}
}



ssize_t IntervalGraph::MakeRoomForEdges(ssize_t readIndex, ssize_t intvEdge, ssize_t intvIndex,
																		ssize_t pathStart, ssize_t pathEnd,
																		std::vector<ssize_t> &newEdges) {
	// This makes room for a range in a path.  The path from pathStart...pathEnd-1
	// gets replaced by newEdges.
	ssize_t newSpace;
	newSpace = newEdges.size() - (pathEnd - pathStart);
	//UNUSED// ssize_t e;
	ssize_t p ;
	//UNUSED// ssize_t curReadPos;
	PathInterval *newPath;
	if (newSpace != 0) {
		// Sanity check.  Make sure that some space is being allocated or deleted.
		if (pathLengths[readIndex] + newSpace == 0) {
			std::cout << "Allocating 0 space for a path, this shouldn't happen since we're adding" 
								<< std::endl;
			std::cout << "to paths. Check it out!" << std::endl;
			exit(0);
		}

		// Allocate the new path. This may be longer or shorter 
		// than the old path, depending on how much is being replaced.
		newPath = new PathInterval[pathLengths[readIndex] + newSpace];

		// Copy the old unmodified parts of the path into the newly allocated space.
		for (p = 0; p < pathStart; p++) {
			newPath[p] = paths[readIndex][p];
		}
		for (p = pathEnd; p < pathLengths[readIndex]; p++ ) {
			newPath[p + newSpace] = paths[readIndex][p];
			if (paths[readIndex][p].edge >= 0 && 
					paths[readIndex][p].index >= 0) 
				(*edges[paths[readIndex][p].edge].intervals)[paths[readIndex][p].index].pathPos += newSpace;
		}

		// The old path is not needed any more.
		delete[] paths[readIndex];
		paths[readIndex] = newPath;
	}
	return newSpace;
}

void IntervalGraph::CopyNewEdges(ssize_t readIndex, ssize_t intvEdge, ssize_t intvIndex,
																 ssize_t pathStart, std::vector<ssize_t> &newEdges) {

	// Copy edges from 'newEdges' into the path 'readIndex' starting at position
	// 'pathStart'. 
	// It is assumed that space has been made in the path for 'newEdges'.

	// Intervals are added to these edges, and initilalized as intervals that cover
	// the entire edge.  They might have to be fixed later.
	ssize_t e;
	ReadInterval newInterval;
	ssize_t curReadPos;
	curReadPos = (*edges[intvEdge].intervals)[intvIndex].readPos;
	//UNUSED// ssize_t lastInterval;
	ssize_t prevEdgeLength, prevVertexSize;
	ssize_t prevEdge, prevInterval;
	prevEdge       = intvEdge;
	if (pathStart > 0) {
		prevEdge       = paths[readIndex][pathStart-1].edge;
		prevInterval   = paths[readIndex][pathStart-1].index;
		prevEdgeLength = (*edges[prevEdge].intervals)[prevInterval].length;
		prevVertexSize = vertices[edges[prevEdge].dest].vertexSize;
	}
	else {
		prevEdge = -1;
		prevInterval = -1;
		prevEdgeLength = 0;
		prevVertexSize = 0;
	}

	for (e = 0; e < newEdges.size(); e++ ) {
		/* 
			 Record which edge this interval is in, and which pos the edge
			 references the interval.

			 Assume the interval uses the entire length of the edge.  This 
			 may be truncated later if it needs to be in order to conform
			 to the positions of the rest of the path.
		*/
		newInterval.read    = readIndex;
		//    newInterval.edge    = -1; // Don't store edges for now.
		newInterval.pathPos = pathStart + e;
		newInterval.length  = edges[newEdges[e]].length;
		newInterval.readPos = curReadPos + prevEdgeLength - prevVertexSize;

		/* This path passes through the beginning of the edge. */
		newInterval.edgePos = 0;

		/* Just make sure something strange isn't happening with the new edge.*/
		assert(newInterval.length > 0);

		/* Add this interval to the edge. */
		edges[newEdges[e]].intervals->push_back(newInterval);
		edges[newEdges[e]].multiplicity++;
		/* 
			 We've added the interval to the edge, so its index is 
			 simply the last one.
		*/

		/* Make the path reference it's spot in the edge.*/
		paths[readIndex][e + pathStart].edge  = newEdges[e];
		paths[readIndex][e + pathStart].index = edges[newEdges[e]].intervals->size()-1;

		/* Store values to use on the next iteration.*/
		prevInterval = edges[newEdges[e]].intervals->size() - 1;
		prevEdge     = newEdges[e];
		prevEdgeLength = (*edges[prevEdge].intervals)[prevInterval].length;
		prevVertexSize = vertices[edges[prevEdge].dest].vertexSize;
		curReadPos     = newInterval.readPos;
	}
}

ssize_t IntervalGraph::ReplacePathRangeForward(ssize_t readIndex,
																					 ssize_t intvEdge, ssize_t intvIndex,
																					 ssize_t pathStart, ssize_t pathEnd,
																					 std::vector<ssize_t> &newEdges) {

	/*
		Given a path 'readIndex', replace a subpath with a different path.  
		The path from (pathStart ... pathEnd] is repalced by the path contained 
		in 'newEdges'.
		The rest of the path is filled in.
		If pathStart == pathEnd, this effectively inserts the path 
		before pathStart.

		After replacing the subpath, the lengths of the new path are updated
		so that the end position of the new subpath is 1- the beginning position
		of succeeding path interval.

		Since the path intervals are given, I don't think it's necessary to provide
		intvEdge and intvIndex as parameters, but that needs to be examined later.

		This can be tested with an assert for:
		paths[pathStart].edge == intvEdge  and
		paths[pathStart].index == intvIndex
	*/
	ssize_t origPathLength = pathLengths[readIndex];

	ReplacePathRange(readIndex, intvEdge, intvIndex, pathStart, pathEnd, newEdges);

	/*
		Now update the read and path positions.
	*/
	/*
		Update the values to be used on the next iteration.
	*/
	/* If this is the last new edge, offset the length so that 
		 the path reaches the next edge on the path.
	*/
	ssize_t lastNewEdge    = paths[readIndex][pathStart + newEdges.size()-1].edge;
	ssize_t lastNewIntv    = paths[readIndex][pathStart + newEdges.size()-1].index;
	ssize_t lastNewIntvLength = (*edges[lastNewEdge].intervals)[lastNewIntv].length;
	ssize_t lastVertexSize = vertices[edges[lastNewEdge].dest].vertexSize;
	ssize_t lastReadPos  = (*edges[lastNewEdge].intervals)[lastNewIntv].readPos;

	ssize_t nextEdge;
	ssize_t nextEdgeInterval;

	if (origPathLength < pathEnd-1) {
		nextEdge         = paths[readIndex][pathStart + newEdges.size()].edge;
		nextEdgeInterval = paths[readIndex][pathStart + newEdges.size()].index;
	}
	else {
		nextEdge = nextEdgeInterval = -1;
		// The entire end of the path was replaced, no need to modify the path length.
		return 1;
	}


	// IF the replaced path joins perfectly with the old path,
	// no interval length adjustments are necessary.
	if (nextEdge != -1 and lastReadPos + lastNewIntvLength - lastVertexSize == 
			(*edges[nextEdge].intervals)[nextEdgeInterval].readPos) {
		return 1;
	}


	// Otherwise, it is necessary to adjust the length of the intervals.


	ssize_t oldPathLength = (*edges[nextEdge].intervals)[nextEdgeInterval].readPos -
		((*edges[intvEdge].intervals)[intvIndex].readPos +
		 (*edges[intvEdge].intervals)[intvIndex].length - 
		 vertices[edges[intvEdge].dest].vertexSize);

	ssize_t newPathLength = 0;
	ssize_t e;
	for (e = 0; e < newEdges.size(); e++ ) {
		newPathLength += edges[newEdges[e]].length;
		newPathLength -= vertices[edges[newEdges[e]].dest].vertexSize;
	}

	ssize_t readPos;
	readPos      = (*edges[intvEdge].intervals)[intvIndex].readPos +
		(*edges[intvEdge].intervals)[intvIndex].length - 
		vertices[edges[intvEdge].dest].vertexSize;


	// now fix the lengths of the intervals
	ssize_t edge=-1, intv=-1;
	//UNUSED// ssize_t  edgeLength;
	ssize_t intervalLength;
	double remainder, exactIntervalLength;
	remainder = 0;

	//UNUSED//	int endReadPos = (*edges[intvEdge].intervals)[intvIndex].readPos;
	for (e = 0; e < newEdges.size() ; e++) {
		edge = paths[readIndex][pathStart + e].edge;
		intv = paths[readIndex][pathStart + e].index;

		exactIntervalLength =
			((double((*edges[edge].intervals)[intv].length
							 - vertices[edges[edge].dest].vertexSize)
				/ newPathLength) * oldPathLength
			 + remainder);


		intervalLength = (ssize_t) floor(exactIntervalLength);

		// Account for accruing round-off error
		remainder = exactIntervalLength - intervalLength;

		// Fix the edge pos to 
		(*edges[edge].intervals)[intv].readPos = readPos;
		assert(readPos >= 0);
		(*edges[edge].intervals)[intv].readPos = readPos;
		(*edges[edge].intervals)[intv].length = intervalLength + vertices[edges[edge].src].vertexSize;
		readPos  += intervalLength;
		assert(readPos >= 0);
	}

	// the hope is that the read pos ends right where the original 
	// path continued.
	if (readPos != (*edges[nextEdge].intervals)[nextEdgeInterval].readPos) {
		(*edges[edge].intervals)[intv].length -= 
			(readPos - (*edges[nextEdge].intervals)[nextEdgeInterval].readPos);

	}
	return 1;
}

ssize_t IntervalGraph::ReplacePathRangeReverse(ssize_t readIndex,
																					 ssize_t intvEdge, ssize_t intvIndex,
																					 ssize_t pathStart, ssize_t pathEnd,
																					 std::vector<ssize_t> &newEdges) {
  assert(newEdges.size() > 0);

	ReplacePathRange(readIndex, intvEdge, intvIndex, pathStart, pathEnd, newEdges);

	/*
		Now update the read and path positions.
	*/
	/* Update the values to be used on the next iteration.*/
	/* If this is the last new edge, offset the length so that 
		 the path reaches the next edge on the path.
	*/
	ssize_t lastNewIntv     = paths[readIndex][pathStart+newEdges.size()-1].index;
	ssize_t lastNewEdge     = paths[readIndex][pathStart+newEdges.size()-1].edge;

	ssize_t nextEdge         = paths[readIndex][pathStart + newEdges.size()].edge;
	ssize_t nextEdgeInterval = paths[readIndex][pathStart + newEdges.size()].index;

	// IF the replace path joins perfectly with the old path,
	// no interval length adjustments are necessary.
	// The beginning of the spliced in section of the path is guaranteed
	// to fit with the end of the last part of the old path 
	// before the spliced in section.
	if ((*edges[lastNewEdge].intervals)[lastNewIntv].readPos + 
			(*edges[lastNewEdge].intervals)[lastNewIntv].length -
			vertices[edges[lastNewEdge].dest].vertexSize == 
			(*edges[nextEdge].intervals)[nextEdgeInterval].readPos) {
		return 1;
	}

	// Otherwise, it is necessary to adjust the length of the intervals.
	ssize_t oldPathLength = (*edges[nextEdge].intervals)[nextEdgeInterval].readPos -
		((*edges[intvEdge].intervals)[intvIndex].readPos + 
		 (*edges[intvEdge].intervals)[intvIndex].length - 
		 vertices[edges[intvEdge].dest].vertexSize);

	ssize_t newPathLength = 0;
	ssize_t e;
	for (e = 0; e < newEdges.size(); e++ ) {
		newPathLength += edges[newEdges[e]].length;
		newPathLength -= vertices[edges[newEdges[e]].dest].vertexSize;
	}
	ssize_t readPos;
	readPos      = (*edges[nextEdge].intervals)[nextEdgeInterval].readPos;

	// now fix the lengths of the intervals, starting at 
	// the END of the new path, so that the length used
	// corresponds to the beginning of the old path
	ssize_t edge=-1, intv=-1; // Quiet compiler warnings.  Since newEdges.size()>0, these will be replaced.
	//UNUSED// ssize_t  edgeLength;
	ssize_t intervalLength;
	double remainder, exactIntervalLength;
	remainder = 0;
	//UNUSED//	int endReadPos = (*edges[intvEdge].intervals)[intvIndex].readPos;
	for (e = newEdges.size()-1; e >= 0 ; e--) {
		edge = paths[readIndex][pathStart + e].edge;
		intv = paths[readIndex][pathStart + e].index;

		exactIntervalLength = 
			(double( (*edges[edge].intervals)[intv].length
							 - vertices[edges[edge].src].vertexSize
							 )
			 / newPathLength) * oldPathLength
			+ remainder;

		intervalLength = (ssize_t) floor(exactIntervalLength);

		// Account for accruing round-off error
		remainder = exactIntervalLength - intervalLength;
		readPos  -= intervalLength;

		// Fix the edge pos to 
		(*edges[edge].intervals)[intv].readPos = readPos;
		(*edges[edge].intervals)[intv].length  = intervalLength + vertices[edges[edge].src].vertexSize;
	}
	// back up the last read interval
	readPos -= (*edges[intvEdge].intervals)[intvIndex].length - vertices[edges[intvEdge].dest].vertexSize;
	// the hope is that the read pos ends right where the original 
	// path continued.
	if (readPos != (*edges[intvEdge].intervals)[intvIndex].readPos) {
		ssize_t diff = (*edges[intvEdge].intervals)[intvIndex].readPos- readPos;
		(*edges[edge].intervals)[intv].length -= diff;
		(*edges[edge].intervals)[intv].readPos += diff;
	}
	return 1;
}

ssize_t IntervalGraph::ReplacePathRange(ssize_t readIndex, 
																		ssize_t intvEdge, ssize_t intvIndex,
																		ssize_t pathStart, ssize_t pathEnd,
																		std::vector<ssize_t> &newEdges) {
	/* 
		 Replace part of a path starting at pathStart and ending at pathEnd
		 with the edges listed in newEdges.
		 Also, add the path intervals to the all edges in newEdges.
		 For now, assume that the path passes through the entire edge in newEdges.
	*/

	/*
		First fix the size of the path if need be.
	*/

	ssize_t newSpace;
	/* 
		 newSpace == 0, no new path needed.
		 < 0, must add space to the old path
		 > 0, must delete from the old path.
	*/
	newSpace = MakeRoomForEdges(readIndex, intvEdge, intvIndex,
															pathStart, pathEnd, newEdges);

	CopyNewEdges(readIndex, intvEdge, intvIndex, pathStart, newEdges);

	pathLengths[readIndex] += newSpace;

	return 1;
}


ssize_t IntervalGraph::StoreAlternatePath(ssize_t curPathEdge, ssize_t nextPathEdge, 
																			std::vector<ssize_t> &altPathEdges,
																			ssize_t maxSearchLength) {
	return 0;
	std::list<ssize_t> altPathEdgeList;
	if (SearchAlternatePath(curPathEdge, nextPathEdge, altPathEdgeList, maxSearchLength)) {
		std::list<ssize_t>::iterator listIt;
		altPathEdges.resize(altPathEdgeList.size());
		ssize_t i = 0;
		for (listIt = altPathEdgeList.begin();
				 listIt != altPathEdgeList.end();
				 listIt++) {
			altPathEdges[i] = *listIt;
			++i;
		}
		return 1;
	}
	else
		return 0;
}

ssize_t IntervalGraph::SearchAlternatePath(ssize_t curPathEdge, ssize_t nextPathEdge, 
																			 std::list<ssize_t> &altPathEdges,
																			 ssize_t maxSearchLength,
																			 ssize_t curPathLength) {
	/* Input:
		 curPathEdge, curPathIndex - the edge and index of a read interval on that edge.
		 nextPathEdge, nextPathIndex - the continuation of that path on some edge that
		 may or may not be adjacent to curPathEdge.
	*/


	if (curPathEdge == nextPathEdge) {
		return 1;
	}
	else {
		ssize_t outEdge, outEdgeIndex;
		ssize_t destVertex;
		destVertex = edges[curPathEdge].dest;
		for (outEdgeIndex = vertices[destVertex].FirstOut();
				 outEdgeIndex < vertices[destVertex].EndOut();
				 outEdgeIndex = vertices[destVertex].NextOut(outEdgeIndex)) {
			outEdge = vertices[destVertex].out[outEdgeIndex];

			// No long cycles.
			if (edges[curPathEdge].traversed == GraphEdge::Marked)
				continue;

			// No very short cycles.
			if (vertices[edges[curPathEdge].dest].traversed == GraphVertex::Marked) // TODO: check; was = instead of ==
				continue;

			if (edges[curPathEdge].dest == edges[nextPathEdge].src and // make sure this is reached from the right direction
					outEdge == nextPathEdge) {
				return 1;
			}
			else {
				/* The out edge is not the continuation of the path
					 we are looking for.  Continue with dfs for nextPathEdge.
				*/

				ssize_t nextEdgeLength = edges[outEdge].length - vertices[edges[outEdge].dest].vertexSize;
				/*
					std::cout << "sap " << altPathEdges.size() << " cur: " << curPathLength << " next " 
					<< nextEdgeLength << " max: " << maxSearchLength <<  " " 
					<< edges[outEdge].index << " " << std::endl;
				*/
				if (maxSearchLength > curPathLength  + nextEdgeLength) {
					altPathEdges.push_back(outEdge);
					edges[curPathEdge].traversed = GraphEdge::Marked;
					vertices[edges[curPathEdge].dest].traversed = GraphVertex::Marked;
					if (SearchAlternatePath(outEdge, nextPathEdge, altPathEdges, 
																	maxSearchLength, 
																	curPathLength +  nextEdgeLength)) {
						edges[curPathEdge].traversed = GraphEdge::NotMarked;
						return 1;
					}
					else {
						/* The search didn't find the nextPathEdge within maxSearchLength.
							 Don't use that edge.
						*/
						altPathEdges.pop_back();
						edges[curPathEdge].traversed = GraphEdge::NotMarked;
					}
				}
			}
		}
		/* 
			 All out edges have been tried, and no path works.
			 Return 0 for not found.
		*/
		return 0;
	}
}

void IntervalGraph::UntraverseReadIntervals() {
	ssize_t e, i;
	for (e = 0; e < edges.size(); e++) {
		for (i = 0; i < edges[e].intervals->size(); i++ ) {
			(*edges[e].intervals)[i].traversed = 0;
		}
	}
}

void IntervalGraph::DeleteEdgeListReadIntervals(std::vector<ssize_t> &edgeList) {
	ssize_t edge;
	for (edge = 0; edge < edgeList.size(); edge++) { 
		DeleteEdgeReadIntervals(edgeList[edge]);
	}
}

void IntervalGraph::DeleteEdgeReadIntervals(ssize_t edge) {
	ssize_t i;
	for (i = 0; i < edges[edge].intervals->size(); i++ ) {
		DeleteReadInterval((*edges[edge].intervals)[i].read,
											 (*edges[edge].intervals)[i].pathPos);
	}
}

void IntervalGraph::DeleteReadInterval(ssize_t readIndex, ssize_t pathPos) {
	PathInterval *newPath;
	/*
		std::cout << "deleting read interval " << readIndex << " " 
		<< pathPos << " " << pathLengths[readIndex] << std::endl;
	*/
	assert(pathLengths[readIndex] > 0);
	if (pathLengths[readIndex]-1 > 0) {
		newPath = new PathInterval[pathLengths[readIndex]-1];
		if (pathPos > 0)
			std::copy(&(paths[readIndex][0]), 
								&(paths[readIndex][pathPos]), 
								&(newPath[0]));
		if (pathPos < pathLengths[readIndex]-1) {
			std::copy(&(paths[readIndex][pathPos+1]), 
								&(paths[readIndex][pathLengths[readIndex]]),
								&(newPath[pathPos]));
		}
		/* Update the references from the edge back to the path. */
		ssize_t p;
		for (p = pathPos+1; p < pathLengths[readIndex]; p++) {
			if (!IsIntervalMarkedForRemoval(paths[readIndex][p].edge, paths[readIndex][p].index)) {
				(*edges[paths[readIndex][p].edge].intervals)[paths[readIndex][p].index].pathPos--;
				assert((*edges[paths[readIndex][p].edge].intervals)[paths[readIndex][p].index].pathPos >= 0);
			}
		}
		delete[] paths[readIndex];
		paths[readIndex] = newPath;
	}
	else {
		delete[] paths[readIndex];
		paths[readIndex] = NULL;
	}
	pathLengths[readIndex]--;
}

void IntervalGraph::RerouteSimplePathIntervals(ssize_t vertex, ssize_t inEdge, ssize_t outEdge) {

	assert(vertices[vertex].InDegree() == 1);

	// First extend all read intervals that pass through the in edge.
	// Read intervals are sorted by read than position on the read.

	ssize_t inInterval, outInterval;
	inInterval = 0;
	outInterval = 0;
	ssize_t numInIntervals  = edges[inEdge].intervals->size();
	ssize_t numOutIntervals = edges[outEdge].intervals->size();
	ssize_t numMerged = 0;
	ssize_t numOutIntervalsMarkedForRemoval = 0;

	ssize_t inPathPos;
	ssize_t readIndex;
	ssize_t nextEdge, nextEdgeIntv;
	for (inInterval = 0; inInterval < numInIntervals; inInterval++ ) {
		if (!IsIntervalMarkedForRemoval(inEdge, inInterval)) {

			inPathPos = (*edges[inEdge].intervals)[inInterval].pathPos;
			readIndex = (*edges[inEdge].intervals)[inInterval].read;
			if (inPathPos < pathLengths[readIndex]-1) {

				nextEdge     = paths[readIndex][inPathPos+1].edge;
				nextEdgeIntv = paths[readIndex][inPathPos+1].index;

				if (!IsIntervalMarkedForRemoval(nextEdge, nextEdgeIntv)) {
					assert(nextEdge == outEdge);
					assert((*edges[inEdge].intervals)[inInterval].readPos + 
								 (*edges[inEdge].intervals)[inInterval].length - 
								 vertices[edges[inEdge].dest].vertexSize == 
								 (*edges[nextEdge].intervals)[nextEdgeIntv].readPos);

					(*edges[inEdge].intervals)[inInterval].length += 
						(*edges[nextEdge].intervals)[nextEdgeIntv].length - 
						vertices[edges[inEdge].dest].vertexSize;
					assert((*edges[inEdge].intervals)[inInterval].length >= 0);
					++numMerged;

					// This out interval has been merged into the in interval.
					// We could remove it from the list, but that's too much
					// memory management.  Just simply mark the start pos of the 
					// out interval as 0
					ssize_t read, pathPos;

					// Flag this interval as being merged
					(*edges[nextEdge].intervals)[nextEdgeIntv].readPos = -1;

					read = (*edges[nextEdge].intervals)[nextEdgeIntv].read;
					pathPos = (*edges[nextEdge].intervals)[nextEdgeIntv].pathPos;

					MarkIntervalForRemoval(nextEdge, nextEdgeIntv);
					DeleteReadInterval(read, pathPos);
				}
			}
		}
	}

	// The intervals that remain in 'outEdge' are not part of 'inEdge'
	// They should be added to 'inEdge', with a starting position
	// to offset for the length of 'inEdge'.

	// Count the number of out intervals that will be skipped.
	ssize_t remainingOutIntervals = 0;
	for( outInterval = 0; outInterval < numOutIntervals; ++outInterval) {
		// if read pos < 0, this interval has already been merged
		if (IsIntervalMarkedForRemoval(outEdge, outInterval)) {
			numOutIntervalsMarkedForRemoval++;
		}
		else 
			remainingOutIntervals++;
	}

	assert(numOutIntervals == numOutIntervalsMarkedForRemoval + remainingOutIntervals);

	if (numInIntervals + numOutIntervals == 0) {
		edges[inEdge].intervals->clear();
	}
	else {
		edges[inEdge].intervals->resize(numInIntervals + remainingOutIntervals);
		outInterval = 0;
		ssize_t appendedIntervals = 0;
		ssize_t outIntervalRead, outIntervalPos;
		for( outInterval = 0; outInterval < numOutIntervals; ++outInterval) {
			// if read pos < 0, this interval has already been merged
			if ((*edges[outEdge].intervals)[outInterval].readPos >= 0) {
				if (!IsIntervalMarkedForRemoval(outEdge, outInterval)) {
					(*edges[inEdge].intervals)[numInIntervals + appendedIntervals].read =
						(*edges[outEdge].intervals)[outInterval].read;

					(*edges[inEdge].intervals)[numInIntervals + appendedIntervals].readPos =
						(*edges[outEdge].intervals)[outInterval].readPos;

					(*edges[inEdge].intervals)[numInIntervals + appendedIntervals].length =
						(*edges[outEdge].intervals)[outInterval].length;

					(*edges[inEdge].intervals)[numInIntervals + appendedIntervals].edgePos =
						edges[inEdge].length + (*edges[outEdge].intervals)[outInterval].edgePos - vertexSize;

					(*edges[inEdge].intervals)[numInIntervals + appendedIntervals].pathPos =
						(*edges[outEdge].intervals)[outInterval].pathPos;
					assert((*edges[inEdge].intervals)[numInIntervals + appendedIntervals].pathPos >= 0);

					/* 
						 Update the path to reference inEdge, and the new slot in inEdge.
					*/
					outIntervalRead = (*edges[outEdge].intervals)[outInterval].read;
					outIntervalPos  = (*edges[outEdge].intervals)[outInterval].pathPos;
					paths[outIntervalRead][outIntervalPos].edge  = inEdge;
					paths[outIntervalRead][outIntervalPos].index = numInIntervals + appendedIntervals;

					++appendedIntervals;
				}
				else {
					std::cout << "Caution!!! interval: " << outEdge << " " << outInterval 
										<< " is marked for removal, and I should have removed everything like that " 
										<< std::endl;
				}
			}
		}
		// sanity check
		//		assert(appendedIntervals == numOutIntervals - numMerged - numOutIntervalsMarkedForRemoval);
		assert(appendedIntervals == remainingOutIntervals);
	}
}

void IntervalGraph::UpdatePathEdges(std::vector<ssize_t> &edgesToRemove) {
	/* 
		 Input: a list of edges to remove.
		 Result: paths that reference an edge that is not removed
		 have that edge packed, so if there are edges 1, 2, and 3, 
		 and edge 2 is deleted, paths that reference edge 3 are updated
		 to reference 2.
	*/

	// Nothing needs to be done if nothing is removed.
	if (edgesToRemove.size() == 0) 
		return;

	ssize_t edge;
	ssize_t numRemovedEdges = 0;
	for (edge = 0; edge < edges.size(); edge++) {
		if (numRemovedEdges < edgesToRemove.size() and
				edge == edgesToRemove[numRemovedEdges]) {
			++numRemovedEdges;
		}
		else {
			ssize_t i;
			ssize_t read;
			ssize_t pathPos;
			for (i = 0; i < edges[edge].intervals->size(); i++ ) {
				if (!IsIntervalMarkedForRemoval(edge,i)) {
					read = (*edges[edge].intervals)[i].read;
					pathPos = (*edges[edge].intervals)[i].pathPos;
					assert(read >= 0);
					assert(read < paths.size());
					assert(pathPos < pathLengths[read]);
					assert(pathPos >= 0);
					assert(paths[read][pathPos].edge >= 0);
					assert(paths[read][pathPos].index >= 0);
					assert(paths[read][pathPos].edge == edge);
					assert(paths[read][pathPos].edge - numRemovedEdges >= 0);
					paths[read][pathPos].edge -= numRemovedEdges;
				}
			}
		}
	}
}

void IntervalGraph::Prune(std::vector<ssize_t> &verticesToRemove,
													std::vector<ssize_t> &edgesToRemove) {

	std::cout << CurTimeString() << ": Prune(verticesToRemove[list of " << verticesToRemove.size() << "] , edgesToRemove[list of " << edgesToRemove.size() << "])" << std::endl;

	std::sort(edgesToRemove.begin(), edgesToRemove.end());
	std::sort(verticesToRemove.begin(), verticesToRemove.end());

	// Do a sanity check on the edges to remove. If any single edge
	// is in the to remove list, then its balanced edge should be in as well
	ssize_t re;
	ssize_t balancedEdge;
	std::vector<ssize_t>::iterator reIt;
	for (re = 0; re < edgesToRemove.size(); re++ ){
		balancedEdge = edges[edgesToRemove[re]].balancedEdge;
		// We don't want to do a sanity check on edges that
		// do not have balanced parts.

		if (balancedEdge >= 0 and
				std::binary_search(edgesToRemove.begin(), edgesToRemove.end(), balancedEdge) == 0) {
			std::cout << "ERROR! Removing edge "<< edgesToRemove[re] 
								<< " (" << edges[edgesToRemove[re]].index << ") "
								<< " but not balanced edge " 
								<< balancedEdge << std::endl;
		}
	}
	if (verticesToRemove.size() > 0) {
		/*
			std::cout << "Deleting " << verticesToRemove.size() 
			<< " vertices out of " << vertices.size()
			<< std::endl;
		*/
		RemoveVertexList(verticesToRemove);
		//    DeleteElementsFromList(vertices, verticesToRemove);
	}
	if (edgesToRemove.size() > 0) {
		RemoveEdgeList(edgesToRemove);
	}
	//	 RemoveTruncatedPathIntervals();
	assert(CheckEdges(vertices, edges));
	//  assert(CheckBalance());
	std::cout << CurTimeString() << ": Exit Prune()" << std::endl;
}

void IntervalGraph::RemoveVertexList(std::vector<ssize_t> &verticesToRemove) {
	std::cout << CurTimeString() << ": RemoveVertexList(verticesToRemove[list of " << verticesToRemove.size() << "])" << std::endl;

	ssize_t v;
	ssize_t newV, removedV;

	newV = 0;
	removedV = 0;
	ssize_t inEdgeIndex, inEdge, outEdgeIndex, outEdge;
	ssize_t removedVertices = 0;
	for (v = 0; v < vertices.size(); v++) {
		if (removedV < verticesToRemove.size() and
				v == verticesToRemove[removedV]) {
			while (v < vertices.size() and 
						 removedV < verticesToRemove.size() and 
						 v == verticesToRemove[removedV]) {
				assert(removedV < verticesToRemove.size());
				vertices[v].out.clear();
				vertices[v].in.clear();
				++removedV;
			}
			++removedVertices;
		}
		else {
			vertices[newV] = vertices[v];
			for (inEdgeIndex = vertices[newV].FirstIn();
					 inEdgeIndex != vertices[newV].EndIn();
					 inEdgeIndex = vertices[newV].NextIn(inEdgeIndex)) {
				inEdge = vertices[newV].in[inEdgeIndex];
				edges[inEdge].dest = newV;
			}
			for (outEdgeIndex = vertices[newV].FirstOut();
					 outEdgeIndex != vertices[newV].EndOut();
					 outEdgeIndex = vertices[newV].NextOut(outEdgeIndex)) {
				outEdge = vertices[newV].out[outEdgeIndex];
				edges[outEdge].src = newV;
			}
			newV++;
		}
	}
	if (removedVertices < vertices.size())
		vertices.resize(vertices.size() - removedVertices);
	else
		vertices.clear();
	std::cout << CurTimeString() << ": Exit RemoveVertexList()" << std::endl;
}


void IntervalGraph::RemoveUnlinkedEdges() {
	vector<ssize_t> vlist, elist;
	ssize_t e;
	for (e = 0; e < edges.size(); e++ ){
		if (edges[e].src == -1) {
			assert(edges[e].dest == -1);
			ssize_t i;
			for (i = 0; i < (*edges[e].intervals).size(); i++) {
				if (!(*edges[e].intervals)[i].markedForDeletion) {
					//					cout << (*edges[e].intervals)[i].read << " in unlinked edge. ";
						/*
					ssize_t pi;
					ssize_t path = (*edges[e].intervals)[i].read;
				for (pi = 0; pi < pathLengths[path]; pi++ ){
						cout << paths[path][pi].edge << " ";
					}
					cout << endl;
					*/
					RemovePath((*edges[e].intervals)[i].read);
				}
			}
			elist.push_back(e);
		}
	}
	RemoveMarkedIntervalsNoPaths();
	Prune(vlist, elist);
}


// TODO: This currently loops over all edges of the graph (huge list),
// checking them against edgesToRemove (potentially a much smaller list).
// Instead, change it to loop over edgesToRemove, from end to start,
// deleting the edges and moving the then-final edge from the graph
// into the slot of the deleted edge.
void IntervalGraph::RemoveEdgeList(std::vector<ssize_t> &edgesToRemove) {
	ssize_t e, newE, removedE;
	newE = 0;
	removedE = 0;
	ssize_t src, srcOutIndex, dest, destInIndex;
	ssize_t i;
	ssize_t removedEdges = 0;
	for (e = 0; e < edges.size(); e++ ) {
		if (removedE < edgesToRemove.size() and
				e == edgesToRemove[removedE]) {
			++removedEdges;
			while (e < edges.size() and 
						 removedE < edgesToRemove.size() and
						 e == edgesToRemove[removedE]) {
				edges[e].Clear();
				removedE++;
			}
		}
		else {
			if (newE != e) {

				edges[newE] = edges[e];
				src = edges[newE].src;

				srcOutIndex = vertices[src].LookupOutIndex(e);
				assert(srcOutIndex >= 0);
				vertices[src].out[srcOutIndex] = newE;

				dest = edges[newE].dest;
				destInIndex = vertices[dest].LookupInIndex(e);
				assert(destInIndex >= 0);
				vertices[dest].in[destInIndex] = newE;

				ssize_t balE;
				balE = edges[newE].balancedEdge;
				if (balE != e) {
					edges[balE].balancedEdge = newE;
				}
				else {
					edges[newE].balancedEdge = newE;
				}

				// update the paths to index the new edge index
				ReadIntervalList *intvList = edges[newE].intervals;
				for (i= 0; i < intvList->size(); i++) { 
					ReadInterval *intvPtr = &((*intvList)[i]);
					if (! intvPtr->markedForDeletion) {
						assert(intvPtr->read >= 0);
						assert(intvPtr->pathPos >= 0);
						paths[intvPtr->read][intvPtr->pathPos].edge = newE;
					}
				}
				edges[newE].altEdges = edges[e].altEdges;
				edges[e].altEdges.clear();
			}
			++newE;
		}
	}
	if (removedEdges < edges.size()) {
		edges.resize(edges.size() - removedEdges);
	}
	else
		edges.clear();
}

void ErodeLeavesFunctor::operator()(ssize_t vertexIndex) {
	ssize_t inDegree = (*vertices)[vertexIndex].InDegree();
	ssize_t outDegree = (*vertices)[vertexIndex].OutDegree();
	ssize_t destVertex, srcVertex;
	//UNUSED// ssize_t e;
	ssize_t balancedIndex;
	//UNUSED// ssize_t balancedSrcVertex;
	ssize_t balancedDestVertex ;
	//UNUSED// ssize_t i;
	if ( outDegree == 0 and inDegree > 0) {
		ssize_t inEdge, inEdgeIndex;
		inEdge = (*vertices)[vertexIndex].FirstIn();
		ssize_t allShortSources = 1;
		while (inEdge < (*vertices)[vertexIndex].EndIn() and
					 allShortSources) {
			inEdgeIndex = (*vertices)[vertexIndex].in[inEdge];

			if (!(*edges)[inEdgeIndex].IsShort(minLength) or 
					(*vertices)[(*edges)[inEdgeIndex].src].marked == GraphVertex::Marked) {
				allShortSources = 0;
			}
			inEdge = (*vertices)[vertexIndex].NextIn(inEdge);
		}    
		if (allShortSources) {
			inEdge = (*vertices)[vertexIndex].FirstIn();
			while (inEdge < (*vertices)[vertexIndex].EndIn()) {
				inEdgeIndex = (*vertices)[vertexIndex].in[inEdge];
				srcVertex = (*edges)[inEdgeIndex].src;
				(*vertices)[srcVertex].marked = GraphVertex::Marked;
				if ((*edges)[inEdgeIndex].IsShort(minLength) or 
						(*edges)[inEdgeIndex].flagged == GraphEdge::Marked) {
					(*edges)[inEdgeIndex].flagged = GraphEdge::Marked;
					erodedEdgeList.push_back(inEdgeIndex);
					//					(*vertices)[srcVertex].EraseOutEdge(inEdgeIndex);
					if ((*vertices)[vertexIndex].flagged != GraphVertex::Marked)
						erodedVertexList.push_back(vertexIndex);

					(*vertices)[vertexIndex].flagged = GraphVertex::Marked;


					if (graph->isBalanced) {
						ssize_t balEdge = (*edges)[inEdgeIndex].balancedEdge;
						ssize_t balEdgeSrc = (*edges)[balEdge].src;
						ssize_t balEdgeDest = (*edges)[balEdge].dest;
						if ((*vertices)[balEdgeSrc].flagged != GraphVertex::Marked)
							erodedVertexList.push_back(balEdgeSrc);
						(*vertices)[balEdgeSrc].flagged = GraphVertex::Marked;
						// Only delete one edge from the dest per iteration.
						(*vertices)[balEdgeDest].marked = GraphVertex::Marked;

						(*edges)[balEdge].flagged = GraphEdge::Marked;
						erodedEdgeList.push_back(balEdge);

					}
				}
				inEdge = (*vertices)[vertexIndex].NextIn(inEdge);
			}
		}
	}

	if ( inDegree == 0 and outDegree > 0) {
		ssize_t outEdge, outEdgeIndex;
		ssize_t allShortSinks = 1;
		outEdge = (*vertices)[vertexIndex].FirstOut();
		while (outEdge < (*vertices)[vertexIndex].EndOut() and
					 allShortSinks ) {
			outEdgeIndex = (*vertices)[vertexIndex].out[outEdge];
			if (!(*edges)[outEdgeIndex].IsShort(minLength) or 
					(*vertices)[(*edges)[outEdgeIndex].dest].marked == GraphVertex::Marked) {
				allShortSinks = 0;
			}
			outEdge = (*vertices)[vertexIndex].NextOut(outEdge);
		}
		if (allShortSinks) {
			outEdge = (*vertices)[vertexIndex].FirstOut();
			while (outEdge < (*vertices)[vertexIndex].EndOut()) {
				outEdgeIndex = (*vertices)[vertexIndex].out[outEdge];
				destVertex = (*edges)[outEdgeIndex].dest;
				(*vertices)[destVertex].marked = GraphVertex::Marked;
				if ((*edges)[outEdgeIndex].IsShort(minLength) or
						(*edges)[outEdgeIndex].flagged == GraphEdge::Marked) {
					// Signal these edges and vertices as GONE!
					(*edges)[outEdgeIndex].flagged = GraphEdge::Marked;
					erodedEdgeList.push_back(outEdgeIndex);
					//					(*vertices)[destVertex].EraseInEdge(outEdgeIndex);
					if ((*vertices)[vertexIndex].flagged != GraphVertex::Marked) 
						erodedVertexList.push_back(vertexIndex);
					(*vertices)[vertexIndex].flagged = GraphVertex::Marked;

					if (graph->isBalanced) {
						ssize_t balEdge = (*edges)[outEdgeIndex].balancedEdge;
						ssize_t balEdgeDest = (*edges)[balEdge].dest;
						ssize_t balEdgeSrc  = (*edges)[balEdge].src;
						if ((*vertices)[balEdgeDest].flagged != GraphVertex::Marked)
							erodedVertexList.push_back(balEdgeDest);
						(*vertices)[balEdgeDest].flagged = GraphVertex::Marked;
						(*vertices)[balEdgeSrc].marked = GraphVertex::Marked;
						(*edges)[balEdge].flagged = GraphEdge::Marked;
						erodedEdgeList.push_back(balEdge);
					}
				}
				outEdge = (*vertices)[vertexIndex].NextOut(outEdge);
			}
		}

		// Do a sanity check on the balanced edges
		balancedIndex = (*edges)[outEdgeIndex].balancedEdge;
		balancedDestVertex = (*edges)[balancedIndex].dest;
		/*
			if (balancedIndex != outEdgeIndex and balancedDestVertex > vertexIndex)
			assert((*vertices)[balancedDestVertex].OutDegree() == 0 );
		*/


	}
	if (inDegree == 0 and outDegree == 0) {
		(*vertices)[vertexIndex].flagged = GraphVertex::Marked;
		erodedVertexList.push_back(vertexIndex);
	}
}


void IntervalGraph::RemoveEdgeAndMarkPathsForRemoval(ssize_t edge,
																										 std::vector<ssize_t> &removedVertices) {

	MarkPathsThroughEdgeForRemoval(edge);
	RemoveEdge(edge, removedVertices);
}


void IntervalGraph::RemoveEdgeAndMarkIntervalsForRemoval(std::vector<ssize_t> &edgeList,
																												 std::vector<ssize_t> &removedVertices) {
	ssize_t e;
	for (e = 0; e < edgeList.size(); e++ ) {
		MarkIntervalsInEdgeForRemoval(edgeList[e]);
		RemoveEdge(edgeList[e], removedVertices);
	}
}

void IntervalGraph::RemoveEdgeAndMarkIntervalsForRemoval(ssize_t edge,
																												 std::vector<ssize_t> &removedVertices) {

	MarkIntervalsInEdgeForRemoval(edge);
	RemoveEdge(edge, removedVertices);
}


void IntervalGraph::RemoveEdge(ssize_t edge, std::vector<ssize_t> &removedVertices) {
	ssize_t sourceVertex, destVertex;
	ssize_t outIndex, inIndex;
	sourceVertex = edges[edge].src;
	destVertex   = edges[edge].dest;
	// Detach this edge if necessary
	if (sourceVertex >= 0) {
		outIndex = vertices[sourceVertex].LookupOutIndex(edge);
		assert(outIndex >= 0);
		vertices[sourceVertex].out[outIndex] = -1;
	}

	// Detach this edge from dest if necessary
	if (destVertex >= 0) {
		inIndex  = vertices[destVertex].LookupInIndex(edge);
		assert(inIndex >= 0);
		vertices[destVertex].in[inIndex] = -1;
	}

	if (sourceVertex >= 0 and
			vertices[sourceVertex].OutDegree() == 0 and 
			vertices[sourceVertex].InDegree() == 0) {
		//    std::cout << "adding " << sourceVertex  << std::endl;
		removedVertices.push_back(sourceVertex);
	}

	if (destVertex >= 0 and 
			vertices[destVertex].InDegree() == 0 and
			vertices[destVertex].OutDegree() == 0) {
		//    std::cout << "adding " << destVertex << std::endl;
		removedVertices.push_back(destVertex);
	}
	edges[edge].src = -1;
	edges[edge].dest = -1;
	/*  RemoveAllPathsThroughEdge(edge);*/
}

void IntervalGraph::RemoveAllButMST() {
	std::cout << CurTimeString() << ": RemoveAllButMST()" << std::endl;

	ssize_t e;
	//UNUSED// ssize_t v;
	GraphAlgo::CalcMST(vertices, edges);  
	std::vector<ssize_t> mstEdges;
	std::cout << "mst: ";
	for (e = 0; e < edges.size(); e++ ) {
		if (edges[e].mst == GraphAlgo::MSTIn and 
				edges[edges[e].balancedEdge].mst != GraphAlgo::MSTIn) {
			// Fix the graph to add back balanced edges
			edges[edges[e].balancedEdge].mst = GraphAlgo::MSTIn;
		}
		if (edges[e].mst == GraphAlgo::MSTIn) {
			mstEdges.push_back(e);
		}
	}
	std::cout << std::endl;
	if (CheckBalancedEdgeList(edges, mstEdges) ) {
		std::cout << "mst edges are BALANCED " << std::endl;
	}
	else {
		std::cout << "mst edges are NOT balanced " << std::endl;
	}

	std::vector<ssize_t> removedVertices, removedEdges;
	for (e = 0; e < edges.size(); e++ ) {
		if (edges[e].mst != GraphAlgo::MSTIn) {
			RemoveEdge(e, removedVertices);
			removedEdges.push_back(e);
		}
	}
	Prune(removedVertices, removedEdges);
	CondenseSimplePaths();
	Erode(60);
	std::cout << CurTimeString() << ": Exit RemoveAllButMST()" << std::endl;
}



void IntervalGraph::Erode(ssize_t minEdgeLength) {
	std::cout << CurTimeString() << ": Erode(minEdgeLength=" << minEdgeLength << ")" << std::endl;

	ErodeLeavesFunctor erodeLeaves;

	erodeLeaves.vertices = &vertices;
	erodeLeaves.edges    = &edges;
	erodeLeaves.graph    = this;
	erodeLeaves.minLength = minEdgeLength;
	std::vector<ssize_t> erodedEdgeList, erodedVertexList;
	ssize_t iteration = 0;
	std::vector<ssize_t> skippedEdges, skippedVertices;
	//UNUSED// ssize_t e;
	ssize_t numCondensedEdges = 0;
	do {
		++iteration;
		erodeLeaves.Clear();
		erodedVertexList.clear();
		erodedEdgeList.clear();
		Unmark();
		CheckGraphStructureBalance();
		TraverseRandomAccess(vertices, edges, erodeLeaves);
		//UNUSED// ssize_t i;
		ssize_t e ;

		std::sort(erodeLeaves.erodedEdgeList.begin(),
							erodeLeaves.erodedEdgeList.end());
		//UNUSED// ssize_t e2;

		assert(CheckEdges(vertices, edges));
		for (e = 0; e < erodeLeaves.erodedEdgeList.size(); e++ ) {
			if (e == 0 or 
					(e > 0 and erodeLeaves.erodedEdgeList[e] != erodeLeaves.erodedEdgeList[e-1])) {
/*				cout << "removing edge (erode ): " << e << " " << edges[e].src << " -> " << edges[e].dest
						 << "  [" << vertices[edges[e].src].index << "] ->  [" 
						 << vertices[edges[e].dest].index << "] " << endl;
*/
				MarkIntervalsInEdgeForRemoval(erodeLeaves.erodedEdgeList[e]);
				//				 DeleteEdgeReadIntervals(erodeLeaves.erodedEdgeList[e]);
				RemoveEdge(erodeLeaves.erodedEdgeList[e], erodedVertexList);
			}
		}
		assert(CheckEdges(vertices, edges));
		assert(CheckGraphStructureBalance());
		std::cout << "removing " << erodeLeaves.erodedEdgeList.size() << " short edges." << std::endl;
		Prune(erodedVertexList, erodeLeaves.erodedEdgeList);
		assert(CheckEdges(vertices, edges));

		skippedVertices.clear();
		skippedEdges.clear();
		//    CheckBalance(edges);

		assert(CheckGraphStructureBalance());
		//		RemoveMarkedIntervals();
		//		 numCondensedEdges+= CondenseSimplePaths();

		// Do a sanity check.  I might turn this off in non-debugging mode.
		CheckBalance();
		//    ClearMST();
	} while (erodedEdgeList.size() > 0);
	/*
	if (containsIntervals)
		RemoveEmptyEdges();
	*/
	numCondensedEdges = CondenseSimplePaths();
	assert(CheckAllPathsBalance(1));
	std::cout << "End of erode, condensed " << numCondensedEdges << " edges." << std::endl;
	assert(CheckGraphStructureBalance());
	std::cout << CurTimeString() << ": Exit Erode(minEdgeLength=" << minEdgeLength << ")" << std::endl;
	RemoveErasedPaths();
}

ssize_t IntervalGraph::CheckBalance() {
	ssize_t e;
	ssize_t balE;
	for (e = 0; e < edges.size(); e++ ) {
		balE = edges[e].balancedEdge;
		if (edges[balE].balancedEdge != e) {
			std::cout << "edges " << e <<" " << balE << " are not balanced " << std::endl;
			return 0;
		}
		else if (edges[balE].intervals->size() != edges[e].intervals->size()) {
			std::cout << "edges " << e << " (" << edges[e].index << ") " 
								<< balE << " ( " << edges[balE].index 
								<< ")  should have equal multiplicity " 
								<< edges[balE].intervals->size() << " " <<  edges[e].intervals->size() 
								<< " indices: " << edges[e].index << " " << edges[balE].index 
								<< std::endl;
			return 0;
		}
	}
	return 1;
}

void IntervalGraph::ClearMST() {
	ssize_t e;
	for (e = 0; e < edges.size(); e++ ) {
		edges[e].mst = GraphAlgo::MSTOut;
	}
}

void IntervalGraph::RemoveAllPathsThroughEdge(ssize_t edge) {
	//UNUSED// ssize_t interval = 0;
	ssize_t prevSize;
	std::vector<ssize_t> pathEdges, pathIntervalIndices;
	//  std::cout << "removing paths through " << edge << std::endl;

	while (edges[edge].intervals->size() > 0)  {  
		// clean up from previous run
		pathEdges.clear();
		pathIntervalIndices.clear();

		// find the path that contains the 0'th interval in edges[e].intervals  
		prevSize = edges[edge].intervals->size();
		TracePath(edge, 0, pathEdges, pathIntervalIndices);  
		// remove it  
		MarkPathForRemoval(pathEdges,pathIntervalIndices);
		// sanity check, the size of the interval set shrank by at least 1
		assert(edges[edge].intervals->size() < prevSize);
	}
}


ssize_t IntervalGraph::TracePathReverse(ssize_t curEdge, ssize_t curIndex, ssize_t &prevEdge, ssize_t &prevIndex) {
	ssize_t inEdgeIndex, inEdge;
	prevEdge = -1;
	prevIndex = -1;
	ssize_t sourceVertex;
	sourceVertex = edges[curEdge].src;
	//  std::cout << "tracing path back through " << vertices[sourceVertex].InDegree() << std::endl;
	if (IsIntervalMarkedForRemoval(curEdge, curIndex))
		return 0;

	for (inEdgeIndex = vertices[sourceVertex].FirstIn();
			 inEdgeIndex < vertices[sourceVertex].EndIn();
			 inEdgeIndex = vertices[sourceVertex].NextIn(inEdgeIndex)) {
		inEdge = vertices[sourceVertex].in[inEdgeIndex];
		prevIndex = TraceReadIntervalReverse(edges[curEdge], curIndex, edges[inEdge], vertexSize);
		if ( prevIndex >= 0) {
			prevEdge = inEdge;
			return 1;
		}
	}
	return (inEdgeIndex < vertices[sourceVertex].EndIn());
}

ssize_t IntervalGraph::RemoveMarkedIntervalsNoPaths() {
	// Simply remove intervals marked for deletion.  This assumes
	// the paths corresponding to these are removed by another 
	// method.

	ssize_t e, i;
	ssize_t numRemoved = 0;
	for (e = 0; e < edges.size(); e++) {
		ssize_t numRemovedCurEdge = 0;
		for (i = 0 ; i < edges[e].intervals->size(); i++ ) {
			if (IsIntervalMarkedForRemoval(e, i)) {
				numRemoved++;
				numRemovedCurEdge++;
			} else {
				if (numRemovedCurEdge>0) {
					(*edges[e].intervals)[i-numRemovedCurEdge] = (*edges[e].intervals)[i];
				}
			}
		}
		if (numRemovedCurEdge < (*edges[e].intervals).size())
			(*edges[e].intervals).resize(edges[e].intervals->size() - numRemovedCurEdge);
		else
			edges[e].intervals->clear();
		UpdatePathIndices(e);
	}
	return numRemoved;
}


ssize_t IntervalGraph::RemoveMarkedIntervals() {
	ssize_t e, i;
	ssize_t numRemoved = 0;
	/*
		This is called after various paths in the graph have been
		marked for removal.  This removes those paths from the graph,
		and packs the remaining lists of paths.
	*/
	for (e = 0; e < edges.size(); e++) {
		ssize_t numRemovedCurEdge = 0;
		for (i = 0 ; i < edges[e].intervals->size(); i++ ) {
			if (IsIntervalMarkedForRemoval(e, i)) {
				numRemoved++;
				numRemovedCurEdge++;
			} else {
				ssize_t pathPos = (*edges[e].intervals)[i].pathPos;
				assert(pathPos >= 0);
				paths[(*edges[e].intervals)[i].read][pathPos].index = i - numRemovedCurEdge;
				if (numRemovedCurEdge > 0) { // TODO: can we move this 2 lines earlier?
					(*edges[e].intervals)[i-numRemovedCurEdge] = (*edges[e].intervals)[i];
				}
			}
		}
		if (numRemovedCurEdge < (*edges[e].intervals).size())
			(*edges[e].intervals).resize(edges[e].intervals->size() - numRemovedCurEdge);
		else
			edges[e].intervals->clear();
	}
	return numRemoved;
}

ssize_t IntervalGraph::IsEdgeInDisjoint(ssize_t edge, ssize_t minSpanningPaths) {
	ssize_t numIntv = edges[edge].intervals->size();
	ssize_t i;
	ssize_t numSpanningPaths = 0;
	for (i = 0; i < numIntv; i++ ) {
		if ((*edges[edge].intervals)[i].pathPos > 0)
			++numSpanningPaths;
	}
	if (numSpanningPaths >= minSpanningPaths)
		return 0;
	else
		return 1;
}

ssize_t IntervalGraph::IsEdgeOutDisjoint(ssize_t edge, ssize_t minSpanningPaths) {
	ssize_t numIntv = edges[edge].intervals->size();
	ssize_t i;
	ssize_t pathIndex, pathPos;
	ssize_t numSpanningPaths = 0;
	for (i = 0; i < numIntv; i++) {
		pathIndex = (*edges[edge].intervals)[i].read;
		pathPos   = (*edges[edge].intervals)[i].pathPos;
		if (pathPos < pathLengths[pathIndex] - 1)
			++numSpanningPaths;
	}
	if (numSpanningPaths >= minSpanningPaths)
		return 0;
	else
		return 1;
}

void IntervalGraph::RemovePath(ssize_t path) {
	MarkPathForRemoval(path);
	assert(paths[path] != NULL);
	pathLengths[path] = 0;
	delete[] paths[path];
	paths[path] = NULL;
}

void IntervalGraph::MarkPathIntervalsForRemoval(ssize_t path) {
	MarkPathForRemoval(path);
}

void IntervalGraph::MarkPathForRemoval(ssize_t path) {
	ssize_t i;
	for (i = 0; i < pathLengths[path]; i++ ){ 
		MarkIntervalForRemoval(paths[path][i].edge, paths[path][i].index);
	}
}


ssize_t IntervalGraph::MarkPathRangeForRemoval(ssize_t pathIndex, ssize_t start, ssize_t end) {

	ssize_t i;
	for (i = start; i < end; i++ ){ 
		(*edges[paths[pathIndex][i].edge].intervals)[paths[pathIndex][i].index].markedForDeletion = 1;
		//		MarkIntervalForRemoval(paths[pathIndex][i].edge, paths[pathIndex][i].index);
	}
	return 0;
}


void IntervalGraph::MarkPathForRemoval(std::vector<ssize_t> &pathEdges,
																			 std::vector<ssize_t> &pathIndices) {
	ssize_t i;
	for (i = 0; i < pathEdges.size(); i++ ) {
		assert(pathIndices[i] < edges[pathEdges[i]].intervals->size());
		MarkIntervalForRemoval(pathEdges[i], pathIndices[i]);
	}
}

void IntervalGraph::MarkIntervalsOnPathForRemoval(ssize_t path) {
	ssize_t pi;
	for (pi = 0; pi < pathLengths[path]; pi++) {
		if (paths[path][pi].edge >= 0) {
			edges[paths[path][pi].edge].MarkIntervalForRemoval(paths[path][pi].index);
		}
	}
}


void IntervalGraph::RemoveMarkedPathIntervals() {
	ssize_t p, pi;
	for (p = 0; p < paths.size(); p++ ) {
		ssize_t curPi = 0;
		for (pi = 0; pi < pathLengths[p]; pi++ ){
			if (paths[p][pi].index != -1) {
				paths[p][curPi].index = paths[p][pi].index;
				paths[p][curPi].edge = paths[p][pi].edge;
				(*edges[paths[p][pi].edge].intervals)[paths[p][pi].index].pathPos = curPi;
				++curPi;
			}
		}				
		pathLengths[p] = curPi;
	}
}

void IntervalGraph::RemovePath(std::vector<ssize_t> &pathEdges,
															 std::vector<ssize_t> &pathIndices) {
	ssize_t i, j;
	for (i = 0; i < pathEdges.size(); i++ ) {
		assert(pathIndices[i] < edges[pathEdges[i]].intervals->size());
		edges[pathEdges[i]].intervals->erase(edges[pathEdges[i]].intervals->begin() + pathIndices[i]);
		// Now fix any indices of intervals that pass through this edge again.  If they do
		// the size of the interval list is smaller, so we delete one earlier.

		// This is a slightly slow way of doing this, but oh well, these lists are small
		for (j = i+1; j < pathEdges.size(); j++ ) {
			if (pathEdges[i] == pathEdges[j] and pathIndices[j] > pathIndices[i])
				pathIndices[j]--;
		}
	}
}

ssize_t IntervalGraph::CutDisjointEdges(ssize_t minPathSpan) {
	ssize_t e;
	ssize_t numCut = 0;
	TVertexList newVertices;
	ssize_t appVertexIndex = vertices.size();
	ssize_t newVertexIndex = 0;
	TVertex tmpV;
	for (e = 0; e < edges.size(); e++) { 
		if (vertices[edges[e].src].InDegree() > 0) {
			if (IsEdgeInDisjoint(e, minPathSpan)) {
				// cut this edge from the source vertex
				vertices[edges[e].src].EraseOutEdge(e);
				edges[e].src = appVertexIndex;
				newVertices.push_back(tmpV);
				newVertices[newVertexIndex].AddOutEdge(e);
				cout << "cutting out edge " << e << " of length: " << edges[e].length << endl;
				++appVertexIndex;
				++newVertexIndex;
				++numCut;
			}
		}
		if (vertices[edges[e].dest].OutDegree() > 0) {
			if (IsEdgeOutDisjoint(e, minPathSpan)) {
				vertices[edges[e].dest].EraseInEdge(e);
				edges[e].dest = appVertexIndex;
				newVertices.push_back(tmpV);
				newVertices[newVertexIndex].AddInEdge(e);
				cout << "cutting in edge " << e << " of length " << edges[e].length << endl;
				++appVertexIndex;
				++newVertexIndex;
				++numCut;
			}
		}
	}
	ssize_t v;
	appVertexIndex = vertices.size();
	vertices.resize(vertices.size() + newVertices.size());
	for (v = 0; v < newVertices.size(); v++ ) {
		vertices[appVertexIndex] = newVertices[v];
		appVertexIndex++;
	}
	CondenseSimplePaths();
	return numCut;
}


void IntervalGraph::AssignIntervalPathOrder() {
	ssize_t e, i, p;
	std::vector<ssize_t> pathEdges, pathIndices;
	pathLengths.resize(maxReadIndex+ 1);
	std::fill(pathLengths.begin(), pathLengths.end(), 0);
	paths.resize(maxReadIndex+1);
	ssize_t readIndex;
	std::cout << CurTimeString() << ": assigning path orders " << std::endl;
	for (e = 0; e < edges.size(); e++) {
		for (i = 0; i < edges[e].intervals->size(); i++) {
			if ((*edges[e].intervals)[i].pathPos == -1) {

				TracePath(e, i, pathEdges, pathIndices);

				assert(pathEdges.size() > 0);
				readIndex = (*edges[e].intervals)[i].read;
				pathLengths[readIndex] = pathEdges.size();
				paths[readIndex] = new PathInterval[pathEdges.size()];

				for (p = 0; p < pathEdges.size(); p++ ) {
					(*edges[pathEdges[p]].intervals)[pathIndices[p]].pathPos = p;
					paths[readIndex][p].edge = pathEdges[p];
					paths[readIndex][p].index = pathIndices[p];
				}

				pathEdges.clear();
				pathIndices.clear();
			}
		}
	}
	std::cout << CurTimeString() << ": done assigning path orders." << std::endl;
}




void IntervalGraph::TracePath(ssize_t curEdge, ssize_t curIndex, 
															std::vector<ssize_t> &pathEdges,
															std::vector<ssize_t> &pathIndices) {
	ssize_t curEdgeCopy = curEdge;
	ssize_t curIndexCopy = curIndex;
	//UNUSED// ssize_t nextEdge, nextIndex;
	//UNUSED// ssize_t prevEdge = -1;
	//UNUSED// ssize_t prevIndex = -1;
	ssize_t printInterval = 0;
	if (IsIntervalMarkedForRemoval(curEdge, curIndex))
		return; // 0;

	/* 
		 Starting at some random point in a path (all paths through an edge) 
		 follow the path backwards (to the beginning) and forwards (to the end).
		 Store the paths.
	*/		 
	StorePathReverse(curEdge, curIndex, pathEdges, pathIndices);
	curEdge = curEdgeCopy;
	curIndex = curIndexCopy;
	(*edges[curEdge].intervals)[curIndex].traversed = 1;
	pathEdges.push_back(curEdge);
	pathIndices.push_back(curIndex);

	StorePathForwards(curEdge, curIndex, pathEdges, pathIndices);

	if (printInterval) {
		ssize_t i;
		for (i = 0; i < pathEdges.size(); i++) {
			std::cout << pathEdges[i] << " " << pathIndices[i] << std::endl;
		}
	}
}

void IntervalGraph::CalcEdgeMultiplicityStats(double &intervalsPerNucleotide) {

	ssize_t totalLength; 		
	ssize_t totalIntervals;
	//UNUSED// ssize_t i;
	ssize_t e;
	totalIntervals = 0;
	totalLength    = 0;
	for (e = 0; e < edges.size(); e++) {
		totalIntervals += edges[e].multiplicity;
		totalLength    += edges[e].length;
	}
	intervalsPerNucleotide =  totalIntervals / double(totalLength);
}

ssize_t IntervalGraph::RemoveEmptyEdges() {
	ssize_t e;
	std::vector<ssize_t> edgesToRemove, verticesToRemove;
	for (e = 0; e < edges.size() ; e++ ) {
		if (edges[e].intervals->size() == 0 and !RemovingEdgeCutsGraph(e) and 
				edges[edges[e].balancedEdge].intervals->size() and
				!RemovingEdgeCutsGraph(edges[e].balancedEdge)) {
			edgesToRemove.push_back(e);
			RemoveEdge(e, verticesToRemove);
		}
	}
	Prune(verticesToRemove, edgesToRemove);
	return edgesToRemove.size();
}

ssize_t IntervalGraph::RemoveEmptyVertices() {
	std::vector<ssize_t> verticesToRemove, blank;
	ssize_t v;
	for (v = 0; v < vertices.size(); v++ ){ 
		if (vertices[v].InDegree() == 0 and vertices[v].OutDegree() == 0) 
			verticesToRemove.push_back(v);
	}
	Prune(verticesToRemove, blank);
	return verticesToRemove.size();
}

ssize_t IntervalGraph::RemoveLowCoverageEdges(double lowCoverageStddev, ssize_t absoluteCutoff) {
	std::cout << CurTimeString() << ": RemoveLowCoverageEdges(lowCoverageStddev=" << lowCoverageStddev << ", absoluteCutoff=" << absoluteCutoff << ")" << std::endl;

	ssize_t e;
	std::vector<ssize_t> edgesToRemove, verticesToRemove;
	ssize_t removedAnEdge = 1;
	ssize_t iter = 0;
	ssize_t prevRemoveSize;
	/* 
		 Removing some edges lowers the multiplicity of others.  
		 Keep removing edges until none are removed (removedAnEdge==0).
	*/
	//UNUSED// ssize_t srcVertex, destVertex;

	//UNUSED// double stddevMultip;
	double avgMultip ;
	CalcEdgeMultiplicityStats(avgMultip);

	while (removedAnEdge) {
		edgesToRemove.clear();
		verticesToRemove.clear();
		++iter;
		do {
			UntraverseReadIntervals();
			prevRemoveSize = edgesToRemove.size();
			ssize_t removedEdge = 0;
			for (e = 0; e < edges.size(); e++ ) {
				if (edgesToRemove.size() > 0 and 
						removedEdge < edgesToRemove.size() and 
						e == edgesToRemove[removedEdge]) 
					++removedEdge;
				else {
					ssize_t lowCov = 1;
					if (edges[e].IsUnexpectedlyLowCoverage(avgMultip, 
																								 lowCoverageStddev, 
																								 absoluteCutoff,
																								 (ssize_t)(vertexSize*3))
							or edges[e].intervals->size() < lowCov) {
						// The edge doesn't have enough reads mapped to it.
						// Maybe remove it.

						// Don't keep empty edges
						/*
						if (edges[e].intervals->size() == 0) {
														edgesToRemove.push_back(e);
						}
						*/
						// Make sure that removing this edge doesn't
						// cut the graph in two
            ssize_t balrme = edges[e].balancedEdge;
					  if (edges[balrme].IsUnexpectedlyLowCoverage(avgMultip, lowCoverageStddev, absoluteCutoff, vertexSize*3)) {
						   //cout << "bal is low too!" << endl;
            }
   					else {
							cout << "bal is HIGH " << edges[e].intervals->size() << " " << edges[balrme].intervals->size() << endl; }

						if (!RemovingEdgeCutsGraph(e)  || edges[e].intervals->size() < lowCov) {
							// Check to see if this edge creates a short cycle. 
							// We trust those edges, and do special cycle-straightening later on.	
/*							cout << "edge: " << e << " multip: " << edges[e].intervals->size()
	<< " baledge: " << edges[e].balancedEdge << " " << edges[edges[e].balancedEdge].intervals->size() << endl;
              cout << endl;
*/
							ssize_t doe, doi;
							ssize_t formsCycle = 0;
							for (doi = vertices[edges[e].dest].FirstOut();
									 doi < vertices[edges[e].dest].EndOut();
									 doi = vertices[edges[e].dest].NextOut(doi)) {
								doe = vertices[edges[e].dest].out[doi];
								if (edges[doe].dest == edges[e].src) {
									//std::cout << "edge " << e << " forms a short cycle"
									// << edges[e].length << std::endl;
									formsCycle = 1;
								}
							}
							if (!formsCycle) {
								edgesToRemove.push_back(e);
							}
						}
					}
				}
			}
			std::sort(edgesToRemove.begin(), edgesToRemove.end());
		} while (prevRemoveSize < edgesToRemove.size());
		std::cout << "Remove low coverage, iter: " << iter << " " 
							<< edgesToRemove.size() << " edges " << std::endl;

		for (e = 0; e < edgesToRemove.size(); e++) {

			RemoveEdgeAndMarkPathsForRemoval(edgesToRemove[e], verticesToRemove);
		}

		// Remove all intervals that were marked
		// when an edge was axed.  The paths have already been
		// removed, so don't bother removing those here.
		//		RemoveMarkedIntervalsNoPaths();
		Prune(verticesToRemove, edgesToRemove);
		//		assert(CheckPathContinuity());
		//    int numCondensed = CondenseSimplePaths();
		assert(CheckAllPathsBalance(1));

		//		UntraverseReadIntervals();
		// All the paths should have been removed already.
		//		RemoveMarkedIntervals();
/*
		if (vertices.size() > 0) 
			Erode(vertices[0].vertexSize * 3);
*/
		if (edgesToRemove.size() > 0)
			removedAnEdge = 1;
		else
			removedAnEdge = 0;
	}
	ssize_t numCondensed = CondenseSimplePaths();
	
	std::cout << "After removing low coverage edges, condensed " 
						<< numCondensed << " edges." << std::endl;
	std::cout << CurTimeString() << ": Exit RemoveLowCoverageEdges(lowCoverageStddev=" << lowCoverageStddev << ", absoluteCutoff=" << absoluteCutoff << ")" << std::endl;
	return 1;
}


ssize_t IntervalGraph::TraceBackEdges(ssize_t startVertex, ssize_t endVertex, ssize_t backEdges[], 
																	std::vector<ssize_t> &path) {
	ssize_t curVertex = endVertex;
	while (curVertex != startVertex) {
		path.push_back(backEdges[curVertex]);
		curVertex = edges[backEdges[curVertex]].src;
	}
	// Reverse the order of the path.	
	ssize_t i;
	ssize_t middle = path.size() / 2;
	ssize_t length = path.size();
	ssize_t temp;
	for (i = 0; i < middle; i++) {
		temp = path[i];
		path[i] = path[length - 1 - i];
		path[length - 1 - i] = temp;
	}
	return path.size();
}


ssize_t IntervalGraph::StorePathReverse(ssize_t curEdge, ssize_t curIndex,
																		std::vector<ssize_t> &pathEdges,
																		std::vector<ssize_t> &pathIndices) {
	ssize_t prevEdge, prevIndex;
	do {
		if (IsIntervalMarkedForRemoval(curEdge, curIndex))
			break;
		TracePathReverse(curEdge, curIndex, prevEdge, prevIndex);
		if (prevEdge >= 0) {
			(*edges[prevEdge].intervals)[prevIndex].traversed = 1;
			pathEdges.insert(pathEdges.begin(),prevEdge);
			pathIndices.insert(pathIndices.begin(),prevIndex);
		}
		curEdge = prevEdge;
		curIndex = prevIndex;
	} while (curEdge >= 0 and curIndex >= 0);
	return pathEdges.size();
}

ssize_t IntervalGraph::StorePathForwards(ssize_t curEdge, ssize_t curIndex,
																		 std::vector<ssize_t> &pathEdges,
																		 std::vector<ssize_t> &pathIndices) {
	ssize_t nextEdge, nextIndex;
	do {
		// Don't try and add the current edge if it is already marked
		// for removal.
		if (IsIntervalMarkedForRemoval(curEdge, curIndex))
			break;

		TracePathForwards(curEdge, curIndex, nextEdge, nextIndex);
		if (nextEdge >= 0 and nextIndex >= 0) {
			(*edges[nextEdge].intervals)[nextIndex].traversed = 1;
			pathEdges.push_back(nextEdge);
			pathIndices.push_back(nextIndex);
		}
		curEdge = nextEdge;
		curIndex = nextIndex;
	} while (curEdge >= 0 and curIndex >= 0);
	return pathEdges.size();
}

ssize_t IntervalGraph::TracePathForwards(ssize_t curEdge, ssize_t curIndex, ssize_t &nextEdge, ssize_t &nextIndex) {
	ssize_t outEdgeIndex, outEdge;
	nextEdge = -1;
	nextIndex = -1;
	ssize_t destVertex;
	//  std::cout << "tracing path forwards through " << vertices[destVertex].OutDegree() << std::endl;
	// if the edge is a self-loop, it's possible the path continues on the same edge
	destVertex = edges[curEdge].dest;
	if (IsIntervalMarkedForRemoval(curEdge, curIndex))
		return 0;

	if (curIndex < edges[curEdge].intervals->size()-1 and
			(*edges[curEdge].intervals)[curIndex].read == (*edges[curEdge].intervals)[curIndex+1].read and
			(*edges[curEdge].intervals)[curIndex].readPos + 
			(*edges[curEdge].intervals)[curIndex].length - vertices[destVertex].vertexSize
			== (*edges[curEdge].intervals)[curIndex+1].readPos) {
		nextEdge = curEdge;
		nextIndex = curIndex + 1;
	}

	// the next index isn't the succeeding index, so it's in the next edge
	destVertex = edges[curEdge].dest;
	for (outEdgeIndex = vertices[destVertex].FirstOut();
			 outEdgeIndex < vertices[destVertex].EndOut();
			 outEdgeIndex = vertices[destVertex].NextOut(outEdgeIndex)) {
		outEdge = vertices[destVertex].out[outEdgeIndex];
		nextIndex = TraceReadIntervalForward(edges[curEdge], curIndex, edges[outEdge], vertices[destVertex].vertexSize);
		if (curEdge == outEdge and 
				nextIndex == curIndex) {
			std::cout << "TRACING stayed in the same spot, that's not good!\n";
			exit(1);
			nextEdge = -1;
			return 0;
		}
		if ( nextIndex >= 0) {
			nextEdge = outEdge;
			return 1;
		}
	}
	return outEdgeIndex < vertices[destVertex].EndOut();
}

template<typename V, typename E>
class FindSourcesAndSinksFunctor {
public:
	std::vector<ssize_t> sources;
	std::vector<ssize_t> sinks;
	std::vector<V> *vertices;
	std::vector<E> *edges;
	void Clear() {
		sources.clear();
		sinks.clear();
	}
	void operator()(ssize_t vertex) {
		ssize_t inDegree, outDegree;
		inDegree = (*vertices)[vertex].InDegree();
		outDegree = (*vertices)[vertex].OutDegree();
		if (inDegree == 1 and outDegree == 0)
			sinks.push_back(vertex);

		if (inDegree == 0 and outDegree == 1) 
			sources.push_back(vertex);

		return;
	}
};


void IntervalGraph::FindSourcesAndSinks(std::vector<ssize_t> &sources,
																				std::vector<ssize_t> &sinks) {

	Unmark();
	FindSourcesAndSinksFunctor<TVertex, TEdge> findSourceSink;
	findSourceSink.vertices = &vertices;
	findSourceSink.edges    = &edges;
	ssize_t vertex;
	for (vertex = 0; vertex < vertices.size(); vertex++ ) {
		if (vertices[vertex].marked == GraphVertex::Marked)
			continue;

		// reset any previous result of sources and sinks
		findSourceSink.Clear();
		TraverseDFS(vertices, edges, vertex, findSourceSink);
		sources.insert(sources.end(), 
									 findSourceSink.sources.begin(), findSourceSink.sources.end());
		sinks.insert(sinks.end(),
								 findSourceSink.sinks.begin(), findSourceSink.sinks.end());
	}
}


void IntervalGraph::TraceSourceToSinkPaths(std::vector<ssize_t> &edgeTraversalCount) {
	// reset the markings on the graph
	Unmark();

	ssize_t vertex, edge;

	// Initialize the edge traversal count
	edgeTraversalCount.resize(edges.size());
	for (edge = 0; edge < edgeTraversalCount.size(); edge++) edgeTraversalCount[edge] = 0;

	FindSourcesAndSinksFunctor<TVertex, TEdge> findSourceSink;
	findSourceSink.vertices = &vertices;
	findSourceSink.edges    = &edges;
	//UNUSED// ssize_t source, sink;
	for (vertex = 0; vertex < vertices.size(); vertex++ ) {
		if (vertices[vertex].marked == GraphVertex::Marked)
			continue;

		// reset any previous result of sources and sinks
		findSourceSink.Clear();

		TraverseDFS(vertices, edges, vertex, findSourceSink);

		if (findSourceSink.sources.size() != findSourceSink.sinks.size()) {
			std::cout << "Warning: the graph does not have an even number of sources or sinks." << std::endl;
			std::cout << "The results are not well defined for this " << std::endl;
		}
	}
}


void IntervalGraph::IncrementOptimalPathCount(ssize_t source, 
																							std::vector<ssize_t> &sinks, 
																							std::vector<ssize_t> &edgeTraversalCount) {

	ssize_t sink;
	std::vector<ssize_t> shortestPathIndices;
	//  std::cout << "finding shortest paths from " << source << std::endl;
	SingleSourceMaximumPath(vertices, edges, source, shortestPathIndices);
	ssize_t curVertex;
	ssize_t edge;
	for (sink = 0; sink < sinks.size(); sink++) {
		curVertex = sinks[sink];
		ssize_t its = 0;
		while (curVertex != source and 
					 shortestPathIndices[curVertex] != -1) {
			assert(its < vertices.size());
			edge = vertices[curVertex].in[shortestPathIndices[curVertex]];
			edgeTraversalCount[edge]++;
			curVertex = edges[edge].src;
			its++;
		}
	}
}

void IntervalGraph::FindHighestScoringPathSourceToSink(ssize_t source, ssize_t sink, 
																											 std::vector<ssize_t> &edgeTraversalCount) {
	std::vector<ssize_t> shortestPathIndices;
	SingleSourceMaximumPath(vertices, edges, source, shortestPathIndices);

	TraceOptimalPath(vertices, edges, source, sink, shortestPathIndices, edgeTraversalCount);
}


void IntervalGraph::MarkBalancedEdges() {
	ssize_t edge;
	ssize_t balancedEdge;
	for (edge = 0; edge < edges.size(); edge++) {
		if (edges[edge].marked == GraphEdge::Marked) {
			balancedEdge = edges[edge].balancedEdge;
			if (balancedEdge >= 0) {
				edges[balancedEdge].marked = GraphEdge::Marked;
			}
		}
	}
}

ssize_t IntervalGraph::SearchForDirectedCycle(ssize_t sourceVertex, ssize_t curVertex, std::set<ssize_t> &visitedEdges,
																					ssize_t curPathLength, ssize_t maxCycleLength) {

	std::vector<ssize_t> sortedEdgeList;
	vertices[curVertex].GetSortedOutEdgeList(edges, sortedEdgeList);

	//UNUSED// ssize_t outEdgeIndex;
	ssize_t outEdge ;
	ssize_t e;
	for (e =0 ; e < sortedEdgeList.size(); e++ ){
		outEdge = sortedEdgeList[e];
		if (visitedEdges.find(outEdge) != visitedEdges.end())
			continue;
		//		edges[outEdge].traversed = GraphEdge::Marked;
		if (curPathLength + edges[outEdge].length - vertices[edges[outEdge].src].vertexSize < maxCycleLength) {
			if (SearchForDirectedCycle(sourceVertex, edges[outEdge].dest, visitedEdges,
																 curPathLength + edges[outEdge].length - vertices[edges[outEdge].src].vertexSize, 
																 maxCycleLength)) {
				return 1;
			}
		}
	}
	return 0;
}


ssize_t IntervalGraph::ExtractMin(set<ssize_t> &keyQueue, map<ssize_t,ssize_t> &values, ssize_t &minKey, ssize_t &minValue) {
	set<ssize_t>::iterator queueIt, minIt;
	// Find the closest vertex.
	assert(keyQueue.size() > 0);
	assert(values.size() > 0);
	minValue = -1;
	for (queueIt = keyQueue.begin(); queueIt != keyQueue.end(); ++queueIt) {
		if (minValue == -1 or values[*queueIt] < minValue) {
			minValue = values[*queueIt];
			minKey   = *queueIt;
			minIt    = queueIt;
		}
	}
	// pop this element.
	keyQueue.erase(minIt);
	return keyQueue.size();
}



// Search for bulges (two parallel directed paths) emanating from same vertex.
// 'vertex' is not passed into this routine, but it exists.
//
// inputs:
// redVertex, blackVertex: there is one edge from 'vertex' to 'redVertex',
//         and another edge from 'vertex' to 'blackVertex'
// bulgeLength: searching for bulges of length <= bulgeLength
//         TODO: figure out how cycle length is defined
//          (ambiguous whether it includes initial & final vertices,
//           initial red & black edges, or not)
// redNodeList, blackNodeList:
//         List of vertices encountered so far in red (or black) search
//         that are within bulgeLength from the start of the search
// redNodeLengths, blackNodeLengths:
//         redNodeLengths[v] = length in nucleotides from 'vertex' to 'v'
//         TODO: figure out ambiguities described above
// redInv, blackInv:
//         redInv[v] = invocation on which 'v' last encountered in red search.
//         Only vertices encountered on this invocation have been visited
//         This way we don't have to form a new array of visited vertices
//         every time this routine is called.
//
// maxDegree: if >0 then prune search from vertices that have outDegree > maxDegree
//         (and increment hasAboveMax when we do this)
//         if =0 then do not prune search for that reason
//
// return:
//   return value: 1 if cycle found, 0 if not found
//
//   hasAboveMax: incremented by the number of vertices that were skipped
//         for having outDegree > maxDegree
//
//   if a cycle is not found, values below are not meaningful.
//   if a cycle is found, they are meaningful:
//      minVertex: bulge goes from 'vertex' to 'minVertex'
//      minCycleLength: total size of the bulge (TODO: ambiguities above)
//      redPath, blackPath: two parallel paths from 'vertex' to 'minVertex'
//          (specified as lists of edges) that are supposed to achieve
//          the discovered cycle length
//          TODO: they don't seem to consistently achieve it
//


ssize_t IntervalGraph::SearchTwoDirectedCycle(ssize_t redVertex, ssize_t blackVertex, ssize_t maxLength,
																					ssize_t maxDegree, ssize_t &hasAboveMax,
																					BinomialHeapNode<ssize_t,ssize_t>* redNodeList[],
																					BinomialHeapNode<ssize_t,ssize_t>* blackNodeList[],
																					ssize_t redNodeLengths[],
																					ssize_t blackNodeLengths[],
																					ssize_t redInv[], ssize_t blackInv[], ssize_t invocation,
																					ssize_t &minVertex, ssize_t &minCycleLength,
																					ssize_t redPath[], ssize_t blackPath[]) {
	//
	// Initialize the queues, should this assume red != black?
	//
	BinomialHeap<ssize_t,ssize_t> redQueue;
	BinomialHeap<ssize_t,ssize_t> blackQueue;
	

	redNodeList[redVertex]             = redQueue.Insert(0, redVertex);
	blackNodeList[blackVertex]         = blackQueue.Insert(0, blackVertex);
	redNodeList[redVertex]->extRef     = &redNodeList[redVertex];
	blackNodeList[blackVertex]->extRef = &blackNodeList[blackVertex];

	std::set<ssize_t> emptyBlacklist;
	/*	std::map<int, BinomialHeapNode<int, int>* > redVertexSubset, blackVertexSubset;
	redVertexSubset[redVertex] = redNodeList[redVertex];
	blackVertexSubset[blackVertex] = blackNodeList[blackVertex];
	*/
	//	cout << "Searching for cycle from " << redVertex << " " << blackVertex << endl;
	ssize_t foundCycle = 0;
	BinomialHeapNode<ssize_t,ssize_t> *minRed = NULL, *minBlack = NULL;
	while (!redQueue.Empty() and
				 !blackQueue.Empty() and 
				 !foundCycle) {
		
		minRed   = redQueue.ExtractMin();
		minBlack = blackQueue.ExtractMin();

		// Record the distance to these nodes
		redNodeLengths[minRed->value] = minRed->key;
		blackNodeLengths[minBlack->value] = minBlack->key;

		// keep track of first visited vertices.
		redInv[minRed->value]     = invocation;
		blackInv[minBlack->value] = invocation;
		
		ssize_t cycleLength;
		assert(minRed != NULL);
		assert(minBlack != NULL);

		// Look to see if the search from the red vertex has 
		// crossed the search from the black vertex.
		if (blackInv[minRed->value] == invocation) {
			// The two searches have crossed, look to see if they form a cycle.
			cycleLength = minRed->key + blackNodeLengths[minRed->value];

			if (minCycleLength == -1 or cycleLength < minCycleLength) {
				minCycleLength = cycleLength;
				minVertex = minRed->value;
//				cout << "found a red->black cycle of length: " << minCycleLength << " " << minVertex << endl;
				foundCycle = 1;
				break;
			}
		}

		// Similar check from a red vertex.
		//

		if (redInv[minBlack->value] == invocation) {
			cycleLength = minBlack->key + redNodeLengths[minBlack->value];

			if (minCycleLength == -1 or cycleLength <minCycleLength) {
				minCycleLength = cycleLength;
				minVertex = minBlack->value;			
//				cout << "found a black->red cycle of length: " << minCycleLength << " " << minVertex << endl;
				foundCycle = 1;
				break;
			}
		}

		// The two paths have not crossed.
		// Add the edges emanating from each vertex to the queue.
		if (maxDegree > 0
				and vertices[minRed->value].OutDegree() > maxDegree) {
			hasAboveMax++;
		} else {
			AddOutEdgeLengths(minRed->value, minRed->key, maxLength,
												redQueue, redNodeList, redInv, invocation, 0, redPath);
		}

		if (maxDegree > 0
				and vertices[minBlack->value].OutDegree() > maxDegree) {
			hasAboveMax++;
		} else {
			AddOutEdgeLengths(minBlack->value, minBlack->key, maxLength,
												blackQueue, blackNodeList, blackInv, invocation, 0, blackPath);
		}
		
		
		redNodeList[minRed->value] = NULL;
		blackNodeList[minBlack->value] = NULL;
		delete minRed; minRed = NULL;
		delete minBlack; minBlack = NULL;
	}
	//	cout << "search space: " << redVertexSubset.size() << " " <<  blackVertexSubset.size() << endl;
	// Free up the space searched in the graph.
	
	if (minRed != NULL)	delete minRed;
	if (minBlack != NULL) delete minBlack; 

	redQueue.Free();
	blackQueue.Free();
	/*	map<int, BinomialHeapNode<int, int>* >::iterator subsetIt, subsetEnd;
	subsetEnd = redVertexSubset.end();
	for (subsetIt = redVertexSubset.begin(); subsetIt != subsetEnd; ++subsetIt) {
		delete subsetIt->second;
	}
	subsetEnd = blackVertexSubset.end();
	for (subsetIt = blackVertexSubset.begin(); subsetIt != subsetEnd; ++subsetIt) {
		delete subsetIt->second;
	}																																		 
	*/

	if (foundCycle) {
		return 1;
	}
	else {
		// No cycles found.
		return 0;
	}
}


void IntervalGraph::AddOutEdgeLengths(ssize_t curVertex, ssize_t lengthToCurVertex,
																			ssize_t maxLength,
																			BinomialHeap<ssize_t,ssize_t> &lengthPQueue,
																			BinomialHeapNode<ssize_t,ssize_t> *nodeRef[],
																			ssize_t invocations[], ssize_t inv,
																			ssize_t useFlagged, ssize_t path[]) {

	// Now store or relax min edge lengths for edges reaching this vertex.
	ssize_t outEdge, outEdgeIndex;
	ssize_t dest = curVertex;
	for (outEdgeIndex = vertices[curVertex].FirstOut();
			 outEdgeIndex != vertices[curVertex].EndOut();
			 outEdgeIndex = vertices[curVertex].NextOut(outEdgeIndex)) {
		outEdge = vertices[curVertex].out[outEdgeIndex];

		// If only using marked edges (for example marked in the DMST)
		// check and continue if not marked.
		if (useFlagged and edges[outEdge].marked != GraphEdge::Marked)
			continue;
		dest = edges[outEdge].dest;

		// Only try and add edges if they are close enough.
		ssize_t arrivalLength = (lengthToCurVertex 
												 + edges[outEdge].length 
												 - vertices[edges[outEdge].dest].vertexSize);
		assert(arrivalLength >= 0);
		//		BinomialHeapNode<int,int> *newNode;

		// 
		// maxLength is the maximum cycle length.  Do not look
		// for a cycle if the length out of this edge will create
		// a cycle longer than maxLength.
		//

		if (arrivalLength <= maxLength) {

			//
			// If the dest node has not been visited this search for a cycle
			// (denoted by 'inv'), add this node (and the time to get to it)
			// to the priority queue.
			//
			if (invocations[dest] != inv) {
				nodeRef[dest]   = lengthPQueue.Insert(arrivalLength, dest);
				//				graphVertexSubset[dest] = nodeRef[dest];
				nodeRef[dest]->extRef = &nodeRef[dest];
				path[dest]            = outEdge;
				assert(edges[outEdge].src != edges[outEdge].dest);
				invocations[dest]     = inv;
			}
			else if (nodeRef[dest] != NULL and 
							 nodeRef[dest]->key  > arrivalLength) {
				//
				// move this key in the tree.
				/*				cout << "decreasing " << dest << " (" << nodeRef[dest]->key << ") " 
						 << " to " << arrivalLength << " " << inv << endl; 
				*/
				/*
				cout << "decreasing key of " << nodeRef[dest] << " " 
						 << nodeRef[dest]->key << " " << nodeRef[dest]->value 
						 << " node inv "
						 << invocations[dest] << " inv " << inv << " " << dest << endl;
				*/
				assert(edges[outEdge].src != edges[outEdge].dest);				
				lengthPQueue.DecreaseKey(nodeRef[dest], arrivalLength);
				//				cout << "key of 2747 " << nodeRef[2747] << endl;
				path[dest]      = outEdge;
			}
		}
	}
}



ssize_t IntervalGraph::FindMinimalMatchedArrivalLength(map<ssize_t,ssize_t> &setA, map<ssize_t,ssize_t> &setB, ssize_t &minElement, ssize_t&minLength) {
	
	map<ssize_t,ssize_t>::iterator itA, endA, itB, endB;
	endA = setA.end();
	endB = setB.end();
	minLength = -1;
	minElement = -1;
	for (itA = setA.begin(); itA != endA; ++itA) {
		itB = setB.find((*itA).first);
		if (itB != endB) {
			if (minLength == -1 or ((*itA).second + (*itB).second) < minLength) {
				minLength = (*itA).second + (*itB).second;
				minElement = (*itA).first;
			}
		}
	}
	if (minLength == -1) {
		return 0;
	}
	else {
		return 1;
	}
}

ssize_t IntervalGraph::SearchForUndirectedCycle2(ssize_t sourceVertex, ssize_t curVertex, std::set<ssize_t> &curPath, 
																						 ssize_t curPathLength, ssize_t maxCycleLength, ssize_t isUndirected,
																						 std::set<ssize_t> &visitedEdges, std::set<ssize_t> &visitedVertices) {

	// Make sure the path is not following a cycle to an internal
	// vertex (that will be dealt with later).
  if (visitedVertices.find(curVertex) != visitedVertices.end())
    return 0;


	// This is visiting a new vertex in the graph.  Record
	// that so we know the path.
  visitedVertices.insert(curVertex);


	// Visit in edges accor
  //UNUSED// ssize_t inEdgeIndex;
  ssize_t inEdge ;
  std::vector<ssize_t> sortedEdgeList;
  vertices[curVertex].GetSortedInEdgeList(edges, sortedEdgeList);
  ssize_t e;
  for (e = 0; e < sortedEdgeList.size(); e++ ) {
    inEdge = sortedEdgeList[e];

    // If this edge has already been visited (perhaps in the other
		// direction), don't traverse it.
    if (visitedEdges.find(inEdge) != visitedEdges.end())
      continue;

		//
		// Check to see if this in-edge completes a small cycle
		// 
    if (sourceVertex == edges[inEdge].src and 
				curPathLength + edges[inEdge].length - vertices[curVertex].vertexSize < maxCycleLength ) {
      return 1;
    }
		
		// No small cycle is found.  Record this edge as visited, and continue searching.
    visitedEdges.insert(inEdge);

		// Traverse this edge in the opposite orientation, if it is possible to find
		// a short cycle after traversing this edge.
    if (curPathLength + ( edges[inEdge].length - vertices[edges[inEdge].dest].vertexSize) < maxCycleLength) {
      curPath.insert(inEdge);
      if (SearchForUndirectedCycle2(sourceVertex, edges[inEdge].src, curPath, 
																		curPathLength + edges[inEdge].length - vertices[edges[inEdge].dest].vertexSize,
																		maxCycleLength, 1, visitedEdges, visitedVertices)) {
				return 1;
      }
    }
  }


	// Visit in edges accor
  //UNUSED// ssize_t outEdgeIndex;
  ssize_t outEdge ;
	sortedEdgeList.clear();
  vertices[curVertex].GetSortedOutEdgeList(edges, sortedEdgeList);

  for (e = 0; e < sortedEdgeList.size(); e++ ) {
    outEdge = sortedEdgeList[e];

    // If this edge has already been visited (perhaps in the other
		// direction), don't traverse it.
    if (visitedEdges.find(outEdge) != visitedEdges.end())
      continue;

		//
		// Check to see if this in-edge completes a small cycle
		// 
    if (sourceVertex == edges[outEdge].dest and 
				curPathLength + edges[outEdge].length - vertices[curVertex].vertexSize < maxCycleLength ) {
      return 1;
    }
		
		// No small cycle is found.  Record this edge as visited, and continue searching.
    visitedEdges.insert(outEdge);

		// Traverse this edge in the opposite orientation, if it is possible to find
		// a short cycle after traversing this edge.
    if (curPathLength + ( edges[outEdge].length - vertices[edges[outEdge].dest].vertexSize) < maxCycleLength) {
      curPath.insert(outEdge);
      if (SearchForUndirectedCycle2(sourceVertex, edges[outEdge].dest, curPath, 
																		curPathLength + edges[outEdge].length - vertices[edges[outEdge].src].vertexSize,
																		maxCycleLength, 1, visitedEdges, visitedVertices)) {
				return 1;
      }
    }
  }
  return 0;
}


ssize_t IntervalGraph::SearchForUndirectedCycle(ssize_t sourceVertex, ssize_t cycleEndVertex, ssize_t maxCycleLength) {
  //	std::cout << "searching for a cycle from " << sourceVertex << " to " << cycleEndVertex << std::endl;
  // Initialize the priority queue to hold the source vertex
  RBTree<ScoredVertex<ssize_t> > queue;

  ScoredVertex<ssize_t> closestVertex;

  // Use a set to count the vertices that are visisted.   Although it 
  // uses extra space and runs in n log n time, it lets us not unmark the entire tree
  // after looking for a cycle, which should only check out a very small subgraph.

  closestVertex.vertex = sourceVertex;
  closestVertex.score  = 0;
  // prepare references to scores
  std::map<ssize_t, HeapScoredVertex*> references;
  std::set<ssize_t> traversedEdges;
  references[sourceVertex] = queue.Insert(closestVertex);


  while(queue.size() > 0) {
    if (queue.Pop(closestVertex) == 0) {
      std::cout << "ERROR working with queue " << std::endl;
      exit(0);
    }
		
    //    std::cout << "got vertex " << scoredVertex.vertex << std::endl;
    // First difference from normal Dijkstra: if there are any in edges
    // that are not already part of the dag.

    // Try and relax any of the edges going out of the current
    // vertex.
    
    //UNUSED// ssize_t  inEdgeIndex;
    ssize_t inEdge;
    std::vector<ssize_t> sortedEdgeList;
    ScoredVertex<ssize_t> srcVertex;
    vertices[closestVertex.vertex].GetSortedInEdgeList(edges, sortedEdgeList);
    //		std::cout << "searching in edge list of length: " << sortedEdgeList.size() << std::endl;
    ssize_t in;
    for (in = 0; in < sortedEdgeList.size(); in++ ) {
      inEdge = sortedEdgeList[in];
      if (traversedEdges.find(inEdge) != traversedEdges.end())
				continue;
      traversedEdges.insert(inEdge);
      if (edges[inEdge].length + closestVertex.score < maxCycleLength) {
				if (edges[inEdge].src == cycleEndVertex) {
					/*
						std::cout << "found a cycle ending at edge " << inEdge << " " << edges[inEdge].src 
						<< " " << edges[inEdge].dest << " of length: " 
						<< edges[inEdge].length + closestVertex.score << std::endl;
					*/
					// CYCLE FOUND!
					return 1;
				}
				else {
					srcVertex.vertex = edges[inEdge].src;

					// It is necessary to update the distance in the queue if
					// the vertex doesn't exist in the queue, or it can be reached 
					// from 'inEdge' with less distance.
					if (references.find(srcVertex.vertex) == references.end() or 
							edges[inEdge].length + closestVertex.score < references[srcVertex.vertex]->data.score) {
						
						// it's possible this edge forms a cycle. update it in the queue
						if (references.find(srcVertex.vertex) != references.end()) {
							queue.Delete(references[srcVertex.vertex]);
						}
						
						srcVertex.score = edges[inEdge].length + closestVertex.score;
						/*
							std::cout << "vertex: " << srcVertex.vertex << " is reachable from " 
							<< srcVertex.score << std::endl;
						*/
						references[srcVertex.vertex] = queue.Insert(srcVertex);
						//						std::cout <<" queue size: " << queue.size() << std::endl;
					}
				}
      }
    }

    // process out edges

    //UNUSED// ssize_t  outEdgeIndex;
    ssize_t outEdge;
    ScoredVertex<ssize_t> destVertex;
    sortedEdgeList.clear();
    vertices[closestVertex.vertex].GetSortedOutEdgeList(edges, sortedEdgeList);
    //		std::cout << "searching out edge list of length: " << sortedEdgeList.size() << std::endl;
    ssize_t out;
    for (out = 0; out < sortedEdgeList.size(); out++ ) {
      outEdge = sortedEdgeList[out];
      if (traversedEdges.find(outEdge) != traversedEdges.end())
				continue;
      traversedEdges.insert(outEdge);

      if (edges[outEdge].length + closestVertex.score < maxCycleLength) {
				if (edges[outEdge].dest == cycleEndVertex) {
					//					std::cout << "found a cycle of length: " << edges[outEdge].length + closestVertex.score << std::endl;
					// CYCLE FOUND!
					return 0;
				}
				else {
					destVertex.vertex = edges[outEdge].dest;

					// It is necessary to update the distance in the queue if
					// the vertex doesn't exist in the queue, or it can be reached 
					// from 'outEdge' with less distance.
					if (references.find(destVertex.vertex) == references.end() or 
							edges[outEdge].length + closestVertex.score < references[destVertex.vertex]->data.score) {
						
						// it's possible this edge forms a cycle. update it in the queue
						if (references.find(destVertex.vertex) != references.end()) {
							queue.Delete(references[destVertex.vertex]);
						}
						
						destVertex.score = edges[outEdge].length + closestVertex.score;
						/*						std::cout << "vertex: " << destVertex.vertex << " is reachable from " 
							<< destVertex.score << std::endl;
						*/
						references[destVertex.vertex] = queue.Insert(destVertex);
						//						std::cout <<" queue size: " << queue.size() << std::endl;
					}
				}
      }
    }
    // No cycle has been found from 'curVertex'.
  }
  return 0;
}

void IntervalGraph::RemoveWhirls(ssize_t whirlLength) {
	std::cout << CurTimeString() << ": RemoveWhirls(whirlLength=" << whirlLength << ")" << std::endl;

	if (whirlLength <= vertexSize) {
		std::cout << CurTimeString() << ": whirlLength too short. Exit RemoveWhirls(whirlLength=" << whirlLength << ")" << std::endl;
		return;
	}

	std::vector<ssize_t> edgesToRemove, verticesToRemove;
	ssize_t edge;
	for (edge = 0; edge < edges.size(); edge++) {
		if (edges[edge].src == edges[edge].dest and 
				edges[edge].length < whirlLength) {

			ssize_t i;
			for (i = 0; i < (*edges[edge].intervals).size(); i++ ){ 
				ssize_t path, pathPos;
				path = (*edges[edge].intervals)[i].read;
				pathPos = (*edges[edge].intervals)[i].pathPos;
				SplicePathRange(path, pathPos, pathPos+1);
			}
			RemoveEdge(edge, verticesToRemove);
			edgesToRemove.push_back(edge);
		}
	}
	std::cout << "Found " << edgesToRemove.size() << " whirls." << std::endl;
	//  RemoveEdgeAndMarkIntervalsForRemoval(edgesToRemove, verticesToRemove);
	Prune(verticesToRemove, edgesToRemove);
	std::cout << CurTimeString() << ": Exit RemoveWhirls(whirlLength=" << whirlLength << ")" << std::endl;
}

ssize_t IntervalGraph::VertexMapToEdgePath(map<ssize_t,ssize_t> &vertexPaths, ssize_t startVertex,
																			 list<ssize_t> &edgePath, ssize_t dir) {
	ssize_t pathLength = 0;
	map<ssize_t,ssize_t>::iterator mapIt, mapEnd;
	mapIt = vertexPaths.find(startVertex);
	while ((*mapIt).second != -1) {
		edgePath.push_back((*mapIt).second);
		if (dir > 0) {
			mapIt = vertexPaths.find(edges[(*mapIt).second].dest);
		}
		else {
			mapIt = vertexPaths.find(edges[(*mapIt).second].src);
		}
		assert(mapIt != vertexPaths.end());
		pathLength++;
	}
	return pathLength;
}

// The graph is already balanced.
// Lookup the index of the vertex that balances 'vertex'
// by examining the edges emanating from it.
// returns -1 if can't find it (e.g., if 'vertex' is isolated)
//
// TODO: Determine if there is a similar routine elsewhere, and/or if this should
// be moved to another file

ssize_t IntervalGraph::LookupBalVertex(ssize_t vertex) {
	return LookupBalancedVertex(vertices, edges, vertex);
#if 0
	for (ssize_t outEdgeIndex = vertices[vertex].FirstOut();
			 outEdgeIndex != vertices[vertex].EndOut();
			 outEdgeIndex = vertices[vertex].NextOut(outEdgeIndex)) {
		ssize_t outEdge = vertices[vertex].out[outEdgeIndex];

		// not an edge
		if (edges[outEdge].dest == -1)
			continue;

		// found an edge.
		// vertex is the src on outEdge,
		// and balVertex is the dest on balEdge
		ssize_t balEdge = edges[outEdge].balancedEdge;
		ssize_t balVertex = edges[balEdge].dest;
		return balVertex;
	}

	return -1;
#endif
}


// Removes bulges by merging edges in cycles of length < bulgeLength
// input:
//   useDMST: if true, merges edges onto those in a directed maximum spanning tree
//            if false, selects edges with other criteria
//   iter: iteration number
//   maxDegree: if 0, ignored.
//              if >0, vertices whose outdegree is > maxDegree will be
//              skipped, but if any exist, hasAboveMax will be set
// output:
//   hasAboveMax: number of vertices skipped because OutDegree > maxDegree
//                0 if none skipped for this reason
//  return value: number of removed edges


ssize_t IntervalGraph::RemoveBulgingEdges(ssize_t bulgeLength, ssize_t useDMST, ssize_t iter, ssize_t maxDegree, ssize_t &hasAboveMax) {
	cout << CurTimeString() << ": RemoveBulgingEdges(bulgeLength = " << bulgeLength << ", useDMST = " << useDMST << ", iter = " << iter << ", maxDegree = " << maxDegree << ", hasAboveMax=" << hasAboveMax << ")" << endl; // DEBUG

  ssize_t edge;
  std::vector<ssize_t> edgesToRemove, orphanedVertices;
  std::set<ssize_t> edgesToRemoveSet;
  //	Untraverse();
  std::set<ssize_t> curPath, visitedEdges, visitedVertices;
  Untraverse();
  //UNUSED// int numNotMarked = 0;

	hasAboveMax = 0;

	map<ssize_t,ssize_t> srcVertexLengths, destVertexLengths, srcVertexPaths, destVertexPaths;

	ssize_t vertex;
	BinomialHeapNode<ssize_t, ssize_t> **redNodeList = new BinomialHeapNode<ssize_t,ssize_t>*[vertices.size()];
	BinomialHeapNode<ssize_t, ssize_t> **blackNodeList = new BinomialHeapNode<ssize_t,ssize_t>*[vertices.size()];

	cout << " allocating node lists of size " << vertices.size() << endl;
	ssize_t *redInv = new ssize_t[vertices.size()];
	ssize_t *blackInv = new ssize_t[vertices.size()];
	ssize_t *redPath = new ssize_t[vertices.size()];
	ssize_t *blackPath = new ssize_t[vertices.size()];
	ssize_t *redDistList = new ssize_t[vertices.size()];
	ssize_t *blackDistList = new ssize_t [vertices.size()];
	std::fill(redInv, redInv + vertices.size(), 0);
	std::fill(blackInv, blackInv + vertices.size(), 0);
	std::fill(redPath, redPath + vertices.size(), edges.size());
	std::fill(blackPath, blackPath + vertices.size(), edges.size());
	std::fill(redNodeList, redNodeList + vertices.size(), (BinomialHeapNode<ssize_t,ssize_t>*) NULL);
	std::fill(blackNodeList, blackNodeList + vertices.size(), (BinomialHeapNode<ssize_t,ssize_t>*) NULL);
	std::fill(redDistList, redDistList + vertices.size(), 0);
	std::fill(blackDistList, blackDistList + vertices.size(), 0);
	IntMatrix pathMat;
	IntMatrix mismatchMat;
	IntMatrix scoreMat;
	InitScoreMat(mismatchMat, 0, 1);
	ssize_t invocation = 0;

	ssize_t prevVertexHigh = 0; // DEBUG
	ssize_t prevInvocation = 0; // DEBUG
	time_t prevTime = 0; // DEBUG

	for (vertex = 0; vertex < vertices.size(); vertex++) {

		if (prevVertexHigh) { // DEBUG
			prevInvocation = invocation - prevInvocation;
			prevTime = time((time_t *) 0) - prevTime;
			cout << CurTimeString() << ": End high degree vertex = " << vertex-1 << "   OutDegree = " << vertices[vertex-1].OutDegree() << "  InDegree = " << vertices[vertex-1].InDegree() << "  invocation = " << invocation << "  net invocation = " << prevInvocation << "  net time = " << prevTime << endl; // DEBUG
			prevVertexHigh = 0; // DEBUG
		} // DEBUG


		if (vertices[vertex].OutDegree() > 1) {
			// Bulging edges may exist from this vertex.

			ssize_t balVertex = LookupBalVertex(vertex);

			if (vertices[vertex].OutDegree() > 50) { // DEBUG
				cout << CurTimeString() << ": High degree vertex = " << vertex << "  (dual " << balVertex << ")" << "   OutDegree = " << vertices[vertex].OutDegree() << "  InDegree = " << vertices[vertex].InDegree() << "  invocation = " << invocation << endl; // DEBUG
				prevVertexHigh = 1;
				prevInvocation = invocation;
				prevTime = time((time_t *) 0); // DEBUG

				// DEBUG
				// Print sequences of all edges emanating from this vertex
				std::cout << "Outgoing edges:" << std::endl;
				ssize_t eCount = -1;
				ssize_t eIndex;
				for (eIndex = vertices[vertex].FirstOut();
						 eIndex != vertices[vertex].EndOut();
						 eIndex = vertices[vertex].NextOut(eIndex)) {
					ssize_t eNum = vertices[vertex].out[eIndex];
					ssize_t dest = edges[eNum].dest;
					ssize_t bal_e = edges[eNum].balancedEdge;
					ssize_t bal_dest = edges[bal_e].src;
					eCount++;

					std::cout << "it=" << iter << " v=" << vertex << " eCount=" << eCount << " eIndex=" << eIndex
										<< " /  edge " << eNum << ":" << vertex << ">" << dest
										<< "  (dual " << bal_e << ":" << balVertex << "<" << bal_dest << ")"
										<< "  len=" << edges[eNum].length
										<< "  mult=" << edges[eNum].intervals->size()
										<< " / ";
					for (ssize_t i=0; i<edges[eNum].seq.length; i++) {
						std::cout << edges[eNum].seq.seq[i];
					}
					std::cout << std::endl;
				}
			} // DEBUG

			// skip self-dual vertices due to complexity of graph simplification with them
			if (vertex == balVertex)
				continue;

			ssize_t outEdge, outEdgeIndex;
			ssize_t markedEdgeIndex, markedEdge, numMarkedEdges;
			numMarkedEdges = 0;
			/*			cout << "searching vertex of out degree: " 
				<< vertices[vertex].OutDegree() << endl;*/
			for (outEdgeIndex = vertices[vertex].FirstOut();
					 outEdgeIndex != vertices[vertex].EndOut();
					 outEdgeIndex = vertices[vertex].NextOut(outEdgeIndex)) {
				outEdge = vertices[vertex].out[outEdgeIndex];
				if ((useDMST and 
						 edges[outEdge].marked == GraphVertex::Marked) or 
						(!useDMST and
						 edges[outEdge].length < bulgeLength)) {
					markedEdge = outEdge;
					markedEdgeIndex = outEdgeIndex;
					numMarkedEdges++;
				}
			}

			if (prevVertexHigh) { // DEBUG
				cout << "numMarkedEdges = " << numMarkedEdges << endl; // DEBUG
			} // DEBUG
			
			if ((useDMST and numMarkedEdges != vertices[vertex].OutDegree()) or
					numMarkedEdges > 0) {

				if (maxDegree > 0
						and vertices[vertex].OutDegree() > maxDegree) {
					cout << "Skipping high degree vertex " << vertex
							 << " of degree " << vertices[vertex].OutDegree() << endl;
					hasAboveMax++;
					continue;
				}


				for (outEdgeIndex = vertices[vertex].FirstOut();
						 outEdgeIndex != vertices[vertex].EndOut();
						 outEdgeIndex = vertices[vertex].NextOut(outEdgeIndex)) {
					outEdge = vertices[vertex].out[outEdgeIndex];
					
					// no edge
					if (edges[outEdge].dest == -1)
						continue;

					// skip one-edge whirls
					if (edges[outEdge].src == edges[outEdge].dest)
						continue;

					// Don't remove bulging edges that are part of the dmst.
					if ((useDMST and edges[outEdge].marked == GraphEdge::Marked) or
							(!useDMST and edges[outEdge].length >= bulgeLength))
						continue;
					// length == bulgeLength case: when the other edges are added in, length will be > bulgeLength

					ssize_t nextOutEdgeIndex, nextOutEdge;

					for(nextOutEdgeIndex = vertices[vertex].NextOut(outEdgeIndex);
							nextOutEdgeIndex != vertices[vertex].EndOut();
							nextOutEdgeIndex = vertices[vertex].NextOut(nextOutEdgeIndex)) {
						nextOutEdge = vertices[vertex].out[nextOutEdgeIndex];

						// no edge
						if (edges[nextOutEdge].dest == -1)
							continue;

						// skip one-edge whirls
						if (edges[nextOutEdge].src == edges[nextOutEdge].dest)
							continue;

						ssize_t cycleLength;
						// This is just the part of the cycle length from the initial pair of edges.
						// The full cycle, if it exists, may have more edges and be longer.
						// TODO: vertexSize is now constant, but this should be adjusted
						// to account for source/dest vertices if vertexSize is not constant.
						cycleLength =   edges[outEdge].length 
														- vertices[edges[outEdge].dest].vertexSize
														+ edges[nextOutEdge].length 
														- vertices[edges[nextOutEdge].dest].vertexSize;
						if (cycleLength > bulgeLength) {
							continue;
						}

								//						if (edges[nextOutEdge].length >= bulgeLength) {
								//							continue;
								//						}

						markedEdge = nextOutEdge;
						//...
						ssize_t outBal, nextOutBal;
						// int balVertex;
						nextOutBal = edges[outEdge].balancedEdge;
						outBal     = edges[markedEdge].balancedEdge;
						balVertex = edges[outBal].dest;
						
						if (markedEdge == nextOutBal
								or outEdge == outBal      // redundant w/previous line
								or outEdge == nextOutBal
								or markedEdge == outBal
								or balVertex == vertex    // redundant w/test at top of for loop
								or edges[nextOutBal].src == vertex
								or edges[outBal].src == vertex
								or edges[markedEdge].dest == balVertex
								or edges[outEdge].dest == balVertex
								) 
							// Don't try and merge these edges since the 
							// balance will be messed up.
							continue;

						// This is a whirl, don't remove it here.
						//						if (edges[outEdge].dest == vertex)
						//							continue;

						//...						

						//
						// If using the dMST, and this edge is marked as part of the dMST,
						// then it will not be merged with the outEdge, which is already marked
						// as part of the dMST.
						//
						if (useDMST and edges[markedEdge].marked == GraphEdge::Marked)
							continue;


						
						//					assert(outEdge != markedEdge);
						// Now look to see if a cycle leaves this edge.
						map<ssize_t,ssize_t> markedPathLengths, unmarkedPathLengths, markedPaths, unmarkedPaths;
						
						ssize_t minTwoDirCycleVertex = -1, minTwoDirCycleLength = -1;
						list<ssize_t> twoDirRedPath, twoDirBlackPath;
						

						//
						// Find an undirected cycle that is composed of n fowrward edges followed by m reverse
						// edges, where n,m>= 0, and n+m > 1.
						//
						
						++invocation;
						if (invocation == 0) {
							cout << "INVOCATION MUST BE RESET!" << endl;
							std::fill(redInv, redInv + vertices.size(), 0);
							std::fill(blackInv, blackInv + vertices.size(), 0);
							++invocation;
						}

						//						cout << "invocation = " << invocation << endl; // DEBUG
						//						cout << "outEdgeIndex = " << outEdgeIndex << "  outEdge = " << outEdge << "(" << edges[outEdge].src << " -> " << edges[outEdge].dest << ")" << endl; // DEBUG
						//						cout << "nextOutEdgeIndex = " << nextOutEdgeIndex << "  nextOutEdge = " << nextOutEdge << "(" << edges[nextOutEdge].src << " -> " << edges[nextOutEdge].dest << ")" << endl; // DEBUG
						//						cout << "  markedEdge = " << markedEdge << "(" << edges[markedEdge].src << " -> " << edges[markedEdge].dest << ")" << endl; // DEBUG

						if (SearchTwoDirectedCycle(edges[markedEdge].dest, edges[outEdge].dest, 
																			 bulgeLength,
																			 maxDegree, hasAboveMax,
																			 redNodeList, blackNodeList, 
																			 redDistList, 
																			 blackDistList,
																			 redInv, blackInv, invocation,
																			 minTwoDirCycleVertex, minTwoDirCycleLength,
																			 redPath, blackPath)) {
							
							ssize_t minVertex, minVertexLength;
							
							minVertex = minTwoDirCycleVertex;
							minVertexLength = minTwoDirCycleLength;

							if (minVertex == vertex)
								continue;

							//
							// Trace out the two paths.
							//
							std::vector<ssize_t> redOptPath, blackOptPath;

							//
							// The first edges to the start vertex are initial
							// conditions for the search.
							//
							redPath[edges[markedEdge].dest] = markedEdge;
							blackPath[edges[outEdge].dest]  = outEdge;

							TraceBackEdges(vertex, minVertex, redPath, redOptPath);
							TraceBackEdges(vertex, minVertex, blackPath, blackOptPath);

							SimpleSequence redSequence, blackSequence;

							PathToSequence(redOptPath, redSequence);
							PathToSequence(blackOptPath, blackSequence);
							/*
								cout << redOptPath.size() << " " << blackOptPath.size() << endl;
								redSequence.PrintSeq(cout, "red");
								blackSequence.PrintSeq(cout, "black");
							*/
							//							cout << "red (size,length) = (" << redOptPath.size() << "," << redSequence.length << ")" << endl; // DEBUG
							//							cout << "black (size,length) = (" << blackOptPath.size() << "," << blackSequence.length << ")" << endl; // DEBUG

							DNASequence redDNA, blackDNA;
							redDNA.seq = redSequence.seq;
							redDNA.length = redSequence.length;
							blackDNA.seq = blackSequence.seq;
							blackDNA.length = blackSequence.length;
							ssize_t alignScore;
							ssize_t *locations = NULL;
							ssize_t band = 10;
							ssize_t scoreIsHigh = 1;

							// Theoretical minimum score that BandedAlign could return.
							ssize_t minAlignScore = szabs(blackSequence.length - redSequence.length);

							// Theoretical maximum score that BandedAlign could return,
							// assuming initial and terminal k-mers of the two sequences match.
							// But after edges are merged, they don't.
							// TODO: is there a way around that?
							ssize_t maxAlignScore = std::max(blackSequence.length,redSequence.length) - 2*vertexSize;

							// If < 10% of red sequence is gaps/mismatches,
							// and < 10% of black sequence is gaps/mismatches,
							// and absolute # of gaps/mismatches is < band,
							// then the edges have sequences close enough that we'll merge them.
							// Otherwise, the sequences are quite different and we don't merge them.
							ssize_t alignScoreCutoff = std::min((ssize_t) floor(.10*std::min(redSequence.length,blackSequence.length)),
																							band);

							// If the minimum possible score is above the threshold, or the maximum
							// possible score is below it, the actual score must be too.

							if (minAlignScore >= alignScoreCutoff) {
								scoreIsHigh = 1;
								// DISABLE max because initial/terminal k-mers of the two sequences don't necessarily match
								//							} else if (maxAlignScore < alignScoreCutoff) {
								//								scoreIsHigh = 0;
								//							} else if (maxAlignScore < alignScoreCutoff) {
								//								scoreIsHigh = 0;
							} else {
								// Compute banded alignment score as # of mismatches and gaps.
								// so match=0, mismatch=gap=1.  So small score if similar sequences,
								// high score if dissimilar sequences.
								// Alignment is banded, so we do not consider gaps of length >band (>=band?),
								// instead they would become mismatches
								if (minAlignScore >= maxAlignScore) {
									alignScore = minAlignScore;
								} else {
									// need to compute the score
									alignScore = BandedAlign(redDNA, blackDNA,
																					 0,1, 1, // match, mismatch, gap scores
																					 band,
																					 locations, scoreMat, pathMat, mismatchMat);
								}

									//								if (alignScore < .10*redSequence.length and
									//										alignScore < .10*blackSequence.length and
									//										alignScore < band ) {

								if (alignScore < alignScoreCutoff) {
#ifdef VERBOSE
									cout << "alignScore: " << alignScore << " " << .10*redSequence.length << " " << .10*blackSequence.length << endl;
									redDNA.PrintlnSeq(cout);
									blackDNA.PrintlnSeq(cout);
#endif
									scoreIsHigh = 0;
								}
							}

							delete[] redSequence.seq;
							delete[] blackSequence.seq;
							// The align score is too high, don't merge these edges.
							if (scoreIsHigh) {
								//							cout << "skipping alignment." << endl;
								continue;
							}

							// TODO: vertexSize is now constant, but this should be adjusted
							// to account for source/dest vertices if vertexSize is not constant.
							cycleLength = ( edges[outEdge].length 
															- vertices[edges[outEdge].dest].vertexSize
															+ edges[markedEdge].length 
															- vertices[edges[markedEdge].dest].vertexSize
															+ minVertexLength );


							//							cout // DEBUG
							//								<< "  e1.len = " << edges[outEdge].length
							//								<< "  e2.len = " << edges[markedEdge].length
							//								<< "  e1v.size = " << vertices[edges[outEdge].dest].vertexSize
							//								<< "  e2v.size = " << vertices[edges[markedEdge].dest].vertexSize
							//								<< "  minVertexLength = " << minVertexLength
							//								<< "    cycleLength = " << cycleLength
							//								<< "    red.len = " << redSequence.length
							//								<< "  black.len = " << blackSequence.length
							//								<< endl;

							if (cycleLength > bulgeLength) {
								continue;
							}

							// TODO: redundant, may keep body and remove if, since the tests were moved earlier
							if (markedEdge != outBal and
									outEdge != nextOutBal and
									edges[nextOutBal].src != vertex and
									edges[outBal].src != vertex and
									edges[markedEdge].dest != balVertex and
									edges[outEdge].dest != balVertex) {
								ssize_t outEdgeDest = edges[outEdge].dest;
#ifdef VERBOSE
								cout << "merging to " << outEdge << " from " << markedEdge << " on vertex: " << vertex <<endl;
#endif
								MergeOutEdges(vertex, outEdge, markedEdge);
								edgesToRemoveSet.insert(markedEdge);
								edges[markedEdge].Clear();
							
								if (outEdgeDest != edges[markedEdge].dest)
									orphanedVertices.push_back(edges[markedEdge].dest);
							
								ssize_t nextOutBalSrc = edges[nextOutBal].src;
#ifdef VERBOSE
								cout << "merging in " << nextOutBal << " from " << outBal << " on vertex: " << balVertex << endl;
#endif
								MergeInEdges(balVertex, nextOutBal, outBal);
								edgesToRemoveSet.insert(outBal);
								if (nextOutBalSrc != edges[outBal].src)
									orphanedVertices.push_back(edges[outBal].src);
								edges[outBal].Clear();
							}
						}
					}
				}
			}
		}
	}
	delete[] redNodeList;
	delete[] blackNodeList;
	delete[] redInv;
	delete[] blackInv;
	delete[] redPath;
	delete[] blackPath;
	delete[] redDistList;
	delete[] blackDistList;

	std::cout << CurTimeString() << ": bulge removal compared " << invocation << " pairs" << std::endl; // DEBUG


	//  std::cout << "bulge removal found " << edgesToRemoveSet.size() << " edges in short bulges of " 
	//						<< numNotMarked << " edges in bulges of " << edges.size() << " total edges." << std::endl;

  std::cout << "bulge removal found " << edgesToRemoveSet.size()
						<< " to remove out of " << edges.size() << " total edges." << std::endl;

	ssize_t v;
	for (v = 0; v < vertices.size(); v++) {
		assert(vertices[v].CheckUniqueOutEdges());
		assert(vertices[v].CheckUniqueInEdges());
	}

  for (edge = 0; edge < edges.size(); edge++ ) {
    if (edgesToRemoveSet.find(edge) != edgesToRemoveSet.end()) {
      edgesToRemoveSet.insert(edges[edge].balancedEdge);
    }
  }
  // Turn the set into a vector
  edgesToRemove.insert(edgesToRemove.begin(), edgesToRemoveSet.begin(), edgesToRemoveSet.end());
  std::sort(edgesToRemove.begin(), edgesToRemove.end());
	
  Unflag();
  for (edge = 0; edge < edgesToRemove.size(); edge++) {
    edges[edgesToRemove[edge]].flagged = GraphEdge::Marked;
  }
	cout << "after re-routing intervals." << endl;
	ssize_t p, pi;
	/*
	for (p = 0; p < paths.size(); p++ ){ 
		for (pi = 0; pi < pathLengths[p]; pi++) {
			assert(paths[p][pi].edge != -1);
		}
	}
	*/

	// Unlink edges from the graph that form bulges.
	//  RemoveEdgeAndMarkIntervalsForRemoval(edgesToRemove, orphanedVertices);
	// Remove intervals that were deleted when re-routing intervals through bulges.
	//	RemoveMarkedIntervalsNoPaths();
  Prune(orphanedVertices, edgesToRemove);
	assert(CheckEdges(vertices,edges));
#ifdef NDEBUG
	for (p = 0; p < paths.size(); p++ ){ 
		for (pi = 0; pi < pathLengths[p]; pi++) {
			assert(paths[p][pi].edge != -1);
		}
	}
#endif

	//  int numRerouted = RouteRemovedIntervals(bulgeLength);
	// Remove intervals marked for deletion by the prune operation.
	ssize_t numZeroDegree = RemoveZeroDegreeVertices(vertices, edges);
	cout << "removed " << numZeroDegree << " zero degree vertices." << endl;
	RemoveMarkedIntervals();
	//
	// These aren't needed when the paths link intervals.
	//  SortAllEdgeIntervalsByReadPos();
	//  UpdateAllPathIndices();

  DiscardGappedPaths();
  assert(CheckGraphStructureBalance());
  assert(CheckAllPathsBalance(1));

	cout << "after discarding bulging edges." << endl;
	//	assert(CheckPathContinuity());
  ssize_t numCondensed = CondenseSimplePaths();
	assert(CheckGraphStructureBalance());
	assert(CheckBalance());
  std::cout << "After bulge removal, condensed " << numCondensed <<" edges " << std::endl;



  /* Do some post-processing of the graph to remove some untidy edges.*/
  Erode(bulgeLength/2);
	RemoveWhirls(bulgeLength);
  RemoveMarkedIntervals();
	CondenseEdgeLists();
  assert(CheckAllPathsBalance(1));	
	std::cout << CurTimeString() << ": Exit RemoveBulgingEdges()" << std::endl; // DEBUG
  return edgesToRemove.size();
}

void IntervalGraph::CondenseEdgeLists() {
	ssize_t v;
	for (v = 0; v < vertices.size(); v++ ){ 
		vertices[v].CondenseEdgeLists();
	}
}	

ssize_t IntervalGraph::CheckPathContinuity() {
	ssize_t p, pi;
	for (p = 0; p < paths.size(); p++ ){ 
		for (pi= 0; pi < pathLengths[p]-1; pi++) {
			ssize_t curEdge = paths[p][pi].edge;
			ssize_t nextEdge = paths[p][pi+1].edge;
			ssize_t dest = edges[curEdge].dest;
			if (curEdge != nextEdge and vertices[dest].LookupOutIndex(nextEdge) == -1) {
				cout << "path: " << p << " is not continuous" << endl;
				return 0;
			}
		}
	}
	return 1;
}

// TODO: Procedure appears to be unused
ssize_t IntervalGraph::LocateEdgeOnPath(ssize_t edge, list<ssize_t> &path) {
	list<ssize_t>::iterator pathIt, pathEnd;
	//UNUSED// ssize_t index = 0;
	ssize_t pathPos = 0;
	pathIt = std::find(path.begin(), path.end(), edge);
	for (pathIt = path.begin(), pathEnd = path.end(); 
			 pathIt != pathEnd; ++pathIt, ++pathPos) {
		if (*pathIt == edge) {
			return pathPos;
		}
	}
	return -1;
}

ssize_t IntervalGraph::ListToVector(list<ssize_t> &l, vector<ssize_t> &v) {
	v.clear();
	v.insert(v.begin(), l.begin(), l.end());
	return v.size();
}

void IntervalGraph::ComputeIntervalPositions(list<ssize_t> &path, vector<ssize_t> &positions) {
	positions.resize(path.size());
	list<ssize_t>::iterator pathIt, pathEnd;
	ssize_t curPos = 0;
	ssize_t i = 0;
	for (pathIt = path.begin(); pathIt != path.end(); ++pathIt, ++i) {
		positions[i] = curPos;
		curPos += edges[*pathIt].length - vertices[edges[*pathIt].dest].vertexSize;
	}
}

ssize_t IntervalGraph::GetPathLength(ssize_t path) {
	ssize_t pathLength = 0;
	ssize_t i;
	ssize_t edge, intv;
	ssize_t dest;
	for (i = 0; i < pathLengths[path]-1; i++ ){
		edge = paths[path][i].edge;
		intv = paths[path][i].index;
		dest = edges[edge].dest;
		pathLength += (*edges[edge].intervals)[intv].length -
			vertices[dest].vertexSize;
	}
	edge = paths[path][i].edge;
	intv = paths[path][i].index;
	// last interval doesn't cross a vertex
	pathLength += (*edges[edge].intervals)[intv].length;
	return pathLength;
}

ssize_t IntervalGraph::ComputePathLength(vector<ssize_t> &path) {
	vector<ssize_t>::iterator pathIt, pathEnd;
	pathEnd = path.end();
	--pathEnd;
	ssize_t pathLength = 0;
	if (path.size() == 0) {
		return 0;
	}

	for (pathIt = path.begin(); pathIt != pathEnd; ++pathIt) {
		ssize_t dest = edges[*pathIt].dest;
		pathLength += edges[*pathIt].length - vertices[dest].vertexSize;
	}
	// The last edge is full-length (no overlapping vertices)
	pathLength += edges[*pathEnd].length;
	return pathLength;
}

ssize_t IntervalGraph::PathToSequence(vector<ssize_t> &path,
																	SimpleSequence &seq) {

	// First compute the length of the path.
	ssize_t pathLength = ComputePathLength(path);
	if (pathLength == 0) {
		seq.length = 0;
		return 0;
	}
	assert(pathLength > 0);
	seq.seq = new unsigned char[pathLength];
	seq.length = pathLength;
	
	vector<ssize_t>::iterator pathIt, pathEnd;
	pathEnd = path.end();
	--pathEnd;
	ssize_t curPos = 0;
	ssize_t segLength = 0; // segment length
	for (pathIt = path.begin(); pathIt != pathEnd; ++pathIt) {
		ssize_t dest, destLength;
		dest = edges[*pathIt].dest;
		destLength = vertices[dest].vertexSize;
		segLength = edges[(*pathIt)].length - destLength;
		memcpy((char*) &seq.seq[curPos], edges[(*pathIt)].seq.seq, segLength );
		curPos += segLength;
	}
	segLength = edges[(*pathEnd)].length;
	memcpy((char*) &seq.seq[curPos], edges[(*pathEnd)].seq.seq, segLength );
	return pathLength;
}

void IntervalGraph::SetMultiplicities() {
  ssize_t e;
  for (e = 0; e < edges.size(); e++) {
    edges[e].multiplicity = edges[e].intervals->size();
  }
}

void IntervalGraph::RemoveBulges(ssize_t bulgeLength, ssize_t useDMST) {
	std::cout << CurTimeString() << ": RemoveBulges(bulgeLength=" << bulgeLength << ", useDMST=" << useDMST << ")" << endl; // DEBUG

  ssize_t numRemovedEdges;
  ssize_t iter = 0;
	ssize_t hasAboveMax = 0;
	ssize_t maxDegree = 500;
	ssize_t maxDegreePass = maxDegree;

  std::stringstream pathOutStrm;

	if (containsIntervals) {
		RemoveEmptyEdges();
	}
  do {
		std::cout << CurTimeString() << ": RemoveBulges iteration " << iter << endl; // DEBUG
    Unmark();
    SetMultiplicities();
		if (useDMST) {
			cout << "Calculating the directed Maximal Spanning Tree" << endl;
			CalcDMST();
			//    std::string dmst = "dmst.dot";
			//    GVZPrintBGraph(vertices ,edges, dmst);
		  MarkBalancedEdges();
		}
		cout << "before removing bulges " << endl;
    numRemovedEdges = RemoveBulgingEdges(bulgeLength, useDMST, iter, maxDegreePass, hasAboveMax);
    pathOutStrm.str("");
    pathOutStrm << "removed."<<iter<< ".txt";
    assert(CheckAllPathsBalance(1));
		//    assert(CheckAllPathsContinuity(1));
    ++iter;
    std::cout << "Remove bulges iter: " << iter << " " << numRemovedEdges << " bulging edges " << std::endl;

		maxDegreePass = maxDegree; // Reset max for next pass
		if (hasAboveMax) {
			std::cout << "   but skipped " << hasAboveMax << " vertices of degree > " << maxDegree << std::endl;
			if (numRemovedEdges == 0) {
				// Do another pass, but this time with no limit on the max
				maxDegreePass = 0;
				numRemovedEdges = 1; // Force another pass
			}
		}
  } while (numRemovedEdges != 0);
	//	SetMultiplicities();

	std::cout << CurTimeString() << ": Exit RemoveBulges(bulgeLength=" << bulgeLength << ", useDMST=" << useDMST << ")" << endl; // DEBUG
}


void IntervalGraph::RemoveTruncatedPathIntervals() {
  // Often paths will be truncated at their ends.
  // This function erases the truncated parts, and
  // udpates the indices in edges to reference
  // back to the paths accordingly (so if the first
  // two intervals of a path are deleted, the 
  // edges that contain later parts of the path must
  // have the pathPos values for the interval decremented
  // by 2.
  ssize_t p, pre, suf, pos;
  for (p = 0; p < paths.size(); p++ ) {
    if (pathLengths[p] > 0) {
      pre = 0;
      while (pre < pathLengths[p] and
						 paths[p][pre].edge == -1)
				pre++;
      // The entire path is removed, nothing to do
      if (pre == pathLengths[p]) {
				delete[] paths[p];
				paths[p] = NULL;
				pathLengths[p] = 0;
				continue;
      }
			
      suf = pathLengths[p] - 1;
      while (suf >= pre and paths[p][suf].edge == -1)
				suf--;
			
      if (suf - pre + 1 == pathLengths[p])
				continue;

      if (suf < pre) {
				delete[] paths[p];
				paths[p] = NULL;
				pathLengths[p] = 0;
      }
      // replace paths with the truncated path
      PathInterval *newPath = new PathInterval[suf-pre+1];
      for (pos = pre; pos <= suf; pos++) {
				newPath[pos - pre] = paths[p][pos];
      }
      delete[] paths[p];
      paths[p] = newPath;
			
      pathLengths[p] = suf - pre + 1;
      // now update the position in each edge
			
      for (pos = 0; pos < pathLengths[p]; pos++) {
				if (paths[p][pos].edge != -1) {
					(*edges[paths[p][pos].edge].intervals)[paths[p][pos].index].pathPos = pos;
				}
      }
    }			
  }
}

ssize_t IntervalGraph::SearchForCycle(ssize_t sourceVertex, ssize_t prevVertex, ssize_t curVertex, 
																	ssize_t curPathLength, ssize_t maxPathLength, std::string padding) {
  
  // End case, a cycle is found
  if (curVertex == sourceVertex) {
    std::cout << "search2 found cycle " << std::endl;
    return 1;
  }

  // For now explore all paths in the graph.  Of course this
  // will take a while when allowing large paths
  ssize_t outEdge, outEdgeIndex, inEdge, inEdgeIndex;
  for (outEdgeIndex = vertices[curVertex].FirstOut(); 
       outEdgeIndex < vertices[curVertex].EndOut();
       outEdgeIndex = vertices[curVertex].NextOut(outEdgeIndex)) {
    outEdge = vertices[curVertex].out[outEdgeIndex];
    if (edges[outEdge].dest != prevVertex and
				edges[outEdge].length + curPathLength - vertexSize < maxPathLength) {
      /*
				std::cout << padding << "searching " << curVertex << " - " << outEdge << " > " << edges[outEdge].dest
				<< " " << edges[outEdge].length << " " << curPathLength << std::endl;
      */
      if (SearchForCycle(sourceVertex, curVertex, edges[outEdge].dest, 
												 edges[outEdge].length + curPathLength - vertexSize, maxPathLength, padding + " "))
				return 1;
    }
  }
  for (inEdgeIndex = vertices[curVertex].FirstIn(); 
       inEdgeIndex < vertices[curVertex].EndIn();
       inEdgeIndex = vertices[curVertex].NextIn(inEdgeIndex)) {
    inEdge = vertices[curVertex].in[inEdgeIndex];
    if (edges[inEdge].src != prevVertex and
				edges[inEdge].length + curPathLength - vertexSize < maxPathLength) {
      /*
				std::cout << padding << "searching " << curVertex << " < " << inEdge << " - " << edges[inEdge].src
				<< " " << edges[inEdge].length << " " << curPathLength << std::endl;
      */
      if (SearchForCycle(sourceVertex, curVertex, edges[inEdge].src, 
												 edges[inEdge].length + curPathLength - vertexSize, maxPathLength, padding + " ")) {
				return 1;
      }
    }
  }

  // At the end of the function, no cycle found
  return 0;
}



void IntervalGraph::ProtectEdges(std::string &protectedEdgeFileName,
																 _INT_ edgeTupleLength,
																 std::ostream &report) {
	//	std::cout << "Function ProtectEdges is deprecated." << std::endl;
  ReadPositions protPositions;
  HashValueFunctor calcHashValue;
  calcHashValue.hashLength = 10;
  HashedSpectrum hashTable(calcHashValue);

  SimpleSequenceList seqList;
  ReadSimpleSequences(protectedEdgeFileName, seqList, report);

  CountedReadPos::hashLength = edgeTupleLength;
  CountedReadPos::sequences  = &seqList;
  calcHashValue.sequences    = &seqList;
	
	hashTable.hashFunction = calcHashValue;
	hashTable.StoreSpectrum(seqList);
	//	StorSpectrum(seqList, hashTable, calcHashValue);
	//	HashToReadPositionList(hashTable, seqList, edgeTupleLength, protPositions);
	hashTable.HashToReadPositionList(seqList, protPositions);
	CompareTuples<SimpleSequenceList> comp;
  comp.sequencesPtr = &seqList;
  comp.length = edgeTupleLength;
	std::sort(protPositions.begin(), protPositions.end(), comp);

	std::cout << "there are: " << protPositions.size() << " protected positions." << std::endl;
  ProtectEdges(protPositions, seqList, edgeTupleLength);
}


void IntervalGraph::ProtectEdges(ReadPositions &protectedPositions, 
																 SimpleSequenceList &protectedEdges, 
																 _INT_ edgeTupleLength) {
  ssize_t e;
  ssize_t protIndex;
  ssize_t totalProtected = 0;
  ssize_t bale;
  ssize_t p;
  ssize_t guarded;
  for (e = 0; e < edges.size(); e++ ) {
    guarded = 0;
    for (p = 0; !guarded and p < edges[e].length - edgeTupleLength + 1; p++) {
      protIndex = LocateTuple(protectedEdges, protectedPositions, 
															edgeTupleLength, (char*) &(edges[e].seq.seq[p]));
      if (protIndex >= 0) {
				guarded = 1;
      }
    }

    if (guarded) {
      if (edges[e].guarded != GraphEdge::Marked) {
				edges[e].guarded = GraphEdge::Marked;
				//				std::cout << "protecting edge: " << e << std::endl;
				edges[e].multiplicity += 100;
				++totalProtected;
      }
      bale = edges[e].balancedEdge;
      if (edges[bale].guarded != GraphEdge::Marked) {
				//				std::cout << "protecting balance " << bale << std::endl;
				++totalProtected;
				edges[bale].multiplicity += 100;
      }
    }
  }
  std::cout << "protected: " << totalProtected << " / " << edges.size() << std::endl;
}

ssize_t IntervalGraph::CountReadsContainedInEdge(ssize_t edge){ 
  ssize_t i;
  ssize_t read, pos;
  ssize_t numContainedReads = 0;
  for (i = 0; i < edges[edge].intervals->size(); i++) {
    read = (*edges[edge].intervals)[i].read;
    pos  = (*edges[edge].intervals)[i].readPos;
		
    if (pathLengths[read] == 1) {
      numContainedReads++;
    }
  }
  return numContainedReads;
}

ssize_t IntervalGraph::CountReadsExtendingIntoEdge(ssize_t edge, ssize_t limit) {
  ssize_t i;
  ssize_t read, pos;
  ssize_t numExtendInto = 0;
  for (i = 0; i < edges[edge].intervals->size(); i++) {
    read = (*edges[edge].intervals)[i].read;
    pos  = (*edges[edge].intervals)[i].pathPos;
    if (pathLengths[read] == 1)
      continue;
    if (pos == 0 and 
				edges[edge].length - (*edges[edge].intervals)[i].edgePos > limit) {
      numExtendInto++;
    }
    else if (pos == pathLengths[read]-1 and
						 (*edges[edge].intervals)[i].length > limit) {
      numExtendInto++;
    }
  }
  return numExtendInto;
}


void IntervalGraph::RemoveLowPathEdges(ssize_t minPaths, ssize_t minExtend) {
  ssize_t e;
  std::vector<ssize_t> verticesToRemove, edgesToRemove;
  ssize_t npaths;
  for (e = 0; e < edges.size(); e++) { 
    npaths = CountReadsContainedInEdge(e) +
      CountReadsPassingThroughEdge(e) +
      CountReadsExtendingIntoEdge(e, minExtend);
    if (npaths < minPaths) {
      RemoveEdgeAndMarkPathsForRemoval(e, verticesToRemove);
      edgesToRemove.push_back(e);
    }
  }
  RemoveMarkedIntervals();
  Prune(verticesToRemove, edgesToRemove);
	RemoveEmptyEdges();
}

ssize_t IntervalGraph::CountReadsPassingThroughEdge(ssize_t edge) {
  ssize_t i;
  ssize_t read, pos;
  ssize_t numPassingThrough = 0;
  for (i = 0; i < edges[edge].intervals->size(); i++) {
    read = (*edges[edge].intervals)[i].read;
    pos  = (*edges[edge].intervals)[i].pathPos;

    if (pos > 0 and pos < pathLengths[read] - 1) {
      numPassingThrough++;
    }
  }
  return numPassingThrough;
}

void IntervalGraph::PathToThread(ssize_t p, ThreadPath &path) {
	path.clear();
	ssize_t i, pathEdge, pathInterval;
	ssize_t destVertex, destVertexLength, intvLength, intvPos;
	ssize_t readLength;
	for (i = 0; i < pathLengths[p] - 1; i++) {
		// Look up the edge, interval for this path interval
		pathEdge     = paths[p][i].edge;
		pathInterval = paths[p][i].index;
				
		// Find the length along this path interval
		destVertex       = edges[pathEdge].dest;
		destVertexLength = vertices[destVertex].vertexSize;
		intvLength       = (*edges[pathEdge].intervals)[pathInterval].length -
			destVertexLength;
		intvPos = (*edges[pathEdge].intervals)[pathInterval].edgePos;
		path.push_back(ThreadPathInterval(pathEdge, intvLength, intvPos));
		readLength += intvLength;
	}
	// The last interval includes the entire vertex length
	pathEdge     = paths[p][i].edge;
	pathInterval = paths[p][i].index;
	intvLength   = (*edges[pathEdge].intervals)[pathInterval].length;
	intvPos      = (*edges[pathEdge].intervals)[pathInterval].edgePos;
	readLength += intvLength;
	path.push_back(ThreadPathInterval(pathEdge, intvLength, intvPos));

}

void IntervalGraph::DisconnectEdgesAtSource(std::vector<ssize_t> &edgeList) {
	// Each edge in the list gets a new source vertex
	ssize_t newSrcIndex = vertices.size();
	vertices.resize(newSrcIndex + edgeList.size());

	ssize_t e;
	ssize_t origSrc;
		
	for (e = 0; e < edgeList.size(); e++) {
		origSrc = edges[edgeList[e]].src;
		ssize_t outEdgeIndex = vertices[origSrc].LookupOutIndex(edgeList[e]);
		vertices[origSrc].out[outEdgeIndex] = -1;
		vertices[newSrcIndex].AddOutEdge(edgeList[e]);
		edges[edgeList[e]].src = newSrcIndex;
		vertices[newSrcIndex].vertexSize = vertexSize;
		newSrcIndex++;
	}
}

void IntervalGraph::DisconnectEdgesAtDest(std::vector<ssize_t> &edgeList) {
	// Each edge in the list gets a new source vertex
	ssize_t newDestIndex = vertices.size();
	vertices.resize(newDestIndex + edgeList.size());

	ssize_t e;
	ssize_t origDest;
		
	for (e = 0; e < edgeList.size(); e++) {
		origDest = edges[edgeList[e]].dest;
		ssize_t inEdgeIndex = vertices[origDest].LookupInIndex(edgeList[e]);
		vertices[origDest].in[inEdgeIndex] = -1;
		vertices[newDestIndex].AddInEdge(edgeList[e]);
		edges[edgeList[e]].dest = newDestIndex;
		vertices[newDestIndex].vertexSize = vertexSize;
		newDestIndex++;
	}
}


ssize_t IntervalGraph::CountIntervalsOnSimplePath(ssize_t edge) {
	ssize_t edgeIndex, dest;

	dest = edges[edge].dest;
	ssize_t numIntervals = edges[edge].intervals->size();
	while (vertices[dest].InDegree() == 1 and
				 vertices[dest].OutDegree() == 1) {
		edgeIndex = vertices[dest].FirstOut();
		edge = vertices[dest].out[edgeIndex];
		dest = edges[edge].dest;
		numIntervals+= edges[edge].intervals->size();
	}
	return numIntervals;
}
		

void IntervalGraph::MoveIntervals(ssize_t toEdge, ssize_t fromEdge, ssize_t toEdgeIntvStartIndex, ssize_t lengthOffset) {
	ReadIntervalList *fromIntvList = edges[fromEdge].intervals,
		*toIntvList = edges[toEdge].intervals;

	ssize_t i;
	ssize_t read, index;
	
	ssize_t numFromIntervals = fromIntvList->size();
	//UNUSED// ssize_t dest             = edges[toEdge].dest;

	// 
	// Modify the edge position of the from edge intervals
	// to be relative to the start of the toEdge.
	//
	for (i = 0; i < numFromIntervals; i++) {
		(*fromIntvList)[i].edgePos += lengthOffset;
	}

	//
	// Move the interval objects.
	//
	//	toIntvList->resize(toIntvList->size() + numFromIntervals);
	std::copy(fromIntvList->begin(), fromIntvList->end(), toIntvList->begin() + toEdgeIntvStartIndex );
	
	//
	// Update the paths to reference the new edge and 
	// new position in the edge.
	//
	for (i = 0; i < fromIntvList->size(); i++) {
		if ((*fromIntvList)[i].markedForDeletion == 0) {
			read  = (*fromIntvList)[i].read;
			index = (*fromIntvList)[i].pathPos;
			paths[read][index].edge = toEdge;
			paths[read][index].index += toEdgeIntvStartIndex;
		}
	}
	fromIntvList->clear();
	edges[fromEdge].multiplicity = 0;

}




void IntervalGraph::MoveIntervals(ssize_t toEdge, ssize_t fromEdge) {
	
	// Update the paths to use 'to' where 'from' is used.
	ssize_t fromIndex = 0;
	ssize_t numFromIntervals = edges[fromEdge].intervals->size();
	ssize_t fromPath, fromPathPos;
	ReadIntervalList *intvList = edges[fromEdge].intervals;
	for (fromIndex = 0; fromIndex < numFromIntervals; fromIndex++) {
		if (!(*intvList)[fromIndex].markedForDeletion) {
			fromPath = (*intvList)[fromIndex].read;
			fromPathPos = (*intvList)[fromIndex].pathPos;
			if (fromPath != -1 and fromPathPos != -1) {
				paths[fromPath][fromPathPos].edge = toEdge;
			}
		}
	}
		
	// Append the intervals from the from edge.
	ssize_t numToIntervals = edges[toEdge].intervals->size();
	ssize_t numIntervals = numToIntervals + numFromIntervals;
	edges[toEdge].intervals->resize(numIntervals);
	ssize_t i, f;
	for (f = 0, i = numToIntervals; i < numIntervals; i++, f++) {
		(*edges[toEdge].intervals)[i] = (*edges[fromEdge].intervals)[f];
	}

	SortReadIntervalsByReadPos(*edges[toEdge].intervals);
	UpdatePathIndices(toEdge);

	// Clear the intervals on the from edge.
	edges[fromEdge].intervals->clear();

	// Set the multiplicity of the to edge.
	numToIntervals = edges[toEdge].intervals->size();
	ssize_t mult = numToIntervals;
	for (i = 0; i < numToIntervals; i++ ){ 
		if ((*edges[toEdge].intervals)[i].markedForDeletion)
			--mult;
	}
	edges[toEdge].multiplicity = mult;
}


void IntervalGraph::MergeOutEdges(ssize_t vertex,
																	ssize_t toEdge, ssize_t fromEdge) {

	ssize_t toDest = edges[toEdge].dest;
	ssize_t fromDest = edges[fromEdge].dest;

	// Remove the in-edges if necessary
	if (edges[toEdge].dest != edges[fromEdge].dest) {
		// Merge parallel edges (remove simple bulge)
		
		// Update the connectivity so that the from edge
		// is no longer used.
		ssize_t fromOutEdge, fromOutEdgeIndex;
		for (fromOutEdgeIndex = vertices[fromDest].FirstOut();
				 fromOutEdgeIndex != vertices[fromDest].EndOut();
				 fromOutEdgeIndex = vertices[fromDest].NextOut(fromOutEdgeIndex)) {
			fromOutEdge = vertices[fromDest].out[fromOutEdgeIndex];
			edges[fromOutEdge].src = toDest;
			vertices[toDest].AddOutEdge(fromOutEdge);
			vertices[fromDest].out[fromOutEdgeIndex] = -1;
		}
		ssize_t fromDestInEdge, fromDestInEdgeIndex;
		for (fromDestInEdgeIndex = vertices[fromDest].FirstIn();
				 fromDestInEdgeIndex != vertices[fromDest].EndIn();
				 fromDestInEdgeIndex = vertices[fromDest].NextIn(fromDestInEdgeIndex)) {
			fromDestInEdge = vertices[fromDest].in[fromDestInEdgeIndex];
			if (fromDestInEdge == fromEdge)
				continue;
			else {
				vertices[edges[toEdge].dest].AddInEdge(fromDestInEdge);
				edges[fromDestInEdge].dest = edges[toEdge].dest;
				vertices[fromDest].in[fromDestInEdgeIndex] = -1;
			}
		}
	}
	AppendAlternativeEdge(edges[toEdge].altEdges,
												edges[fromEdge].seq, 0);
#ifdef VERBOSE
	cout << "appending alt edge: " << endl;
	edges[toEdge].seq.PrintSeq(cout, "toEdge (out)");
	edges[fromEdge].seq.PrintSeq(cout, "fromEdge (out)");
#endif

	// Transfer intervals.  Leave lengths intact.
	MoveIntervals(toEdge, fromEdge);

	// Unlink the vertex->fromEdge, since fromEdge is deleted.
	ssize_t fromEdgeIndex = vertices[vertex].LookupOutIndex(fromEdge);
	assert(fromEdgeIndex >= 0);
	vertices[vertex].out[fromEdgeIndex] = -1;

	
	
	// Remove the in-edges if necessary
	if (edges[toEdge].dest != fromDest) {
		// Free the adjacency lists in the dest vertex.
		vertices[fromDest].out.clear();
		// hack for releasing the memory in an stl vector.
		vector<ssize_t>().swap(vertices[fromDest].out);
		vertices[fromDest].in.clear();
		// hack for releasing the memory of an stl vector.
		vector<ssize_t>().swap(vertices[fromDest].in);
	}
	else {
		fromEdgeIndex = vertices[fromDest].LookupInIndex(fromEdge);
		assert(fromEdgeIndex >= 0);
		vertices[fromDest].in[fromEdgeIndex] = -1;
	}


	edges[fromEdge].src = -1;
	edges[fromEdge].dest = -1;
}


void IntervalGraph::MergeInEdges(ssize_t vertex,
																 ssize_t toEdge, ssize_t fromEdge) {


	// Store which edge is short or long.
	ssize_t toSrc = edges[toEdge].src;
	ssize_t fromSrc = edges[fromEdge].src;
	if (edges[toEdge].src != edges[fromEdge].src) {
		// Merge parallel edges (remove simple bulge)

		// Update the connectivity so that the from edge
		// is no longer used.
		ssize_t fromInEdge, fromInEdgeIndex;
		for (fromInEdgeIndex = vertices[fromSrc].FirstIn();
				 fromInEdgeIndex != vertices[fromSrc].EndIn();
				 fromInEdgeIndex = vertices[fromSrc].NextIn(fromInEdgeIndex)) {
			fromInEdge = vertices[fromSrc].in[fromInEdgeIndex];
			edges[fromInEdge].dest = toSrc;
			vertices[toSrc].AddInEdge(fromInEdge);
			vertices[fromSrc].in[fromInEdgeIndex] = -1;
		}

		ssize_t fromSrcEdgeOutEdge, fromSrcEdgeOutEdgeIndex;
		for (fromSrcEdgeOutEdgeIndex = vertices[fromSrc].FirstOut();
				 fromSrcEdgeOutEdgeIndex != vertices[fromSrc].EndOut();
				 fromSrcEdgeOutEdgeIndex = vertices[fromSrc].NextOut(fromSrcEdgeOutEdgeIndex)) {
			fromSrcEdgeOutEdge = vertices[fromSrc].out[fromSrcEdgeOutEdgeIndex];
			// This edge is simply removed.
			if (fromSrcEdgeOutEdge == fromEdge)
				continue;
			else {
				edges[fromSrcEdgeOutEdge].src = toSrc;
				vertices[toSrc].AddOutEdge(fromSrcEdgeOutEdge);
				vertices[fromSrc].out[fromSrcEdgeOutEdgeIndex] = -1;
			}
		}
	}
	ssize_t edgeOffset = edges[toEdge].seq.length - edges[fromEdge].seq.length ;
	if (edgeOffset < 0) edgeOffset = 0;
	/*
	// TODO: Why is edgeOffset computed but not used?
	AppendAlternativeEdge(edges[toEdge].altEdges,
												edges[fromEdge].seq, 0);
	*/
	AppendAlternativeEdge(edges[toEdge].altEdges,
												edges[fromEdge].seq, edgeOffset);

#ifdef VERBOSE
	cout << "edge: " << toEdge << " has edge with offset: " << edges[toEdge].altEdges[edges[toEdge].altEdges.size() - 1].offset << endl;
	edges[toEdge].seq.PrintSeq(cout, string("to"));
	cout << "alternative:" << endl;
	edges[fromEdge].seq.PrintSeq(cout, string("alt"));
#endif
	MoveIntervals(toEdge, fromEdge);

	// Unlink the fromEdge.
	ssize_t fromEdgeIndex = vertices[vertex].LookupInIndex(fromEdge);
	assert(fromEdgeIndex >= 0);
	vertices[vertex].in[fromEdgeIndex] = -1;

	if (edges[toEdge].src != edges[fromEdge].src) {
		// fromSrc is merged with edges[toEdge].src, so 
		// the adjacency lists may be safely removed.
		vertices[fromSrc].in.clear();
		vector<ssize_t>().swap(vertices[fromSrc].in);
		vertices[fromSrc].out.clear();
		vector<ssize_t>().swap(vertices[fromSrc].out);
	}
	else {
		fromEdgeIndex = vertices[fromSrc].LookupOutIndex(fromEdge);
		assert(fromEdgeIndex >= 0);
		vertices[fromSrc].out[fromEdgeIndex] = -1;
	}
		
	//	fromEdgeIndex = vertices[fromSrc].LookupOutIndex(fromEdge);
	//	assert(fromEdgeIndex >= 0);
	//	vertices[fromSrc].out[fromEdgeIndex] = -1;
	edges[fromEdge].dest = -1;
	edges[fromEdge].src  = -1;
}


void IntervalGraph::RemoveErasedPaths() {
	std::cout << CurTimeString() << ": RemoveErasedPaths()" << std::endl;


	ssize_t p;
	ssize_t pi, curP;
	
	for (p = 0; p < paths.size(); p++) {
		curP = 0;
		for (pi = 0; pi < pathLengths[p]; pi++) {
			if (paths[p][pi].edge != -1 and paths[p][pi].index != -1) {
				paths[p][curP].edge = paths[p][pi].edge;
				paths[p][curP].index = paths[p][pi].index;
				(*edges[paths[p][curP].edge].intervals)[paths[p][pi].index].pathPos = curP;
				++curP;
			}
		}
		pathLengths[p] = curP;
	}
	
	std::cout << CurTimeString() << ": Exit RemoveErasedPaths()" << std::endl;
}

void IntervalGraph::Free() {
	
	FreePaths(paths, pathLengths);
	ssize_t e;
	for (e = 0; e < edges.size(); e++) {
			edges[e].Clear();
	}
}

void IntervalGraph::ReadAlternativeEdges(const char *altEdgeInName,
																				 std::ostream &report) {
	std::cout << CurTimeString() << ": ReadAlternativeEdges(altEdgeInName=" << altEdgeInName << ")" << std::endl;

  //
	// Read the entire file into a list of alternative edges stored at 
	// the graph level.  Though the alt edges are created and stored
	// in edges after they are created, later they must be searched
	// all at once (without looping through each edge).  Once 
	// finding an alternative edge, the position in the alternative
	// edge is then mapped back to the original edge.
	//
	ifstream altEdgeIn;
	//	if (OpenIfAvailable(string(altEdgeInName), altEdgeIn) {
	if (openck(string(altEdgeInName), altEdgeIn, std::ios::in, report, -1)) {
		cout << "Reading from " << altEdgeInName << endl;
		DNASequence seq;
		ssize_t offset;
		ssize_t destEdge;
		ssize_t lastEdge = 0;
		while (ReadAltEdge(altEdgeIn, seq, destEdge, offset)) {
			AlternativeEdge altEdge;//(seq, destEdge, offset);
			
			altEdges.push_back(AlternativeEdge());
			altEdges[lastEdge].seq.seq  = new unsigned char[seq.length];
			memcpy(altEdges[lastEdge].seq.seq, seq.seq, seq.length);
			altEdges[lastEdge].seq.length = seq.length;
			altEdges[lastEdge].destEdge = destEdge;
			altEdges[lastEdge].offset = offset;
			lastEdge++;
		}
	}
	std::cout << CurTimeString() << ": Exit ReadAlternativeEdges(altEdgeInName=" << altEdgeInName << ")" << std::endl;
}

void IntervalGraph::WriteAlternativeEdges(string &altEdgeOutName,
																					std::ostream &report) {
	std::cout << CurTimeString() << ": WriteAlternativeEdges(altEdgeOutName=" << altEdgeOutName << ")" << std::endl;
	ofstream seqOut;
	openck(altEdgeOutName, seqOut, std::ios::out, report);
	ssize_t e;
	for (e = 0; e < edges.size(); e++) {
		WriteAltEdges(edges[e].altEdges, e, seqOut);
	}
	std::cout << CurTimeString() << ": Exit WriteAlternativeEdges(altEdgeOutName=" << altEdgeOutName << ")" << std::endl;
}


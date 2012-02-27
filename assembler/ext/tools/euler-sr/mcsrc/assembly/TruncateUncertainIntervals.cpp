/***************************************************************************
 * Title:          TruncateUncertainIntervals.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  10/24/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "IntervalGraph.h"
#include "align/alignutils.h"
#include "DeBruijnGraph.h"
#include "IntervalGraph.h"
#include "SeqReader.h"
#include <string.h>
#include "IntegralTupleStatic.h"


void CreateReverseSequence(unsigned char* seq, ssize_t seqLength, 
													 unsigned char *&reverse);

ssize_t RoutePathEnd(IntervalGraph &graph,
								 ssize_t pathEdge, ssize_t routeEdge, ssize_t intervalIndex, ssize_t routeLength);

ssize_t RoutePathBegin(IntervalGraph &graph, 
									 ssize_t pathEdge, ssize_t routeEdge, ssize_t intervalIndex, ssize_t routeLength);

ssize_t ComputeFittingEditDistance(unsigned char *pathSeq, ssize_t pathSeqLength,
															 unsigned char *edgeSeq, ssize_t edgeSeqLength);


ssize_t RouteIntervalsBeginningAtEdge(IntervalGraph &graph,
																	 ssize_t edge, ssize_t routeLength, ssize_t minEditDistance);

ssize_t RouteIntervalsEndingAtEdge(IntervalGraph &graph,
																ssize_t edge, ssize_t routeLength, ssize_t minEditDistance);

void TruncatePaths(IntervalGraph &graph, ssize_t routeLength, ssize_t minEditDistance);

void PrintUsage() {

	std::cout << "usage: truncateUncertain baseIn baseOut [options] " << std::endl;
	std::cout << "            -endLength len [5] : consider ends of paths of length 'len'. " << std::endl;
	std::cout << "            -minDiff   diff [2]: any paths that end on an edge with an interval " << std::endl
						<< "                                 length 'len', and  'diff' mutations are truncated. " << std::endl; 

	exit(0);
}

int main(int argc, char* argv[]) {

	int vertexSize = 20;
  IntervalGraph graph;
  std::string baseInName, baseOutName, edgeOutName, graphOutName, edgeFileName,
		pathFileName;
  std::string graphFileName, intervalFileName, bGraphOutName, intvOutName, gvzOutName,
		pathOutName;

  ssize_t endLength = 5;
	ssize_t minDiff   = 3;
  int argi = 1;
  if (argc < 2) {
    PrintUsage();
    exit(1);
  }
  baseInName       = argv[argi++];
  baseOutName      = argv[argi++];
	while (argi < argc) {
		if (strcmp(argv[argi], "-vertexSize") == 0) 
			vertexSize = atoi(argv[++argi]);
		++argi;
	}
	
  graphFileName    = baseInName  + ".bgraph";
  intervalFileName = baseInName  + ".intv";
  edgeFileName     = baseInName  + ".edge";
	pathFileName     = baseInName  + ".path";
  intvOutName      = baseOutName + ".intv";
  bGraphOutName    = baseOutName + ".bgraph";
  graphOutName     = baseOutName + ".graph";
  edgeOutName      = baseOutName + ".edge";
  gvzOutName       = baseOutName + ".dot";
	pathOutName      = baseOutName + ".path";
  while(argi < argc) {
    if (strcmp(argv[argi], "-vertexSize") == 0) {
      vertexSize = atoi(argv[++argi]);
    }
    else if (strcmp(argv[argi], "-endLength") == 0) {
			endLength = atoi(argv[++argi]);
    }
    else if (strcmp(argv[argi], "-minDiff") == 0) {
      minDiff = atoi(argv[++argi]);
    }
    else {
      std::cout << "didn't understand " << argv[argi] << std::endl;
      PrintUsage();
      exit(1);
    }
    ++argi;
  }
  graph.vertexSize    = vertexSize;

	graph.ReadIntervalGraph(graphFileName, intervalFileName, pathFileName);
  ReadSequences(edgeFileName, graph.edges);

	Score score;
	score.CalculateEditDistance();
	

	TruncatePaths(graph, endLength, minDiff);
	graph.RemoveTruncatedPathIntervals();
	graph.RemoveMarkedIntervals();
	graph.CheckAllPathsBalance(1);
	graph.RemoveLowCoverageEdges(1, 2);
  PrintGraph(graph.vertices,graph.edges, graphOutName);
	std::cout << "printing edges to " << edgeOutName << std::endl;
	PrintEdges(graph.vertices, graph.edges, edgeOutName);
  graph.PrintIntervalGraph(bGraphOutName, intvOutName);
	WriteReadPaths(pathOutName, graph.paths, graph.pathLengths);
}


void TruncatePaths(IntervalGraph &graph, ssize_t routeLength, ssize_t minEditDistance) {
	//UNUSED+// ssize_t i;
	ssize_t e ;
	ssize_t totalRemoved = 0;
	for (e = 0; e < graph.edges.size(); e++ ) {
		totalRemoved += RouteIntervalsBeginningAtEdge(graph, e, routeLength, minEditDistance);
		totalRemoved += RouteIntervalsEndingAtEdge(graph, e, routeLength, minEditDistance);
	}
	std::cout << "truncated: " << totalRemoved << " paths " << std::endl;
}

void CreateReverseSequence(unsigned char* seq, ssize_t seqLength, 
													 unsigned char *&reverse) {
	reverse = new unsigned char[seqLength];
	ssize_t i;
	for (i = 0; i < seqLength; i++) {
		reverse[seqLength - i - 1] = seq[i];
	}
}

ssize_t RoutePathBegin(IntervalGraph &graph, 
									 ssize_t pathEdge, ssize_t routeEdge, ssize_t intervalIndex, ssize_t routeLength) {
	ssize_t destVertex = graph.edges[pathEdge].src;
	ssize_t destVertexLength = graph.vertices[destVertex].vertexSize;

	unsigned char *pathSeq, *routeSeq;

	ssize_t pathEdgeLength, routeEdgeLength;
	pathEdgeLength = graph.edges[pathEdge].length;
	routeEdgeLength = graph.edges[routeEdge].length;
	
	ssize_t routeSeqLength, pathSeqLength;
	assert(routeLength >= (*graph.edges[pathEdge].intervals)[intervalIndex].length - destVertexLength);
	pathSeqLength = (*graph.edges[pathEdge].intervals)[intervalIndex].length - destVertexLength;

	if (routeEdgeLength < routeLength + destVertexLength) 
		routeSeqLength = routeEdgeLength - destVertexLength;
	else
		routeSeqLength = routeLength;

	pathSeq  = &(graph.edges[pathEdge].seq.seq[pathEdgeLength - destVertexLength - pathSeqLength]);
	routeSeq = &(graph.edges[routeEdge].seq.seq[routeEdgeLength - destVertexLength - routeSeqLength]);
	
	unsigned char *pathSeqRev, *routeSeqRev;
	CreateReverseSequence(pathSeq, pathSeqLength, pathSeqRev);
	CreateReverseSequence(routeSeq, routeSeqLength, routeSeqRev);
	
	ssize_t editDistance;
	editDistance = ComputeFittingEditDistance(pathSeqRev, pathSeqLength,
																						routeSeqRev, routeSeqLength);
	
	delete[] pathSeqRev;
	delete[] routeSeqRev;

	return editDistance;
}

ssize_t RoutePathEnd(IntervalGraph &graph,
								 ssize_t pathEdge, ssize_t routeEdge, ssize_t intervalIndex, ssize_t routeLength) {
	ssize_t sourceVertex = graph.edges[pathEdge].src;
	ssize_t sourceVertexLength = graph.vertices[sourceVertex].vertexSize;
	
	unsigned char *pathSeq, *routeSeq;
	pathSeq  = &(graph.edges[pathEdge].seq.seq[sourceVertexLength]);
	routeSeq = &(graph.edges[routeEdge].seq.seq[sourceVertexLength]);

	ssize_t pathSeqLength, routeSeqLength;
	assert(routeLength >= (*graph.edges[pathEdge].intervals)[intervalIndex].length - sourceVertexLength);
	pathSeqLength = (*graph.edges[pathEdge].intervals)[intervalIndex].length - sourceVertexLength;

	if (routeLength > graph.edges[routeEdge].length - sourceVertexLength)
		routeSeqLength = graph.edges[routeEdge].length - sourceVertexLength;
	else
		routeSeqLength = routeLength;

	ssize_t editDistance;
	editDistance = ComputeFittingEditDistance(pathSeq, pathSeqLength,
																						routeSeq, routeSeqLength);

	return editDistance;
}

ssize_t ComputeFittingEditDistance(unsigned char *pathSeq, ssize_t pathSeqLength,
															 unsigned char *edgeSeq, ssize_t edgeSeqLength) {

	IntMatrix scores;
	CreateMatrix<ssize_t>(scores, pathSeqLength+1, edgeSeqLength+1);

	ssize_t i, j;
	ssize_t matchScore;
	ssize_t gapPath, gapEdge;
	// Initialize gap scores
	for (i = 0; i <= pathSeqLength; i++ ) 
		scores[i][0] = i;
	
	for (j = 0; j <= edgeSeqLength; j++ )
		scores[0][j] = j;

	for (i = 1; i <= pathSeqLength; i++ ) {
		for (j = 1; j <= edgeSeqLength; j++ ) {
			if (pathSeq[i-1] == edgeSeq[j-1])
				matchScore = scores[i-1][j-1];
			else
				matchScore = 1 + scores[i-1][j-1];
			
			gapPath = scores[i-1][j] + 1;
			gapEdge = scores[i][j-1] + 1;
			
			if (matchScore <= gapPath and
					matchScore <= gapEdge) 
				scores[i][j] = matchScore;
			else if (gapPath <= matchScore and gapPath <= gapEdge)
				scores[i][j] = gapPath;
			else 
				scores[i][j] = gapEdge;
		}
	}

	ssize_t minScore = edgeSeqLength;
	for (j = 0; j < edgeSeqLength+1; j++) {
		if (scores[pathSeqLength][j] < minScore) {
			minScore = scores[pathSeqLength][j];
		}
	}
	return minScore;
}


ssize_t RouteIntervalsEndingAtEdge(IntervalGraph &graph,
															 ssize_t edge, ssize_t routeLength, ssize_t minEditDistance) {
	ssize_t i;
	ssize_t readIndex;
	ssize_t pathPos;
	ssize_t srcVertex = graph.edges[edge].src;
	ssize_t srcVertexLength = graph.vertices[srcVertex].vertexSize;
	ssize_t curEdgeIndex = graph.vertices[srcVertex].LookupOutIndex(edge);
	ssize_t outEdge, outEdgeIndex;
	ssize_t routeEditDistance;
	ssize_t totalRemoved = 0;
	for (i = 0; i < graph.edges[edge].intervals->size(); i++ ) {
		readIndex = (*graph.edges[edge].intervals)[i].read;
		pathPos   = (*graph.edges[edge].intervals)[i].pathPos;
		if (graph.IsIntervalMarkedForRemoval(edge, i)) 
			continue;

		if (pathPos == graph.pathLengths[readIndex]-1) {
			// interval 'i' is a path that ends at this edge. 
			// look to see if its sequence is close to that of other
			if ((*graph.edges[edge].intervals)[i].length - srcVertexLength <= routeLength) {
				// This path ends early in the edge, check and see if we trust it.
				for (outEdgeIndex = graph.vertices[srcVertex].FirstOut();
						 outEdgeIndex < graph.vertices[srcVertex].EndOut();
						 outEdgeIndex = graph.vertices[srcVertex].NextOut(outEdgeIndex)) {
					// Don't try and route this interval in the edge it belongs to
					if (outEdgeIndex == curEdgeIndex) 
						continue;

					// otherwise, look and see if the sequence for this edge is similar to that 
					// of anohter edge
					outEdge = graph.vertices[srcVertex].out[outEdgeIndex];
					routeEditDistance = RoutePathEnd(graph, edge, outEdge, i, routeLength);					

					if (routeEditDistance < minEditDistance) {
						graph.MarkIntervalForRemoval(edge, i);
						totalRemoved++;
						/*
						std::cout << "Path: " << readIndex << " edge " << edge << " within " 
											<< routeEditDistance << " to " << outEdge << std::endl;
						*/
					}
				}
			}
		}
	}
	return totalRemoved;
}


ssize_t RouteIntervalsBeginningAtEdge(IntervalGraph &graph,
																	ssize_t edge, ssize_t routeLength, ssize_t minEditDistance) {
	ssize_t i;
	ssize_t readIndex;
	ssize_t pathPos;
	ssize_t destVertex = graph.edges[edge].dest;
	ssize_t destVertexLength = graph.vertices[destVertex].vertexSize;
	ssize_t curEdgeIndex = graph.vertices[destVertex].LookupInIndex(edge);
	ssize_t inEdge, inEdgeIndex;
	ssize_t routeEditDistance;
	ssize_t totalRemoved = 0;
	for (i = 0; i < graph.edges[edge].intervals->size(); i++ ) {
		readIndex = (*graph.edges[edge].intervals)[i].read;
		pathPos   = (*graph.edges[edge].intervals)[i].pathPos;
		if (graph.IsIntervalMarkedForRemoval(edge, i)) 
			continue;
		if (pathPos == 0) {
			// interval 'i' is a path that ends at this edge. 
			// look to see if its sequence is close to that of other
			if ((*graph.edges[edge].intervals)[i].length - destVertexLength <= routeLength) {
				// This path ends early in the edge, check and see if we trust it.
				for (inEdgeIndex = graph.vertices[destVertex].FirstIn();
						 inEdgeIndex < graph.vertices[destVertex].EndIn();
						 inEdgeIndex = graph.vertices[destVertex].NextIn(inEdgeIndex)) {
					// Don't try and route this interval in the edge it belongs to
					if (inEdgeIndex == curEdgeIndex) 
						continue;

					// otherwise, look and see if the sequence for this edge is similar to that 
					// of anohter edge
					inEdge = graph.vertices[destVertex].in[inEdgeIndex];
					routeEditDistance = RoutePathBegin(graph, edge, inEdge, i, routeLength);					
					if (routeEditDistance < minEditDistance) {
						graph.MarkIntervalForRemoval(edge, i);
						++totalRemoved;
						/*
						std::cout << "Path: " << readIndex << " e " << edge << " within " 
											<< routeEditDistance << " to " << inEdge << std::endl;
						*/
					}
				}
			}
		}
	}
	return totalRemoved;
}

/***************************************************************************
 * Title:          Scaffold.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  01/19/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "Scaffold.h"


ssize_t IsEndVertex(IntervalGraph &g, ssize_t vertexIndex) {
	return (g.vertices[vertexIndex].OutDegree() == 0 and 
					g.vertices[vertexIndex].InDegree() == 1);
}	

ssize_t IsBeginVertex(IntervalGraph &g, ssize_t vertexIndex) {
	return (g.vertices[vertexIndex].InDegree() == 0 and 
					g.vertices[vertexIndex].OutDegree() == 1);
}	


ssize_t IsEndEdge(IntervalGraph &g, ssize_t edgeIndex) {
	ssize_t dest = g.edges[edgeIndex].dest;
	return (IsEndVertex(g, dest));
}


ssize_t IsBeginEdge(IntervalGraph &g, ssize_t edgeIndex) {
	ssize_t src = g.edges[edgeIndex].src;
	return IsBeginVertex(g, src);
}

void GrowMatrices(IntMatrix &scoreMat, IntMatrix &pathMat, ssize_t nrows, ssize_t ncols) {
	if (nrows == 0 or ncols == 0) 
		return;

	if (nrows > scoreMat.size() or (scoreMat.size() > 0 and scoreMat[0].size() < ncols)) {
		ClearMatrix(scoreMat);
		ClearMatrix(pathMat);
		CreateMatrix(scoreMat, nrows, ncols);
		CreateMatrix(pathMat, nrows, ncols);
	}
}


void MatePairScaffoldJoinEdges(IntervalGraph &graph,
															 ReadMateList  &mateList,
															 ssize_t           mateType,
															 ssize_t           minMateCount,
															 IntMatrix     &matchMat,
															 vector<ssize_t>   &vToRemove) {
	std::vector<std::map<ssize_t, MateCount> > endToBegin, beginToEnd, sourceToSinkMateCount;
	endToBegin.resize(graph.edges.size());
	beginToEnd.resize(graph.edges.size());
	ssize_t mateSepMean;
	double mateSepStddev;
	ComputeMatePairLengthDistribution(graph, mateList, mateType, mateSepMean, mateSepStddev);
			
	cout << "mate sep mean: " << mateSepMean << " with stddev: " << mateSepStddev << endl;

	std::cout << "looking for mate pair overlaps." << std::endl;
			
	//
	// Count the number of mate-pairs linking each end edge to start edge.
	//
	ssize_t e;
	for (e = 0; e < graph.edges.size(); e++ ) {
				
		if (!IsEndEdge(graph, e)) {
			//			cout << "edge: " << e << " is not an end edge: " << graph.vertices[graph.edges[e].dest].OutDegree() << endl;
			// This edge is an internal edge, don't try and join it.
			continue;
		}

		ReadIntervalList *intervals;
		intervals = graph.edges[e].intervals;
		ssize_t noMate = 0;
		ssize_t intv;
		ssize_t pathIndex, pathPos;
		ssize_t readIndex, mateIndex, mateType;
		for (intv = 0; intv < intervals->size(); intv++) {
			pathIndex = (*intervals)[intv].read;
			pathPos   = (*intervals)[intv].pathPos;
			if (pathIndex % 2 != 0)
				// only process reads in the forward direction since the 
				// reverse complement is handled by the mate.
				continue;
			readIndex = pathIndex / 2;
			mateIndex = mateList[readIndex].mateIndex * 2 + 1;
			mateType  = mateList[readIndex].mateType;
			ssize_t lastEdge, lastIntv;
			ssize_t mateFirstEdge, mateFirstIntv;
			if (GetLastEdgeIndex(graph.paths, graph.pathLengths, pathIndex, lastEdge, lastIntv) == 0) {
				// The path for this read has been deleted, or something bad has happened.
				// I'm not sure if this should ever happen.
				assert(0);
			}
			if (mateIndex >= 0 ) {
				if (GetFirstEdgeIndex(graph.paths, 
															graph.pathLengths, mateIndex, mateFirstEdge, mateFirstIntv)) {

					// This may join two edges.
					if (lastEdge != mateFirstEdge) {
						ssize_t distToEnd   = graph.edges[e].length - ((*intervals)[intv].edgePos +
																											 (*intervals)[intv].length);
						//							mateFirstEdge = graph.paths[mateIndex][0].edge;
						//							mateFirstIntv = graph.paths[mateIndex][0].index;

						ssize_t distToStart = (*graph.edges[mateFirstEdge].intervals)[mateFirstIntv].edgePos;
								
						// Register this link if the mate-first edge is a source edge.
						if (!IsBeginVertex(graph, graph.edges[mateFirstEdge].src))
							continue;
						if (endToBegin[e].find(mateFirstEdge) ==
								endToBegin[e].end()) {
							endToBegin[e][mateFirstEdge].count = 1;
						}
						else {
							endToBegin[e][mateFirstEdge].count++;
							
						}
						if (beginToEnd[mateFirstEdge].find(e) ==
								beginToEnd[mateFirstEdge].end()) {
							beginToEnd[mateFirstEdge][e].count = 1;
						}
						else {
							beginToEnd[mateFirstEdge][e].count++;
						}
						endToBegin[e][mateFirstEdge].beginPos += distToEnd;
						endToBegin[e][mateFirstEdge].endPos += distToStart;
						//
						// Keep track of the number of mate-pairs linking this edge.
						// We don't want to link a source edge to multiple sink edges.
						//
					}
				}
			}
			else {
				++noMate;
			}
		}
	}

	// 
	// Compute edge overlaps.
	//
	std::map<ssize_t, MateCount>::iterator mapIt, next;
	
	IntMatrix pathMat, scoreMat;

	for (e = 0; e < graph.edges.size(); e++ ) {
	//
		// Remove links that have low count.
		//
		mapIt = endToBegin[e].begin();
		while( mapIt != endToBegin[e].end()) {
			if ((*mapIt).second.count < minMateCount) {
				next = mapIt;
				++next;
				endToBegin[e].erase(mapIt);
				mapIt = next;
			}
			else {
				++mapIt;
			}
		}
	
		mapIt = beginToEnd[e].begin();
		while( mapIt != beginToEnd[e].end()) {
			if ((*mapIt).second.count < minMateCount) {
				next = mapIt;
				++next;
				beginToEnd[e].erase(mapIt);
				mapIt = next;
			}
			else {
				++mapIt;
			}
		}
	}

	for (e = 0; e < graph.edges.size(); e++ ) {
		// Not a sink edge, or nothing stored for it.
		if (endToBegin[e].size() == 0)
			continue;
		//		std::cout << e <<" " << graph.edges[e].length << ", ";

		//
		// If this edge maps to multiple source edges
		// don't try and join them now (maybe scaffold later
		// but that's hard.)
		//
		if (endToBegin[e].size() != 1)
			continue;

		mapIt = endToBegin[e].begin();

		ssize_t sourceEdge = (*mapIt).first;

		//
		// If multiple source edges map to this end edge,
		// do not try and scaffold it.
		if (beginToEnd[sourceEdge].size() != 1)
			continue;

		// Compute the average begin/end positions.
		(*mapIt).second.beginPos /= (*mapIt).second.count;
		(*mapIt).second.endPos   /= (*mapIt).second.count;
		ssize_t mateSep = (*mapIt).second.beginPos + (*mapIt).second.endPos;
		//UNUSED// ssize_t edgeOverlap;

		//
		// Compute alignment statistics.
		//

		SimpleSequence source, sink;

		source.seq     = (unsigned char*) graph.edges[sourceEdge].seq.seq;
		source.length  = graph.vertices[graph.edges[sourceEdge].src].vertexSize - 1;
		sink.seq       = &(graph.edges[e].seq.seq[graph.edges[e].length - 
																										 graph.vertices[graph.edges[e].dest].vertexSize +1]);
		sink.length = graph.vertices[graph.edges[e].dest].vertexSize - 1;
		//UNUSED// int sinkVertexSize = graph.vertices[graph.edges[e].dest].vertexSize;
		ssize_t maxScore = 0;
		ssize_t nMisMatch, nIndel;
		GrowMatrices(scoreMat, pathMat, source.length+1, sink.length+1);
		ssize_t maxScoreSourcePos, maxScoreSinkPos;
			
		PrefixSuffixAlign(source, sink, maxScoreSourcePos, maxScoreSinkPos, maxScore, nMisMatch, nIndel,
											scoreMat, pathMat, matchMat);
		
		(*mapIt).second.expOverlap = mateSep - mateSepMean;
		(*mapIt).second.beginAlignPos = maxScoreSourcePos;
		(*mapIt).second.endAlignPos = maxScoreSinkPos;
		(*mapIt).second.alignScore  = maxScore;
		++mapIt;
	}

	//
	// Double check mate-pairs for balance.
	//
	for (e = 0 ;e  < graph.edges.size(); e++ ) {
		if (!IsEndEdge(graph,e ))
			continue;

		if (endToBegin[e].size() != 1)
			continue;

		ssize_t startEdge = -1;
		assert(endToBegin[e].size() == 1);
		startEdge = (*endToBegin[e].begin()).first;

		if (beginToEnd[startEdge].size() != 1)
			continue;

		
		// check balance.
		ssize_t balStart = graph.edges[startEdge].balancedEdge;
		ssize_t balEnd   = graph.edges[e].balancedEdge;
		if (endToBegin[balStart].size() == 0) {
			// somehow this isn't linked, remove forward.
			endToBegin[e].clear();
			//			cout << "deleting " << e << endl;
		}
		if ((*endToBegin[balStart].begin()).first != balEnd) {
			// The balance isn't correct.
			// Do not join either e nor e's balance.
			endToBegin[balStart].clear();
			endToBegin[e].clear();
			//			cout << "deleting " << e << " " << balStart << endl;
		}
	}

	cout << "done with mate-pair balance check." << endl;

	//
	// Join the edges that have been merged by mate-pairs.
	//
	ssize_t numJoined = 0;
	for (e = 0; e < graph.edges.size(); e++) {
		if (!IsEndEdge(graph, e))
			continue;
		
		//		cout << "edge: " << e << " out degree: " << endToBegin[e].size() << endl;
		if (endToBegin[e].size() != 1)
			continue;

		// just one mate exists, use that.
		ssize_t mateEdge;
			
		//
		// There should only be 1 sinkToSource edge.
		//
		mateEdge  = (*endToBegin[e].begin()).first;

		//		cout << "edge " << e << " to mate " << mateEdge << " " << beginToEnd[mateEdge].size() << endl;
		if (beginToEnd[mateEdge].size() != 1)
			continue;

		//
		// Check the balance of the mate edges.
		//
		if (endToBegin[graph.edges[mateEdge].balancedEdge].size() != 1) {
			continue;
		}
		if (beginToEnd[graph.edges[e].balancedEdge].size() != 1) {
			continue;
		}
		assert(endToBegin[graph.edges[mateEdge].balancedEdge].size() == 1);
		assert((*endToBegin[graph.edges[mateEdge].balancedEdge].begin()).first ==
					 graph.edges[e].balancedEdge);
		

		// The mate edge should come from a 0/1 vertex (in/out)
		if (!IsBeginEdge(graph, mateEdge))
			continue;
		//
		// There should only be one sinkEdge paired with mateEdge.
		//

		//
		// Join the sinkEdge and the mate edge.
		//
		ssize_t startEdgeSource = graph.edges[mateEdge].src;
		ssize_t startEdgeSourceIndex;
		startEdgeSourceIndex = graph.vertices[startEdgeSource].LookupOutIndex(mateEdge);
		assert(startEdgeSourceIndex != -1);
		graph.vertices[startEdgeSource].out[startEdgeSourceIndex] = -1;
					 
		++numJoined;
			
		// Link end edge forward
		ssize_t endEdgeDest = graph.edges[e].dest;
		graph.vertices[endEdgeDest].AddOutEdge(mateEdge);
			
		// Link start edge back
		graph.edges[mateEdge].src = endEdgeDest;
		vToRemove.push_back(startEdgeSource);
			
		// Adjust the length of the out edge.
		// Use the balanced edge statistics if they exist.
		ssize_t balancedBegin, balancedEnd;
		balancedBegin = graph.edges[mateEdge].balancedEdge;
		balancedEnd   = graph.edges[e].balancedEdge;

		ssize_t alignLength, alignScore, expOvpLength;
			
		//
		// Use only one alignment score for the forward and 
		// reverse strands to make sure there is a 1-1 correspondence
		// between linked edges and their reverse complements.
		//
		if (endToBegin[balancedBegin][balancedEnd].endAlignPos != -1) {
			alignLength  = endToBegin[balancedBegin][balancedEnd].endAlignPos;
			alignScore   = endToBegin[balancedBegin][balancedEnd].alignScore;
			expOvpLength = endToBegin[balancedBegin][balancedEnd].expOverlap;				
		}
		else {
			//
			// 
			alignLength  = endToBegin[e][mateEdge].endAlignPos;
			alignScore   = endToBegin[e][mateEdge].alignScore;
			expOvpLength = endToBegin[e][mateEdge].expOverlap;
				
			// Make sure the balanced edge uses its own statistics.
			endToBegin[e][mateEdge].endAlignPos = -1;
		}

		// Check to see if the alignment is sensical, or if a gap should be inserted.
		if (expOvpLength < 0 and alignLength - alignScore > 5) {
			// There is a gap in sequence coverage. 
			InsertGap(graph.edges, mateEdge, -expOvpLength);
			//			cout << "gapping " << e << " " << mateEdge << " by: " << -expOvpLength << endl;
			/*			cout << graph.edges[e].balancedEdge << " "
							<< graph.edges[mateEdge].balancedEdge << endl;*/
		}
		else {
			//UNUSED//			int sinkAlignLength = graph.vertices[graph.edges[e].dest].vertexSize
			//UNUSED//				- endToBegin[e][mateEdge].beginAlignPos;
			GrowEdge(graph.edges, e, mateEdge, alignLength, 
							 graph.vertices[graph.edges[e].dest].vertexSize - alignLength);

		}
	}
	cout << "joined " << numJoined << " contigs in scaffolds." << endl;
}


void GrowEdge(IntervalGraph &graph,
							ssize_t edge,
							ssize_t length,
							unsigned char *seq) {
	unsigned char *newSeq = new unsigned char[graph.edges[edge].length + length];
	memcpy(newSeq, 	graph.edges[edge].seq.seq, graph.edges[edge].length);
	ssize_t i;
	for (i = 0; i < length; i++) {
		newSeq[graph.edges[edge].length + i] = seq[i];
	}
	delete[] graph.edges[edge].seq.seq;
	graph.edges[edge].seq.seq = newSeq;
	graph.edges[edge].seq.length = graph.edges[edge].length + length;
	graph.edges[edge].length = graph.edges[edge].seq.length;
}
	

void PrefixSuffixAlign(SimpleSequence &sourceSeq,
											 SimpleSequence &sinkSeq,
											 ssize_t &maxScoreSourcePos,
											 ssize_t &maxScoreSinkPos,
											 ssize_t &maxScore, ssize_t &nMisMatch, ssize_t &nIndel,
											 IntMatrix &scoreMat, IntMatrix &pathMat, IntMatrix &matchMat, ssize_t printAlign) {

	
	ssize_t i, j;
	ssize_t ci, cj;

	// 0 - back diagonally
	// 1 - back vertically
	// 2 - back horizontally

	// init the first column to be gaps
	for (i = 0; i < sourceSeq.length + 1; i++) {
		scoreMat[i][0] = -i;
		pathMat[i][0] = 1;
	}
	// init the first row to be free gaps
	for (j = 0; j < sinkSeq.length + 1; j++) {
		scoreMat[0][j] = 0;
		pathMat[0][j] = 2;
	}

	// compute the scores
	ssize_t matchScore, insSinkScore, insSourceScore;
	unsigned char ni, no;
	for (i = 0; i < sourceSeq.length; i++) {
		for (j = 0; j < sinkSeq.length; j++) {
			ci = i + 1; cj = j + 1;
			no = nucToIndex[sourceSeq.seq[i]];
			ni = nucToIndex[sinkSeq.seq[j]];
			/*			if (printAlign) {
							std::cout << i << " " << j << " " << (char) sourceSeq.seq[i] << " " << (ssize_t) no 
							<< " " << (char) sinkSeq.seq[j] << " " << (ssize_t) ni << " " << matchMat[ni][no] << " " 
							<< scoreMat[ci-1][cj-1] << std::endl;
							}
			*/
			matchScore = scoreMat[ci-1][cj-1] + matchMat[ni][no];
			insSourceScore = scoreMat[ci-1][cj] - 1;
			insSinkScore   = scoreMat[ci][cj-1] - 1;
			if (matchScore >= insSourceScore and
					matchScore >= insSinkScore) {
				scoreMat[ci][cj] = matchScore;
				pathMat[ci][cj]  = 0;
			}
			else if (insSourceScore >= matchScore and
							 insSourceScore >= insSinkScore){ 
				scoreMat[ci][cj] = insSourceScore;
				pathMat[ci][cj] = 1;
			}
			else{
				scoreMat[ci][cj] = insSinkScore;
				pathMat[ci][cj] = 2;
			}
		}
	}
	// Find the maximum score by checkign the last column;
	/*
		sourceSeq.PrintSeq(std::cout, "source");
		sinkSeq.PrintSeq(std::cout, "sink");
		PrintMatrix(scoreMat, std::cout,4);
		std::cout << std::endl;
		PrintMatrix(pathMat, std::cout,4);
		std::cout << std::endl;
	*/
	ssize_t lastColumn = sinkSeq.length;
	ssize_t maxRow = 0;
	maxScore = scoreMat[maxRow][lastColumn];
	for (i = 0; i < sourceSeq.length + 1; i++) {
		if (scoreMat[i][lastColumn] > maxScore) {
			maxScore = scoreMat[i][lastColumn];
			maxRow   = i;
		}
	}
	
	ci = maxRow;
	cj = lastColumn;
	nMisMatch = 0;
	nIndel = 0;
	while (ci > 0) {
		assert(cj >= 0);
		if (pathMat[ci][cj] == 0) {if (ci > 0 and cj > 1 and scoreMat[ci-1][cj-1] > scoreMat[ci][cj]) nMisMatch++; ci--; cj--;}
		else if (pathMat[ci][cj] == 1) {ci--; nIndel++;}
		else {cj--; nIndel++;}
	}
	maxScoreSourcePos = maxRow;
	maxScoreSinkPos   = cj;
	if (printAlign) {
		PrintMatrix(scoreMat, std::cout, 3);
	}
}

void GrowEdge(TEdgeList &edges, ssize_t destEdge, ssize_t sourceEdge, ssize_t sourceOffset, ssize_t sourceLength) {
	
	unsigned char* prevSeq;
	prevSeq  = edges[destEdge].seq.seq;
	ssize_t prevLength = edges[destEdge].length;
	edges[destEdge].seq.seq = new unsigned char[prevLength + sourceLength];
	
	assert(prevLength < prevLength + sourceLength);
	memcpy((char*) edges[destEdge].seq.seq, (char*) prevSeq, prevLength);
	ssize_t i;
	for (i = 0; i < sourceLength; i++) { 
		assert(sourceOffset + i < edges[sourceEdge].length );
		edges[destEdge].seq.seq[prevLength + i] = edges[sourceEdge].seq.seq[sourceOffset + i];
	}
	edges[destEdge].seq.length += sourceLength;
	edges[destEdge].length += sourceLength;

	delete[] prevSeq;
}

void InsertGap(TEdgeList &edges,  ssize_t e, ssize_t gapLength) {
	assert(gapLength > 0);
	unsigned char* prevSeq = edges[e].seq.seq;
	ssize_t newLength = edges[e].length + gapLength;
	edges[e].seq.seq = new unsigned char[newLength];

	memcpy(&edges[e].seq.seq[gapLength], prevSeq, edges[e].length);
	ssize_t i;
	for (i = 0; i < gapLength; i++) {
		edges[e].seq.seq[i] = 'N';
	}

	for (i = 0; i < edges[e].intervals->size(); i++) {
		(*edges[e].intervals)[i].edgePos += gapLength;
	}
	edges[e].length += gapLength;
	edges[e].seq.length += gapLength;
	delete [] prevSeq;
}


void InitScaffoldMatchMatrix(IntMatrix &matchMat) {

	CreateMatrix(matchMat, 4, 4);
	
	matchMat[0][0] = 1; matchMat[1][1] = 1; matchMat[2][2] = 1; matchMat[3][3] = 1;
	
	matchMat[0][1] = -1; matchMat[0][2] = -1; matchMat[0][3] = -1;
	matchMat[1][0] = -1; matchMat[1][2] = -1;	matchMat[1][3] = -1;
	matchMat[2][0] = -1; matchMat[2][1] = -1; matchMat[2][3] = -1;
	matchMat[3][0] = -1; matchMat[3][1] = -1;	matchMat[3][2] = -1;
}

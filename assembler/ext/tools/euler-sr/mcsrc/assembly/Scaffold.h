/***************************************************************************
 * Title:          Scaffold.h 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  04/13/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef SCAFFOLD_H_
#define SCAFFOLD_H_
#include "MateLibrary.h"
#include "IntervalGraph.h"

class MateCount {
public:
	ssize_t count;
	ssize_t beginPos;
	ssize_t endPos;
	ssize_t overlap;
	ssize_t expOverlap;
	ssize_t alignScore;
	ssize_t beginAlignPos, endAlignPos;
	MateCount() {
		count = 0;
		beginPos = 0;
		endPos   = 0;
		overlap  = 0;
		expOverlap = 0;
		beginAlignPos = endAlignPos = -1;
		alignScore = 0;
	}
};

ssize_t IsEndEdge(IntervalGraph &g, ssize_t edgeIndex);
ssize_t IsBeginEdge(IntervalGraph &g, ssize_t edgeIndex);
ssize_t IsEndVertex(IntervalGraph &g, ssize_t vertexIndex);
ssize_t IsBeginVertex(IntervalGraph &g, ssize_t vertexIndex);
void GrowMatrices(IntMatrix &scoreMat, IntMatrix &pathMat, ssize_t nrows, ssize_t ncols);
void GrowEdge(TEdgeList &edges, ssize_t destEdge, ssize_t sourceEdge, ssize_t sourceOffset, ssize_t sourceLength);

void GrowEdge(IntervalGraph &graph,
							ssize_t edge,
							ssize_t length,
							unsigned char *seq);


void PrefixSuffixAlign(SimpleSequence &sourceSeq,
											 SimpleSequence &sinkSeq,
											 ssize_t &maxScoreSourcePos,
											 ssize_t &maxScoreSinkPos,
											 ssize_t &maxScore, ssize_t &nMisMatch, ssize_t &nIndel,
											 IntMatrix &scoreMat, IntMatrix &pathMat, IntMatrix &matchMat,
											 ssize_t printAlign=0);


void InsertGap(TEdgeList &edges,  ssize_t e, ssize_t gapLength);


void MatePairScaffoldJoinEdges(IntervalGraph &graph,
															 ReadMateList  &mateList,
															 ssize_t           mateType,
															 ssize_t           minMateCount,
															 IntMatrix     &scoreMat,
															 vector<ssize_t>   &vToRemove);

void InitScaffoldMatchMatrix(IntMatrix &matchMat);

#endif

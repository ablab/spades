/***************************************************************************
 * Title:          graphalign.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  10/01/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef _GRAPH_ALIGN_H_
#define _GRAPH_ALIGN_H_

#include "compatibility.h"
#include "DNASequence.h"
#include <limits.h>
#include <vector>
#include <mctypes.h>

#define FORWARD_DIR 1
#define REVERSE_DIR -1

ssize_t AlignRegion(DNASequence &refSeq, ssize_t refStart, ssize_t refEnd,
									DNASequence &qrySeq, ssize_t qryStart, ssize_t qryEnd,
									IntMatrix &scoreMat, ssize_t gapOpen, ssize_t gapExtend,
									ssize_t &refAlignStart, ssize_t &qryAlignStart, 
									ssize_t *&alignment, ssize_t &length);


class AlignVertex {
public:
  // Use for storing optimal scores/pointers.
  ssize_t matchScore, refGapScore, qryGapScore;
  AlignVertex *matchPrev, *refGapPrev, *qryGapPrev;
  ssize_t matchIndex, refGapIndex, qryGapIndex;
  // Use for storing griddata structure.
  AlignVertex *north, *south, *east, *west, *northwest, *southeast;
  ssize_t scored;
  ssize_t refPos, qryPos;
  AlignVertex(ssize_t refPosP, ssize_t qryPosP) {
    refPos     = refPosP;
    qryPos     = qryPosP;
    north      = NULL;
    south      = NULL;
    northwest  = NULL;
    east       = NULL;
    west       = NULL;
    southeast  = NULL;
    matchPrev  = NULL;
    refGapPrev = NULL;
    qryGapPrev = NULL;
    matchIndex = none; 
    refGapIndex= none;
    qryGapIndex= none;
    matchScore = infty;
    refGapScore= infty;
    qryGapScore= infty;
    scored = 0;
  }
  AlignVertex* GetNorth();
  AlignVertex* GetWest();
  AlignVertex* GetNorthWest();
  AlignVertex* GetSoutheast();
  AlignVertex* GetSouth();
  AlignVertex* GetEast();

  ssize_t ComputeDiagonalScore(AlignVertex *diag, AlignVertex *afRef, AlignVertex *afQry, 
			     IntMatrix &scoreMat, char refChar, char qryChar);

  ssize_t ComputeAffineGapScore(AlignVertex *afGap, AlignVertex *diag, 
			      ssize_t affineGapOpen,
			      ssize_t affineGapExtend);

  ssize_t ScoreQryGap(ssize_t gapOpenCost, ssize_t gapExtendCost, ssize_t maxScore, ssize_t minScore);
  ssize_t ScoreRefGap(ssize_t gapOpenCost, ssize_t gapExtendCost, ssize_t maxScore, ssize_t minScore);
  ssize_t ScoreMatch(DNASequence &refSeq, DNASequence &qrySeq, 
		   ssize_t refPos, ssize_t qryPos,
		   ssize_t refStartPos, ssize_t qryStartPos, ssize_t extDir,
		   IntMatrix &scoreMat, ssize_t maxScore, ssize_t minScore);

  // use these to ensure a path is not taken
  void StopMatchAlign() { matchScore = INT_MAX; }
  void StopRefGapAlign() { refGapScore = INT_MAX; }
  void StopQryGapAlign() { qryGapScore = INT_MAX; }

  AlignVertex* InitializeSouthEast(AlignVertex *southeast);
  AlignVertex* InitializeSouth(AlignVertex *south);
  AlignVertex* InitializeEast(AlignVertex *east);

  ssize_t Perfect(ssize_t maxScore);
  ssize_t Capable(ssize_t maxScore);
  ssize_t ValidRefGap(ssize_t gapOpenScore, ssize_t maxScore) {
    return refGapScore < maxScore || gapOpenScore + matchScore < maxScore;
  }

  ssize_t ValidQryGap(ssize_t gapOpenScore, ssize_t maxScore) {
    return qryGapScore < maxScore || gapOpenScore + matchScore < maxScore;
  }

  ssize_t ValidMatch(ssize_t maxScore) {
    return matchScore < maxScore;
  }
  
  void Score(DNASequence &refSeq, DNASequence &qrySeq, 
	     ssize_t refPos, ssize_t qryPos, 
	     ssize_t refStartPos, ssize_t qryStartPos, ssize_t extDir,
	     IntMatrix &scoreMat, ssize_t gapOpen, ssize_t gapExtend, ssize_t maxScore, ssize_t minScore) {
    ScoreRefGap(gapOpen, gapExtend, maxScore, minScore);
    ScoreQryGap(gapOpen, gapExtend, maxScore, minScore);
    ScoreMatch(refSeq, qrySeq, 
	       refPos, qryPos, 
	       refStartPos, qryStartPos, extDir,
	       scoreMat, maxScore, minScore);
    scored = 1;
  }

  static _INT_ match, refGap, qryGap, none;
  static _INT_ infty;
};


#endif

/***************************************************************************
 * Title:          SAP.h 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  11/26/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef SAP_H_
#define SAP_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include "DNASequence.h"
#include "IntegralTuple.h"
#include "SeqReader.h"
#include "SeqUtils.h"
#include "utils.h"
#include "mctypes.h"
#include "FixErrorsStats.h"

using namespace std;


class SEdge {
public:
	// prevNuc, prevLevel, and prevPosition serve as the
	// back pointers in the dynamic programming table (3-d grid).
  char prevNuc;
  unsigned char prevLevel;
	unsigned char solidStretch;
	
  ssize_t prevPosition;
	
	// The current tuple
	CountedIntegralTuple tuple;

	// The score of this point in the lattice.
  ssize_t score;
  static int tupleSize;
  SEdge() {
    prevNuc = 0;
    prevLevel = 0;
    prevPosition = 0;
		//    tuple = MultTuple("");
    score = 0;
		solidStretch = 0;
  }
  SEdge& operator=(const SEdge& rhs) {
		if (this != &rhs) {
			prevNuc   = rhs.prevNuc;
			prevLevel = rhs.prevLevel;
			prevPosition = rhs.prevPosition;
			tuple = rhs.tuple;
			score = rhs.score;
			solidStretch = rhs.solidStretch;
		}
    return *this;
  }

};

class FixParams {
public:
  ssize_t maxGap;
  ssize_t gapOpen;
  ssize_t gapExtend;
  ssize_t misMatch;
  ssize_t scoreThreshold;
  ssize_t span;
  ssize_t maxTrim;
  ssize_t edgeLimit;
	ssize_t startScore;
	ssize_t stepScore;
	ssize_t maxScore;
};

typedef std::map<CountedIntegralTuple, SEdge> Cell;
typedef Cell::iterator EdgeIterator;


typedef Cell* Column;
typedef Column* Grid;
typedef std::vector<Grid> Cube;
void ReverseSeq(DNASequence &seq, DNASequence &rev);
void PatchSeq(DNASequence &seq, ssize_t pos, DNASequence &patch, ssize_t replaceLength);
void CreateGrid(ssize_t dim, Grid& grid);
void DeleteGrid(ssize_t dim, Grid &grid);

ssize_t FindSolidPosition(DNASequence &seq, 
											CountedIntegralTupleDict &spectrum,
											ssize_t span, ssize_t &lastSolidPos, ssize_t& pos);

ssize_t StoreEdge(ssize_t prevPos, ssize_t prevLevel, 
							ssize_t prevNuc, CountedIntegralTuple prevTuple, 
							Cell& cell, CountedIntegralTuple tuple, ssize_t score, 
							ssize_t prevScore, ssize_t prevSolid);

ssize_t FindMinimumScoreEdge(Cube &matrix, FixParams &p, EdgeIterator &edgeIt, 
												 ssize_t &minN, ssize_t &minK, ssize_t& numMin, ssize_t &minScore, 
												 CountedIntegralTuple &minTuple, ssize_t &minTrim, ssize_t &trim);

void Backtrack(Cube &matrix, ssize_t pos, ssize_t level, ssize_t nuc, 
							 CountedIntegralTuple &tuple, 
							 DNASequence &seq, Stats &fixStats, FixParams &params,
							 ssize_t &firstEdit, ssize_t &end);

ssize_t SolidifyRead(DNASequence &seq,
								 //								 CountedIntegralTupleList &spectrum,
								 CountedIntegralTupleDict &spectrum,
								 IntMatrix &scoreMat,
								 FixParams &params, Stats &fixStats, ssize_t &readWasModified);

ssize_t SolidifyUntilFixed(DNASequence &seq,
											 CountedIntegralTupleList &origEnumSeq, ssize_t origEnumPos,
											 //											 CountedIntegralTupleList &spectrum,
											 CountedIntegralTupleDict &spectrum,
											 IntMatrix &scoreMat,
											 int tupleSize, FixParams &params,
											 DNASequence &fixedSeq, ssize_t &replaceLength, 
											 Stats &fixStats,
											 ssize_t &lastEdit, ssize_t &end);

#endif

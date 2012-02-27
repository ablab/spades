/***************************************************************************
 * Title:          FixErrorsLibrary.h 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  11/26/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef FIX_ERRORS_H_
#define FIX_ERRORS_H_

#include "DNASequence.h"
#include "RuleList.h"
#include "SeqReader.h"
#include "SeqUtils.h"
#include "utils.h"
#include "hash/HashUtils.h"
#include "IntegralTuple.h"


#include <vector>
#include <iostream>
#include <ext/hash_map>
#include <map>


class Stats {
public:
  ssize_t numEdge;
  ssize_t numIns;
  ssize_t numDel;
  ssize_t numMut;
  ssize_t numNoSolid;
  ssize_t numNoPathFound;
  ssize_t numMultiplePaths;
	ssize_t numErrorAtEnd;
  Stats() {
    Reset();
  }
  void Reset() {
    numEdge = numIns = numDel = numMut = 0;
    numNoSolid = 0;
    numNoPathFound = 0;
    numMultiplePaths = 0;
		numErrorAtEnd = 0;
  }
	Stats &Append(Stats &s) {
    numEdge += s.numEdge;
    numIns += s.numIns;
    numDel += s.numDel;
    numMut += s.numMut;
    numNoSolid += s.numNoSolid;
    numNoPathFound  += s.numNoPathFound;
    numMultiplePaths += s.numMultiplePaths;
		numErrorAtEnd += s.numErrorAtEnd;
    return *this;
	}
  Stats &operator+=(Stats &s) {
		Append(s);
  }
	friend std::ostream &operator<<(std::ostream &strm, const Stats &s) {
		strm << s.numEdge <<" " << s.numIns 
				 <<" " << s.numDel << " " 
				 << s.numMut <<" " 
				 << s.numNoSolid << " " << s.numNoPathFound 
				 << " " << s.numMultiplePaths << std::endl;
	}
};


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

int SEdge::tupleSize = 0;
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

ssize_t FindSolidPosition(DNASequence &seq, CountedIntegralTupleList &spectrum,
											ssize_t span, ssize_t &lastSolidPos, ssize_t& pos);

ssize_t StoreEdge(ssize_t prevPos, ssize_t prevLevel, 
							ssize_t prevNuc, CountedIntegralTuple prevTuple, 
							Cell& cell, CountedIntegralTuple tuple, ssize_t score, 
							ssize_t prevScore, ssize_t prevSolid);

ssize_t FindMinimumScoreEdge(Cube &matrix, FixParams &p, EdgeIterator &edgeIt, 
												 ssize_t &minN, ssize_t &minK, ssize_t& numMin, ssize_t &minScore, 
												 CountedIntegralTuple &minTuple, ssize_t &minTrim, ssize_t &trim);

ssize_t Backtrack(Cube &matrix, ssize_t pos, ssize_t level, ssize_t nuc, CountedIntegralTuple &tuple, 
							DNASequence &seq, Stats &fixStats, FixParams &params,
							ssize_t &firstEdit, ssize_t &end);

ssize_t SolidifyRead(DNASequence &seq,
								 CountedIntegralTupleList &spectrum,
								 FixParams &params, Stats &fixStats, ssize_t &readWasModified);

ssize_t SolidifyUntilFixed(DNASequence &seq,
											 CountedIntegralTupleList &origEnumSeq, ssize_t origEnumPos,
											 CountedIntegralTupleList &spectrum,
											 int tupleSize, FixParams &params,
											 DNASequence &fixedSeq, ssize_t &replaceLength, 
											 Stats &fixStats,
											 ssize_t &lastEdit, ssize_t &end);



#endif FIX_ERRORS_H_

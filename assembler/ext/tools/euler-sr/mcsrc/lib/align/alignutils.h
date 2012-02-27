/***************************************************************************
 * Title:          alignutils.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  03/17/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef _ALIGNUTILS_
#define _ALIGNUTILS_

#include "DNASequence.h"
#include "mctypes.h"
#include <iomanip>
#include <iostream>
#include <vector>

// Constants for storing the direction of optimal paths (not 
// used for computing the scores of paths!).

class Score {
public:
  Score(std::string &scoreMatName, ssize_t open, ssize_t extend);
  Score();
	Score(ssize_t match, ssize_t mismatch, ssize_t gapOpen, ssize_t gapExtend);
	void CalculateEditDistance();
  IntMatrix scoreMat;
  ssize_t gapOpen, gapExtend;
};


static ssize_t MATCH  = 0;
static ssize_t GAP_A  = 1;
static ssize_t GAP_B  = 2;

static ssize_t GAP_EXTEND = 0;
static ssize_t GAP_OPEN   = 1;
static ssize_t CLOSE_GAP_A = 3;
static ssize_t CLOSE_GAP_B = 4;

static ssize_t BAND_MATCH  = 5;
static ssize_t BAND_GAP_A  = 6;
static ssize_t BAND_GAP_B  = 7;

static ssize_t LOCAL_START = 8;

static ssize_t DO_LOCAL_ALIGN = 1;
static ssize_t DO_GLOBAL_ALIGN = 0;
// Infinity enough!
static ssize_t INF = 9999999;

static char nuc_index[256] = {
  0,1,2,3,4,4,4,4,4,4,   // 0 9
  4,4,4,4,4,4,4,4,4,4,   // 10 19
  4,4,4,4,4,4,4,4,4,4,   // 20 29
  4,4,4,4,4,4,4,4,4,4,   // 30 39
  4,4,4,4,4,4,4,4,4,4,   // 40 49
  4,4,4,4,4,4,4,4,4,4,   // 50 59
  4,4,4,4,4,1,4,2,4,4,   // 60 69
  4,0,4,4,4,4,4,4,4,4,   // 70 79
  4,4,4,4,3,4,4,4,4,4,   // 80 89
  4,4,4,4,4,4,4,1,4,2,   // 90 99
  4,4,4,0,4,4,4,4,4,4,   // 100 109
  4,4,4,4,4,4,3,4,4,4,   // 110 119
  4,4,4,4,4,4,4,4,4,4,   // 120 129
  4,4,4,4,4,4,4,4,4,4,   // 130 139
  4,4,4,4,4,4,4,4,4,4,   // 140 149
  4,4,4,4,4,4,4,4,4,4,   // 150 159
  4,4,4,4,4,4,4,4,4,4,   // 160 169
  4,4,4,4,4,4,4,4,4,4,   // 170 179
  4,4,4,4,4,4,4,4,4,4,   // 180 189
  4,4,4,4,4,4,4,4,4,4,   // 190 199
  4,4,4,4,4,4,4,4,4,4,   // 200 209
  4,4,4,4,4,4,4,4,4,4,   // 210 219
  4,4,4,4,4,4,4,4,4,4,   // 220 229
  4,4,4,4,4,4,4,4,4,4,   // 230 239
  4,4,4,4,3,2,1,0,3,2,   // 240 249
  1,0,3,2,1,0};          // 250 255



/*
  static way of initializing this:
  
  ssize_t a, c, t, g;
  a => 0; c => 1; g => 2; t => 3;

  nuc_index[0] = g;
  nuc_index[1] = a;
  nuc_index[2] = c;
  nuc_index[3] = t;
  nuc_index[(unsigned char)-1] = 2;  255
  nuc_index[(unsigned char)-2] = 0;  254
  nuc_index[(unsigned char)-3] = 1;  253
  nuc_index[(unsigned char)-4] = 3;  252
  nuc_index[(unsigned char)-5] = 2;  251
  nuc_index[(unsigned char)-6] = 0;  250
  nuc_index[(unsigned char)-7] = 1;  249
  nuc_index[(unsigned char)-8] = 3;  248
  nuc_index[(unsigned char)-9] = 2;  247
  nuc_index[(unsigned char)-10] = 0; 246
  nuc_index[(unsigned char)-11] = 1; 245
  nuc_index[(unsigned char)-12] = 3; 244

  nuc_index['g'] = 0; 103
  nuc_index['G'] = 0;  71
  nuc_index['a'] = 1;  97
  nuc_index['A'] = 1;  65
  nuc_index['c'] = 2;  99
  nuc_index['C'] = 2;  67
  nuc_index['t'] = 3; 116
  nuc_index['T'] = 3;  84
  nuc_index['n'] = 4; 110
  nuc_index['N'] = 4;  78
  nuc_index['x'] = 4; 120
  nuc_index['X'] = 4;  88
*/


void ResizeScoreMat(IntMatrix &scoreMat);

ssize_t CalcNumMatches(DNASequence &seqa, DNASequence &seqb, ssize_t *locations);

ssize_t CalcPercentIdentity(DNASequence &seqa, DNASequence &seqb, 
			  ssize_t *locations);

ssize_t CalcRepatPercentIdentity(DNASequence &seqa, DNASequence &seqb, 
			       ssize_t *locations);

ssize_t CalcNonRepeatPercentIdentity(DNASequence &seqa, DNASequence &seqb, 
				   ssize_t *locations);

ssize_t AffineAlign(DNASequence &seqa, DNASequence &seqb, 
		  ssize_t match, ssize_t mismatch, ssize_t gap, 
		  ssize_t gapOpen, ssize_t gapExtend,
		  ssize_t *&locations, IntMatrix &scoreMat);

ssize_t AffineAlign(DNASequence &seqa, DNASequence &seqb, 
		  ssize_t match, ssize_t mismatch, ssize_t gap, 
		  ssize_t gapOpen, ssize_t gapExtend,
		  ssize_t *&locations);


ssize_t FitAlign(DNASequence &seqa, DNASequence &seqb, 
							 ssize_t match, ssize_t mismatch, ssize_t gap, 
							 ssize_t *&locations, IntMatrix &scorMat,
							 IntMatrix &score, IntMatrix &path);

ssize_t Align(DNASequence &seqa, DNASequence &seqb, 
	    ssize_t match, ssize_t mismatch, ssize_t gap, 
	    ssize_t *&locations, IntMatrix & scoreMat);

ssize_t BandedAlign(DNASequence &seqa, DNASequence &seqb, 
									ssize_t match, ssize_t mismatch, ssize_t gap, 
									ssize_t k,
									ssize_t *&locations,  
									IntMatrix &score,
									IntMatrix &path,
									IntMatrix &scoreMat,
									ssize_t *optScores = NULL);


ssize_t BandedAlign(DNASequence &seqa, DNASequence &seqb, 
									ssize_t match, ssize_t mismatch, ssize_t gap, 
									ssize_t k,
									ssize_t *&locations, ssize_t *optScores=NULL);

ssize_t CoreAlign(DNASequence &seqa, DNASequence &seqb,
		ssize_t match, ssize_t mismatch, ssize_t gap, 
		ssize_t *&locations, IntMatrix &score, IntMatrix &path,
		ssize_t alternative, // 0 for smith-waterman, +infty for needleman wunsch
		ssize_t &globalMinRow, ssize_t &globalMinCol);

ssize_t LocalAlign(DNASequence &seqa, DNASequence &seqb, 
		 ssize_t match, ssize_t mismatch, ssize_t gap, 
		 ssize_t *&locations,
		 IntMatrix &scorMat);


ssize_t LocalAlign(DNASequence &seqa, DNASequence &seqb, 
		 Score &scoreMat, ssize_t *&locations);

ssize_t OverlapAlign(DNASequence &seqa, DNASequence &seqb,
									 ssize_t match, ssize_t mismatch, ssize_t gap, 
									 ssize_t *&locations, IntMatrix &scoreMat);

ssize_t AffineLocalAlign(DNASequence &seqa, DNASequence &seqb, 
		       ssize_t match, ssize_t mismatch, ssize_t gap, 
		       ssize_t gapOpen, ssize_t gapExtend,
		       ssize_t *&locations, IntMatrix &scoreMat);

ssize_t AffineLocalAlign(DNASequence &seqa, DNASequence &seqb,
		       Score &scoreMat, ssize_t *&locations);

void ParseScoreMatStr(std::string scoreMatStr, IntMatrix &scoreMat);

void ReadScoreMatFile(std::string scoreMatFileName, IntMatrix & scoreMat);

class AlignVertex;

void GetGraphAlignment(AlignVertex *vertex, 
		       ssize_t *&locations, ssize_t &length, 
		       ssize_t refStartPos, ssize_t qryStartPos, ssize_t extDir);


ssize_t ScoreBandedAffineAlign(DNASequence &refSeq, DNASequence &qrySeq, // sequences to compare 
			    ssize_t initialScore,  // score of seed
			    IntMatrix &scoreMat, ssize_t gapOpen, ssize_t gapExtend, // scoring parameters
			    ssize_t maxScore,  // negative score is good.
			    ssize_t *&locations, ssize_t &length, 
			    ssize_t &refAlignStart, ssize_t &qryAlignStart, // used for storing the resulting map
			    ssize_t refStartPos=0, ssize_t qryStartPos=0, 
			    ssize_t extDir=1
			    );


void AssignScoreMat(IntMatrix &sm1, IntMatrix &sm2);


void DeleteGraph(std::vector<AlignVertex*> &graph);

void InitScoreMat(IntMatrix &scoreMat, ssize_t match, ssize_t mismatch);

template<typename T>
T** CreateMatrix(ssize_t rows, ssize_t cols, ssize_t init=0) {
  ssize_t i, j;
  T **mat;
  mat = new T*[rows];
  for (i = 0; i < rows; i++) {
    mat[i] = new T[cols];
    for (j = 0; j < cols; j++) {
      mat[i][j] = init;
    }
  }
  return mat;
}


template <typename T>
void DeleteMatrix(T **mat, ssize_t rows) {
  ssize_t r;
  for (r = 0; r < rows; r++)
    delete [] mat[r];
  delete[] mat;
}


template <typename T>
void PrintMatrix(T **mat, ssize_t rows, ssize_t cols, ssize_t width=3) {
  ssize_t r, c;
  std::cout << std::setw(width) << 1;
  for (c = 0; c < cols; c++) 
    std::cout << std::setw(width) << c;
  std::cout << std::endl;
  for (r = 0; r < rows; r++) {
    std::cout << std::setw(width) << r;
    for (c = 0; c < cols; c++ ) {
      std::cout << std::setw(width) << mat[r][c];
    }
    std::cout << std::endl;
  }
}



#endif

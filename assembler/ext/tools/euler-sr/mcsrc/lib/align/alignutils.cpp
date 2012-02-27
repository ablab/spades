/***************************************************************************
 * Title:          alignutils.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <algorithm>
#include <iostream>
#include <istream>
#include <sstream>

#include "mctypes.h"
#include "DNASequence.h"
#include "alignutils.h"
#include "compatibility.h"

using namespace std;
Score::Score(ssize_t match, ssize_t mismatch, ssize_t gapOpenP, ssize_t gapExtendP) {
	CreateMatrix(scoreMat, 5, 5);
	ssize_t i, j;
	for ( i= 0; i < 5; i++) {
		for (j = 0; j < 5; j++) {
			scoreMat[i][j] = mismatch;
		}
		scoreMat[i][i] = match;
	}
	
	gapOpen = gapOpenP;
	gapExtend = gapExtendP;
}

Score::Score() {
  /* 
     Initialize score matrix with UCSC human vs mouse calibrated 
     matrix.
  */
  char *home = getenv("HOME");
  std::string scoreMatName = std::string(home) + "/projects/mcsrc/lib/align/data/scoremat.txt";
  ReadScoreMatFile(scoreMatName, scoreMat);
  gapOpen = 400;
  gapExtend = 30;
}

void Score::CalculateEditDistance() {
	ssize_t i, j;
	ResizeScoreMat(scoreMat);
	// Use a score matrix that has the simple interpretation
	// of score = edit distance.
	for (i = 0; i < 5; i++) {
		for (j = 0; j < 5; j++ ){
			if (i == j) {
				scoreMat[i][j] = 0;
			}
			else {
				scoreMat[i][j] = 1;
			}
		}
	}
	gapOpen   = 1;
	gapExtend = 1;
}

Score::Score(std::string &scoreMatName, ssize_t open, ssize_t extend) {
  ReadScoreMatFile(scoreMatName, scoreMat);
  gapOpen = open;
  gapExtend = extend;
}

ssize_t AffineAlign(DNASequence &seqa, DNASequence &seqb, 
									ssize_t match, ssize_t mismatch, ssize_t gap, 
									ssize_t gapOpen, ssize_t gapExtend,
									ssize_t *&locations) {
  IntMatrix empty;
  return AffineAlign(seqa, seqb, match, mismatch, gap, gapOpen, gapExtend,
										 locations, empty);
}

ssize_t BandedAlign(DNASequence &seqa, DNASequence &seqb, 
									ssize_t match, ssize_t mismatch, ssize_t gap, 
									ssize_t k,
									ssize_t *&locations, ssize_t *optScores) {
  IntMatrix empty;
	IntMatrix emptyScores;
	IntMatrix emptyScoreMat;
  return BandedAlign(seqa, seqb, match, mismatch, gap, k, locations, emptyScores, empty, emptyScoreMat, optScores);
}

ssize_t RepeatMasked(char n) {
  return (n == 'a' || n == 'c' || n == 'g' || n == 't');
}

ssize_t CalcRepatPercentIdentity(DNASequence &seqa, DNASequence &seqb, ssize_t *locations) {
  ssize_t i;
  ssize_t length = seqa.length;
  ssize_t front, back;

  // find the front
  i = 0;
  while (i < length && locations[i] == -1)
    i++;
  front = i;

  // find the back
  i = length-1;
  while (i >= front && locations[i] == -1) 
    i--;
  back = i+1;

  if (front >= back)
    return 0;

  ssize_t countEqual = 0;
  ssize_t numRepeats = 0;
  for (i = front; i < back; i++) {
    if (locations[i] != -1) {
      if (RepeatMasked(seqa.seq[i]) && RepeatMasked(seqb.seq[i])) {
				numRepeats++;
				if (nuc_index[(unsigned char)seqa.seq[i]] == nuc_index[(unsigned char)seqb.seq[locations[i]]])
					++countEqual;
      }
    }
  }
  if (numRepeats > 0)
    return (ssize_t) (double(countEqual)/numRepeats);
  else 
    return 0;
}

ssize_t CalcNonRepeatPercentIdentity(DNASequence &seqa, DNASequence &seqb, ssize_t *locations) {
  ssize_t i;
  ssize_t length = seqa.length;
  ssize_t front, back;

  // find the front
  i = 0;
  while (i < length && locations[i] == -1)
    i++;
  front = i;

  // find the back
  i = length-1;
  while (i >= front && locations[i] == -1) 
    i--;
  back = i+1;

  if (front >= back)
    return 0;

  ssize_t countEqual = 0;
  ssize_t numNonRepeats = 0;
  for (i = front; i < back; i++) {
    if (locations[i] != -1) {
      if (!RepeatMasked(seqa.seq[i]) && !RepeatMasked(seqb.seq[i])) {
				numNonRepeats++;
				if (nuc_index[(unsigned char)seqa.seq[i]] == nuc_index[(unsigned char)seqb.seq[locations[i]]])
					++countEqual;
      }
    }
  }
  if (numNonRepeats > 0)
    return (ssize_t) (double(countEqual)/numNonRepeats);
  else 
    return 0;
}

ssize_t CalcNumMatches(DNASequence &seqa, DNASequence &seqb, ssize_t *locations) {
	ssize_t i, nMatches;
	nMatches = 0;
  for (i = 0; i < seqa.length; i++ ) {
    if (locations[i] != -1) {
      if (nuc_index[(unsigned char)seqa.seq[i]] == nuc_index[(unsigned char)seqb.seq[locations[i]]])
				nMatches++;
    }
  }
	return nMatches;
}


ssize_t CalcPercentIdentity(DNASequence &seqa, DNASequence &seqb, ssize_t *locations) {
  ssize_t i;
  ssize_t length = seqa.length;
  ssize_t front, back;
  // find the front
  i = 0;
  ssize_t equal = 0;
  for (i = 0; i < seqa.length; i++ ) {
    if (locations[i] != -1) {
      if (nuc_index[(unsigned char)seqa.seq[i]] == nuc_index[(unsigned char)seqb.seq[locations[i]]])
				equal+=2;
    }
  }
  return (ssize_t) (double(equal) / (seqa.length + seqb.length));
  
  while (i < length && locations[i] == -1)
    i++;

  front = i;
  i = length-1;
  while (i >= front && locations[i] == -1) 
    i--;
  
  if (front >= back)
    return 0;

  back = i+1;
  ssize_t countEqual = 0;
  for (i = front; i < back; i++) {
    if (locations[i] != -1)
      if (nuc_index[(unsigned char)seqa.seq[i]] == nuc_index[(unsigned char)seqb.seq[locations[i]]])
				++countEqual;
  }
  
  return (ssize_t)  (double(countEqual)/(back-front));
}

void AssignScoreMat(IntMatrix &sm1, IntMatrix &sm2) {
  sm1 = sm2;
}


void InitScoreMat(IntMatrix &scoreMat, ssize_t match, ssize_t mismatch) {
  ssize_t i, j;
  CreateMatrix(scoreMat, 5,5);
  for (i = 0; i < scoreMat.size(); i++) {
    for (j = 0; j < scoreMat.size(); j++) {
      scoreMat[i][j] = mismatch;
    }
    scoreMat[i][i] = match;
  }
  scoreMat[4][4] = mismatch; // N and N should not be rewarded
}

void ReadScoreMatFile(std::string scoreMatFileName, 
											IntMatrix &scoreMat) {
  std::ifstream in;
  in.open(scoreMatFileName.c_str());
  if (!in.good()) {
    std::cout << "could not open score mat file: " 
							<< scoreMatFileName << std::endl;
    exit(0);
  }
  CreateMatrix(scoreMat, 5,5);

  std::string scoreMatStr;
  std::string tmp;
  std::getline(in, tmp);
  while (tmp.length() > 0 and tmp.c_str()[0] == '#') {
    std::getline(in, tmp);
  }
  if (tmp.length() > 0) {
    scoreMatStr += tmp + " ";
  }

  while (in.good()) {
    in >> tmp;
    scoreMatStr += tmp + " ";
  }
  ParseScoreMatStr(scoreMatStr, scoreMat); 
}

void ResizeScoreMat(IntMatrix &scoreMat) {
  CreateMatrix(scoreMat, 5, 5);
}

void ParseScoreMatStr(std::string scoreMatStr, IntMatrix & scoreMat) {
  std::stringstream matstream(scoreMatStr);
  ssize_t i, j;
  //UNUSED// ssize_t v;
  for (i = 0; i < 5; i++) 
    for (j = 0; j < 5; j++) {
      matstream >> scoreMat[i][j];
    }
}

void GetAlignedLocations(IntMatrix &matchPath, IntMatrix &gapAPath, IntMatrix &gapBPath, 
												 ssize_t rows, ssize_t cols, ssize_t *locations, ssize_t localAlign,
												 ssize_t initialRow=-1, ssize_t initialCol=-1);

void init(ssize_t match, ssize_t mismatch, IntMatrix &scoreMat) {
  // Have mapping for binary encoding of nucleotides that is used 
  // by DNASequence.h, and by SeqReader.h
  CreateMatrix(scoreMat, 6,6);
  scoreMat[0][0] = match; 
  scoreMat[0][1] = mismatch; 
  scoreMat[0][2] = mismatch; 
  scoreMat[0][3] = mismatch; 
  scoreMat[0][4] = mismatch; 
  scoreMat[1][0] = mismatch; 
  scoreMat[1][1] = match; 
  scoreMat[1][2] = mismatch; 
  scoreMat[1][3] = mismatch; 
  scoreMat[1][4] = mismatch; 
  scoreMat[2][0] = mismatch; 
  scoreMat[2][1] = mismatch;
  scoreMat[2][2] = match;    
  scoreMat[2][3] = mismatch; 
  scoreMat[2][4] = mismatch; 
  scoreMat[3][0] = mismatch; 
  scoreMat[3][1] =mismatch; 
  scoreMat[3][2] = mismatch; 
  scoreMat[3][3] = match; 
  scoreMat[3][4] = mismatch;   
  scoreMat[4][0] = mismatch; 
  scoreMat[4][1] = mismatch; 
  scoreMat[4][2] = mismatch; 
  scoreMat[4][3] = mismatch; 
  scoreMat[4][4] = match;   
}

void GetBandedAlignedLocations(IntMatrix &alignScores,
															 IntMatrix &matchPath, 
															 ssize_t rows, ssize_t cols, 
															 ssize_t *locations,
															 ssize_t *scores,
															 ssize_t k,
															 ssize_t initialRow=-1, ssize_t initialCol=-1) {

  ssize_t r, c;
  r = initialRow;
  c = initialCol;
  while (r > 0 || c > k+1) {
    if (matchPath[r][c] == BAND_MATCH) {
      r = r -1;
      locations[r] = c + r - k - 1; // reverse transform from 
			// alignment routine.
			if (scores != NULL)
				scores[r] = alignScores[r][c];
    }
    else if (matchPath[r][c] == BAND_GAP_A) {
      r = r - 1;
      c = c + 1;
    }
    else if (matchPath[r][c] == BAND_GAP_B) {
      c = c - 1;
    }
		else {
			cout << "match path pos: " << r << " " << c << " is not right: " << matchPath[r][c] << " should be: " << BAND_MATCH << ", " << BAND_GAP_A<< ", or " << BAND_GAP_B << endl;
			PrintMatrix(matchPath, cout);
			assert(0);
		}
  }
}

void GetAlignedLocations(IntMatrix &matchPath, // path computed from dyn prog mat
												 IntMatrix &gapAPath,  // affine gap mat
												 IntMatrix &gapBPath,  // ditto
												 ssize_t rows, ssize_t cols, // size of mats
												 ssize_t *locations, // store locations of b that a maps to here.
												 ssize_t localAlign, // trace path of a local align
												 ssize_t initialRow, ssize_t initialCol) {
  
  ssize_t row, col;
  // Matrix references which grid we are in. 0 = match, 1 = gap A, 2 = gapB
  ssize_t matrix; 
  if (initialRow == -1)
    row = rows;
  else
    row = initialRow;

  if (initialCol == -1)
    col = cols;
  else
    col = initialCol;

  matrix = 0; // Start out on match.
  while ((row > 0 || col > 0) && (!localAlign || (matchPath[row][col] != LOCAL_START))) {
    if (matrix == 0) {
			// Traversing match matrix
      if (localAlign && matchPath[row][col] == LOCAL_START) {
				row = 0;
				col = 0;
      }
      else if (matchPath[row][col] == MATCH) {
				assert(row-1>=0);
				locations[row-1] = col-1;
				row = row - 1;
				col = col - 1;
      }
      else if (matchPath[row][col] == GAP_A) {
				row = row - 1;
      }
      else if (matchPath[row][col] == GAP_B) {
				col = col - 1;
      }
      else if (matchPath[row][col] == CLOSE_GAP_A) {
				matrix = 2;
      }
      else if (matchPath[row][col] == CLOSE_GAP_B) {
				matrix = 1;
      }
    }
    else if (matrix == 1) {
      // Traversing gap in B matrix
      if (gapBPath[row-1][col-1] == GAP_EXTEND) {
				col = col - 1;
      }
      else if (gapBPath[row-1][col-1] == GAP_OPEN) {
				matrix = 0;
				col = col - 1;
      }
    }
    else if (matrix == 2) {
      if (gapAPath[row-1][col-1] == GAP_EXTEND) {

				row = row - 1;
      }
      else if (gapAPath[row-1][col-1] == GAP_OPEN) {
				row = row - 1;
				matrix = 0;
      }
    }
    //    std::cout << std::endl;
  }
}

ssize_t AffineAlign(DNASequence &seqa, DNASequence &seqb, 
									ssize_t match, ssize_t mismatch, ssize_t gap, 
									ssize_t gapOpen, ssize_t gapExtend,
									ssize_t *&locations, IntMatrix & scoreMat ){

  if (locations == NULL) {
    locations = new ssize_t[seqa.length];
  }
  //  init(match, mismatch);
  if (scoreMat.size() == 0)
    InitScoreMat(scoreMat, match, mismatch);

  ssize_t rows = seqa.length;
  ssize_t cols = seqb.length;

  IntMatrix gapAScore, gapBScore, matchScore;
  IntMatrix gapAPath, gapBPath, matchPath;
  ssize_t matchScoreVal, gapAScoreVal, gapBScoreVal, 
    gapACloseScore, gapBCloseScore, gapOpenScore, gapExtendScore,
    gapScore;


  CreateMatrix(gapAScore, rows+1, cols+ 1);
  CreateMatrix(gapBScore, rows+1, cols+ 1);
  CreateMatrix(matchScore, rows+1, cols+ 1);
  
  CreateMatrix(gapAPath, rows+1, cols+ 1);
  CreateMatrix(gapBPath, rows+1, cols+ 1);
  CreateMatrix(matchPath, rows+1, cols+ 1);


  ssize_t r, c;
  ssize_t rdx, cdx;
  // Initialize the first row/column of all matricies.
  gapAScore[0][0]  = INF;
  gapAPath[0][0] = 0;
  gapBScore[0][0]  = INF;
  gapBPath[0][0] = 0;

  
  matchScore[0][0] = 0;
  ssize_t gapClose, minScore;
  
  // Initialize the first row  
  for (c = 0; c < cols; c++) {
    cdx = c + 1;
    gapExtendScore = gapBScore[0][cdx-1] + gapExtend;
    gapOpenScore   = matchScore[0][cdx-1] + gapOpen;
    
    if (gapExtendScore < gapOpenScore) {
      gapBScore[0][cdx] = gapExtendScore;
      gapBPath[0][cdx] = GAP_EXTEND;
    }
    else {
      gapBScore[0][cdx] = gapOpenScore;
      gapBPath[0][cdx]  = GAP_OPEN;
    }
    
    gapScore   = matchScore[0][cdx] = matchScore[0][cdx-1] + gap;
    gapClose   = gapBScore[0][cdx];
    
    if (gapScore < gapClose) {
      matchScore[0][cdx] = gapScore;
      matchPath[0][cdx] =  CLOSE_GAP_B;
    }
    else {
      matchScore[0][cdx] = gapClose;
      matchPath[0][cdx]  = GAP_B;
    }

    gapAScore[0][cdx] = INF;
    gapAPath[0][cdx] = CLOSE_GAP_A;
  }

  // Initialize the first column
  for (r = 0; r < rows; r++) {
    rdx = r + 1;
    gapExtendScore = gapAScore[rdx-1][0] + gapExtend;
    gapOpenScore   = matchScore[rdx-1][0] + gapOpen;
    if (gapExtendScore < gapOpenScore) {
      gapAScore[rdx][0] = gapExtendScore;
      gapAPath[rdx][0] = GAP_EXTEND;
    }
    else {
      gapAScore[rdx][0] = gapOpenScore;
      gapAPath[rdx][0]  = GAP_OPEN;
    }

    gapScore   = matchScore[rdx][0] = matchScore[rdx-1][0] + gap;
    gapClose   = gapAScore[rdx][0];
    
    if (gapScore < gapClose) {
      matchScore[rdx][0] = gapScore;
      matchPath[rdx][0] = CLOSE_GAP_A;
    }
    else {
      matchScore[rdx][0] = gapClose;
      matchPath[rdx][0]  = GAP_A;
    }

    gapBScore[rdx][0] = INF;
    gapBPath[rdx][0] = CLOSE_GAP_B;
  }
  

  // Perform alignment
  for (r = 0; r < rows; r++) {
    for (c = 0; c < cols; c++) {
      rdx = r + 1;
      cdx = c + 1;

      // Gaps in sequence A (horizontal gaps)
      gapExtendScore = gapAScore[rdx-1][cdx] + gapExtend;
      gapOpenScore   = matchScore[rdx-1][cdx] + gapOpen;
      if (gapExtendScore < gapOpenScore) {
				gapAScore[rdx][cdx] = gapExtendScore;
				gapAPath[rdx][cdx] = GAP_EXTEND;
      }
      else {
				gapAScore[rdx][cdx] = gapOpenScore;
				gapAPath[rdx][cdx]  = GAP_OPEN;
      }
      
      // Gaps in sequence B (vertical gaps)
      gapExtendScore = gapBScore[rdx][cdx-1] + gapExtend;
      gapOpenScore   = matchScore[rdx][cdx-1] + gapOpen;
      
      if (gapExtendScore < gapOpenScore) {
				gapBScore[rdx][cdx] = gapExtendScore;
				gapBPath[rdx][cdx] = GAP_EXTEND;
      }
      else {
				gapBScore[rdx][cdx] = gapOpenScore;
				gapBPath[rdx][cdx]  = GAP_OPEN;
      }

      // Gaps in match matrix
      gapACloseScore = gapAScore[rdx][cdx];
      gapBCloseScore = gapBScore[rdx][cdx];
      matchScoreVal  = matchScore[rdx-1][cdx-1] + 
				scoreMat[nuc_index[(unsigned char)seqa[r]]][nuc_index[(unsigned char)seqb[c]]];
      gapAScoreVal   = matchScore[rdx-1][cdx] + gap;
      gapBScoreVal   = matchScore[rdx][cdx-1] + gap;

      minScore = std::min(gapACloseScore, std::min(gapBCloseScore, std::min(matchScoreVal, std::min(gapAScoreVal, gapBScoreVal))));

      matchScore[rdx][cdx] = minScore;

      if (minScore == gapACloseScore) 
				matchPath[rdx][cdx] = CLOSE_GAP_A;
      if (minScore == gapBCloseScore)
				matchPath[rdx][cdx] = CLOSE_GAP_B;
      if (minScore == matchScoreVal)
				matchPath[rdx][cdx] = MATCH;
      if (minScore == gapAScoreVal)
				matchPath[rdx][cdx] = GAP_A;
      if (minScore == gapBScoreVal)
				matchPath[rdx][cdx] = GAP_B;

    }
  }
  for (r = 0; r < rows; r++) {
    locations[r] = -1;
  }
  ssize_t score = matchScore[rows][cols];
  GetAlignedLocations(matchPath, gapAPath, gapBPath, rows, cols, locations, 0);
  return score;
}

ssize_t Align(DNASequence &seqa, DNASequence &seqb, 
						ssize_t match, ssize_t mismatch, ssize_t gap, 
						ssize_t *&locations, IntMatrix &scoreMat) {
	if (locations == NULL) {
		locations = new ssize_t[seqa.length];
	}
  if (scoreMat.size() == 0)
    InitScoreMat(scoreMat, match, mismatch);

  ssize_t rows = seqa.length + 1;
  ssize_t cols = seqb.length + 1;
  // Create the score and path arrays
  IntMatrix score;
  IntMatrix path;
  CreateMatrix(score, rows, cols);
  CreateMatrix(path, rows, cols);

  //UNUSED// ssize_t i, j;
  
  ssize_t r, c;
  ssize_t ind1, ind2;
  ssize_t gapA, gapB, mut;
  ssize_t minScore;
  for (r = 1; r < rows; r++) {
    score[r][0] = score[r-1][0] + gap;
    path[r][0] = GAP_A;
  }
  for (c = 1; c < cols; c++) {
    score[0][c] = score[0][c] + gap;
    path[0][c] = GAP_B;
  }
  for (r =1; r < rows; r++) {
    for (c = 1; c < cols; c++) {
      ind1 = (ssize_t) seqa.seq[r-1];
      ind2 = (ssize_t) seqb.seq[c-1];
      ind1 = nuc_index[(unsigned char)ind1];
      ind2 = nuc_index[(unsigned char)ind2];
      mut  = score[r-1][c-1] + scoreMat[ind1][ind2];
      gapA = score[r-1][c] + gap;
      gapB = score[r][c-1] + gap;
      minScore = std::min(mut, std::min(gapA, gapB));
      if (minScore == mut) {
				score[r][c] = mut;
				path[r][c]  = MATCH;
      }
      else if (minScore == gapA) {
				score[r][c] = gapA;
				path[r][c]  = GAP_A;
      }
      else { 
				score[r][c] = gapB;
				path[r][c] = GAP_B;
      }
    }
  }

  for (r = 0; r < seqa.length; r++) {
    locations[r] = -1;
  }

  IntMatrix nogap;
  GetAlignedLocations(path, nogap, nogap, rows, cols, locations, 0, rows-1, cols-1);
  minScore = score[rows-1][cols-1];
  return minScore;

}

ssize_t FitAlign(DNASequence &seqa, DNASequence &seqb, 
							 ssize_t match, ssize_t mismatch, ssize_t gap, 
							 ssize_t *&locations, IntMatrix &scoreMat,
							 IntMatrix &score, IntMatrix &path) {

	// Fit a into b

  if (scoreMat.size() == 0)
    InitScoreMat(scoreMat, match, mismatch);


  ssize_t rows = seqa.length + 1;
  ssize_t cols = seqb.length + 1;
  // Create the score and path arrays

  // perform fitting alignment;
  //UNUSED+// ssize_t j;
  ssize_t i ;
  for (i = 0; i < cols; i++) {
    score[0][i] = 0;
		path[0][i]  = LOCAL_START;
  }
  for (i = 1; i < rows; i++) {
    score[i][0] = 0;
		path[i][0]  = LOCAL_START;
  }
  
  ssize_t r, c;
  ssize_t ind1, ind2;
  ssize_t gapA, gapB, mut;
  ssize_t minScore;
	ssize_t localMinRow, localMinCol;
	ssize_t localMinScore;
	localMinRow = localMinCol = 0;
	localMinScore = 0;
  for (r = 1; r < rows; r++) {
    for (c = 1; c < cols; c++) {
      ind1 = (ssize_t) seqa.seq[r-1];
      ind2 = (ssize_t) seqb.seq[c-1];
      ind1 = nuc_index[(unsigned char)ind1];
      ind2 = nuc_index[(unsigned char)ind2];
      mut  = score[r-1][c-1] + scoreMat[ind1][ind2];
      gapA = score[r-1][c] + gap;
      gapB = score[r][c-1] + gap;
      minScore = std::min(mut, std::min(gapA, gapB));
			if (minScore >= 0) {
				score[r][c] = 0;
				path[r][c]  = LOCAL_START; // jump to start
			}
      else if (minScore == mut) {
				score[r][c] = mut;
				path[r][c]  = MATCH;
      }
      else if (minScore == gapA) {
				score[r][c] = gapA;
				path[r][c]  = GAP_A;
      }
      else { 
				score[r][c] = gapB;
				path[r][c] = GAP_B;
      }
			if (localMinScore > minScore) {
				localMinScore = minScore;
				localMinCol   = c;
				localMinRow   = r;
			}
    }
  }
  minScore = INF;
	ssize_t minCol = -1;
	ssize_t minRow = rows-1;
	for (c = 0; c < cols; c++) {
    if (minScore > score[rows-1][c]) {
			minCol   = c;
			minScore = score[rows-1][c];
		}
  }

	for (r = 0; r < rows; r++ ){
		if (minScore > score[r][cols-1]) {
			minRow = r;
			minCol = cols-1;
			minScore = score[r][cols-1];
		}
	}
	if (minScore > localMinScore) {
		minRow = localMinRow;
		minCol = localMinCol;
		minScore = localMinScore;
	}
  //for (r = 0; r < rows; r++) {
	for (r = 0; r < rows-1; r++) {
    locations[r] = -1;
  }

  IntMatrix nopath;
  GetAlignedLocations(path, nopath, nopath, rows, cols, locations, DO_LOCAL_ALIGN, minRow, minCol);

  minScore = score[minRow][minCol];
  
  return minScore;
}

ssize_t CoreAlign(DNASequence &seqa, DNASequence &seqb,
								ssize_t match, ssize_t mismatch, ssize_t gap, 
								ssize_t *&locations, IntMatrix &score, IntMatrix &path,
								ssize_t alternative, // 0 for smith-waterman, +inf for needleman wunsch
								ssize_t &globalMinRow, ssize_t &globalMinCol) {
  IntMatrix scoreMat;
  InitScoreMat(scoreMat, match, mismatch);
  if (locations == NULL) {
    locations = new ssize_t[seqa.length];
  }
  ssize_t rows = seqa.length + 1;
  ssize_t cols = seqb.length + 1;

  // The path and match arrays should have been initialized 
  // outside of this function.
  ssize_t r, c;
  ssize_t ind1, ind2;
  ssize_t gapA, gapB, mut;
  ssize_t minScore = 0;
  ssize_t globalMinScore;
  globalMinScore= 0; // TODO: check if this should be +inf
  globalMinRow  = 0;
  globalMinCol  = 0;
  for (r =1; r < rows; r++) {
    for (c = 1; c < cols; c++) {
      ind1 = (ssize_t) seqa.seq[r-1];
      ind2 = (ssize_t) seqb.seq[c-1];
      ind1 = nuc_index[(unsigned char)ind1];
      ind2 = nuc_index[(unsigned char)ind2];
      mut  = score[r-1][c-1] + scoreMat[ind1][ind2];
      gapA = score[r-1][c] + gap;
      gapB = score[r][c-1] + gap;
      minScore = std::min( alternative, std::min(mut, std::min(gapA, gapB)));
      if (minScore == mut) {
				score[r][c] = mut;
				path[r][c]  = MATCH;
      }
      else if (minScore == gapA) {
				score[r][c] = gapA;
				path[r][c]  = GAP_A;
      }
      else if (minScore == gapB) { 
				score[r][c] = gapB;
				path[r][c] = GAP_B;
      } 
      else {
				score[r][c] = 0;
				path[r][c]  = LOCAL_START;
      }
      if (minScore < globalMinScore) {
				globalMinScore = minScore;
				globalMinCol = c;
				globalMinRow = r;
      }
    }
  }
  for (r = 0; r < seqa.length; r++) {
    locations[r] = -1;
  }
  return minScore;
}

ssize_t AffineLocalAlign(DNASequence &seqa, DNASequence &seqb, 
											 ssize_t match, ssize_t mismatch, ssize_t gap, 
											 ssize_t gapOpen, ssize_t gapExtend,
											 ssize_t *&locations, IntMatrix & scoreMat) {
  if (scoreMat.size() == 0)
    init(match, mismatch,scoreMat);

  ssize_t rows = seqa.length;
  ssize_t cols = seqb.length;

  IntMatrix gapAScore, gapBScore, matchScore;
  IntMatrix gapAPath, gapBPath, matchPath;
  ssize_t matchScoreVal, gapAScoreVal, gapBScoreVal, 
    gapACloseScore, gapBCloseScore, gapOpenScore, gapExtendScore;
	//UNUSED// int gapScore;


  CreateMatrix(gapAScore, rows+1, cols+ 1);
  CreateMatrix(gapBScore, rows+1, cols+ 1);
  CreateMatrix(matchScore, rows+1, cols+ 1);
  
  CreateMatrix(gapAPath, rows+1, cols+ 1);
  CreateMatrix(gapBPath, rows+1, cols+ 1);
  CreateMatrix(matchPath, rows+1, cols+ 1);


  ssize_t r, c;
  ssize_t rdx, cdx;
  // Initialize the first row/column of all matricies.
  gapAScore[0][0]  = INF;
  gapAPath[0][0]   = 0;
  gapBScore[0][0]  = INF;
  gapBPath[0][0]   = 0;

  
  matchScore[0][0] = 0;
  //UNUSED+// ssize_t  gapClose;
  ssize_t minScore;
  
  // Initialize the first row  
  for (c = 0; c < cols; c++) {
    cdx = c + 1;
    gapBScore[0][cdx] = 0;
    gapBPath[0][cdx]  = GAP_EXTEND;
    
    matchScore[0][cdx] = 0;
    matchPath[0][cdx] =  CLOSE_GAP_B;

    gapAScore[0][cdx] = INF;
    gapAPath[0][cdx] = CLOSE_GAP_A;
  }

  // Initialize the first column
  for (r = 0; r < rows; r++) {
    rdx = r + 1;

    gapAScore[rdx][0] = 0;
    gapAPath[rdx][0] = GAP_EXTEND;

    matchScore[rdx][0] = 0;
    matchPath[rdx][0] = CLOSE_GAP_A;

    gapBScore[rdx][0] = INF;
    gapBPath[rdx][0] = CLOSE_GAP_B;
  }
  

  // Perform alignment
  for (r = 0; r < rows; r++) {
    for (c = 0; c < cols; c++) {
      rdx = r + 1;
      cdx = c + 1;

      // Gaps in sequence A (horizontal gaps)
      gapExtendScore = gapAScore[rdx-1][cdx] + gapExtend;
      gapOpenScore   = matchScore[rdx-1][cdx] + gapOpen;
      if (gapExtendScore < gapOpenScore) {
				gapAScore[rdx][cdx] = gapExtendScore;
				gapAPath[rdx][cdx] = GAP_EXTEND;
      }
      else {
				gapAScore[rdx][cdx] = gapOpenScore;
				gapAPath[rdx][cdx]  = GAP_OPEN;
      }
      
      // Gaps in sequence B (vertical gaps)
      gapExtendScore = gapBScore[rdx][cdx-1] + gapExtend;
      gapOpenScore   = matchScore[rdx][cdx-1] + gapOpen;
      
      if (gapExtendScore < gapOpenScore) {
				gapBScore[rdx][cdx] = gapExtendScore;
				gapBPath[rdx][cdx] = GAP_EXTEND;
      }
      else {
				gapBScore[rdx][cdx] = gapOpenScore;
				gapBPath[rdx][cdx]  = GAP_OPEN;
      }

      // Gaps in match matrix
      gapACloseScore = gapAScore[rdx][cdx];
      gapBCloseScore = gapBScore[rdx][cdx];
      matchScoreVal  = matchScore[rdx-1][cdx-1] + 
				scoreMat[nuc_index[(unsigned char)seqa[r]]][nuc_index[(unsigned char)seqb[c]]];
      gapAScoreVal   = matchScore[rdx-1][cdx] + gap;
      gapBScoreVal   = matchScore[rdx][cdx-1] + gap;

      minScore = std::min((ssize_t)0.0,
													std::min(gapACloseScore, 
																	 std::min(gapBCloseScore, matchScoreVal)));


      /*            std::cout << rdx << " " << cdx << " " << seqa[r] << " " << seqb[c] << " "
										<< scoreMat[nuc_index[(unsigned char)seqa[r]]][nuc_index[(unsigned char)seqb[c]]] 
										<< " " << matchScoreVal 
										<< " " << gapACloseScore 
										<< " " << gapBCloseScore 
										<< " gap: " << gapAScore[rdx][cdx] 
										<< " " << gapBScore[rdx][cdx] << std::endl;
      */
      matchScore[rdx][cdx] = minScore;
      if (minScore == 0) 
				matchPath[rdx][cdx] = LOCAL_START;
      else if (minScore == gapACloseScore) 
				matchPath[rdx][cdx] = CLOSE_GAP_A;
      else if (minScore == gapBCloseScore)
				matchPath[rdx][cdx] = CLOSE_GAP_B;
      else if (minScore == matchScoreVal)
				matchPath[rdx][cdx] = MATCH;
    }
    /*      std::cout << "match " << std::endl;
						PrintMatrix(matchScore, std::cout, 2);
						std::cout << "gapa: " << std::endl;
						PrintMatrix(gapAScore, std::cout, 2);
						std::cout << "gapb: " << std::endl;
						PrintMatrix(gapBScore, std::cout, 2 );
    */
  }

  // Find minimum score.
  ssize_t minRow, minCol;
  minScore = 1;
  minRow = 0; minCol = 0;
  for (r = 1; r <= rows; r++ )
    for (c = 1; c <= cols; c++) 
      if (minScore > matchScore[r][c]) {
				minRow = r; minCol = c; minScore = matchScore[r][c];
      }
      
  std::cout << "minscore: " << minScore << " " << minRow << " " << minCol << " of " << rows << " " << cols << std::endl;
  for (r = 0; r < rows; r++) {
    locations[r] = -1;
  }
  ssize_t score = matchScore[minRow][minCol];
	std::cout << "match " << std::endl;
	PrintMatrix(matchScore, std::cout, 2);
	PrintMatrix(matchPath, std::cout, 2);
	std::cout << "gapa: " << std::endl;
	PrintMatrix(gapAPath, std::cout, 2);
	std::cout << "gapb: " << std::endl;
	PrintMatrix(gapBPath, std::cout, 2 );
  GetAlignedLocations(matchPath, gapAPath, gapBPath, rows, cols, locations, 1, minRow, minCol);
  ssize_t loc;
  for (loc = 0; loc < rows-1; loc++ )
    std::cout << locations[loc] << " ";
  std::cout << std::endl;
  return score;
}

ssize_t LocalAlign(DNASequence &seqa, DNASequence &seqb, 
								 ssize_t match, ssize_t mismatch, ssize_t gap, 
								 ssize_t *&locations, IntMatrix &scoreMat) {

  if (scoreMat.size() == 0)
    InitScoreMat(scoreMat, match, mismatch);
  ssize_t rows = seqa.length + 1;
  ssize_t cols = seqb.length + 1;
  // Create the score and path arrays
  IntMatrix score;
  IntMatrix path;
  CreateMatrix(score, rows, cols);
  CreateMatrix(path, rows, cols);
  //UNUSED+// ssize_t j;
  ssize_t i ;
  ssize_t minScore;
  // Initialize gaps.
  for (i = 0; i < rows; i++) {
    score[i][0] = 0;
    path[i][0]  = GAP_A;
  }

  for (i = 0; i < cols; i++) {
    score[0][i] = 0;
    path[0][i] = GAP_B;
  }
  ssize_t globalMinRow, globalMinCol;
  CoreAlign(seqa, seqb, 
						match, mismatch, gap, 
						locations, score, path, 0, globalMinRow, globalMinCol);
  IntMatrix noPath;
  GetAlignedLocations(path, noPath, noPath, 
											rows, cols, 
											locations, 1, globalMinRow, globalMinCol);

  minScore = score[globalMinRow][globalMinCol];
  return minScore;
}

ssize_t KBandRow(ssize_t r, ssize_t c, ssize_t k) {
  return r;
}

ssize_t KBandCol(ssize_t r, ssize_t c, ssize_t k) {
  return c + k - r + 1;
}

ssize_t OverlapAlign(DNASequence &seqa, DNASequence &seqb,
									 ssize_t match, ssize_t mismatch, ssize_t gap, 
									 ssize_t *&locations, IntMatrix &scoreMat) {
  
  if (scoreMat.size() == 0)
    InitScoreMat(scoreMat,match, mismatch);

  ssize_t rows = seqa.length + 1;
  ssize_t cols = seqb.length + 1;
  // Create the score and path arrays
  IntMatrix score;
  IntMatrix path;
  CreateMatrix(score, rows, cols);
  CreateMatrix(path, rows, cols);
  //UNUSED+// ssize_t j;
  ssize_t i ;
  ssize_t minScore;
  // Initialize gaps.
  for (i = 0; i < rows; i++) {
    score[i][0] = 0;
    path[i][0]  = GAP_A;
  }

  for (i = 0; i < cols; i++) {
    score[0][i] = gap * i;
    path[0][i] = GAP_B;
  }
  ssize_t globalMinRow, globalMinCol;
  CoreAlign(seqa, seqb, 
						match, mismatch, gap, 
						locations, score, path, INF, globalMinRow, globalMinCol);

  // To get the overlap alignment, find the minimum score along the 
  // rightmost column.
  ssize_t globalMinScore;
  globalMinScore = INF;
  globalMinRow   = -1;
	globalMinCol   = -1;
  for (i = 0; i < cols; i++) {
    if (score[rows-1][i] < globalMinScore) {
      globalMinScore = score[rows-1][i];
      globalMinCol   = i;
    }
  }
  
  IntMatrix nopath;
  GetAlignedLocations(path, nopath, nopath,
											rows, cols, 
											locations, 0, rows-1, globalMinCol);

  minScore = score[rows-1][globalMinCol];
  return minScore;
}

ssize_t BandedAlign(DNASequence &seqa, DNASequence &seqb, 
									ssize_t match, ssize_t mismatch, ssize_t gap, 
									ssize_t k,
									ssize_t *&locations, 
									IntMatrix &score, 
									IntMatrix &path,
									IntMatrix &scoreMat,
									ssize_t *optScores
									) {
	if ( szabs(seqa.length - seqb.length) > k) {
		return   szabs(seqa.length - seqb.length) * gap;
	}

  if (scoreMat.size() == 0)
    InitScoreMat(scoreMat,match, mismatch);

  ssize_t rows, cols, bandedCols;
  ssize_t r, c;

	// Determine the dimensionality of the alignment.

	if (seqa.length < seqb.length) {
		cols = std::min(seqa.length + 1 + k, seqb.length + 1);
		rows = seqa.length + 1;
	}
	else {
		rows = std::min(seqb.length + 1 + k, seqa.length + 1);
		cols = seqb.length;
	}

  bandedCols = 2*k+2+1;
	
	
	// Resize the matrices if necessary.  These are stored
	// as static variables so that they do not need to be
	// reallocated in between function calls.
	if (score.size() < rows or 
			(score.size() > 0 and  score[0].size() < bandedCols)) {
		CreateMatrix(score, rows, bandedCols);
		CreateMatrix(path, rows, bandedCols);
	}

  // Initialize so I can see what is being assigned
  for (r = 0; r < rows; r++) 
    for (c = 0; c < bandedCols; c++)
      score[r][c] = 0; 

  ssize_t cb;
  // Initialize columns with gap score
  for (c = 1; c < k+1; c++) {
    cb = KBandCol(0, c, k);
    score[0][cb] = c * gap;
    path[0][cb]  = BAND_GAP_B;
  }

  for (r = 1; r < rows; r++) {
    c = KBandCol(r, -1, k);
    if (c >= 0) {
      score[r][c] = r*gap;
      path[r][c] = BAND_GAP_A;
    }
		c = KBandCol(r, r+k+1, k);
		if (c < 2*k+3){
			score[r][c] = r*gap;
			path[r][c] = BAND_GAP_A;
		}
  }
  // Initialize the diagonal boundaries of the score matrix.
  for (r = k+1; r < rows; r++) {
    c = KBandCol(r, r-k-1, k);
    if (c >= 0)
      score[r][c] = INF;
  }
  for (r = 0; r < rows - k+1; r++) {
    c = KBandCol(r, r+k+1, k);
    if (c < bandedCols)
      score[r][c] = INF ;
  }
  // Perform banded alignment
  //UNUSED// ssize_t bc; // banded column
  ssize_t matchScore, gapAScore, gapBScore, minScore;
  ssize_t rdx, cdx, bcdx;
  for (r = 0; r < rows-1; r++) { // iterate over all chars in a
    for (c = std::max(r - k, (_SSZT_) 0);  // TODO: remove cast
				 c < std::min((_SSZT_) (r+k+1), seqb.length);  // TODO: remove cast
				 c++) {
      rdx = r + 1;
      cdx = KBandCol(r+1,c + 1,k);
      if (c < 0) continue;
      if (c >= seqb.length) continue;
      bcdx = KBandCol(rdx-1, c+1-1, k);
      matchScore = score[rdx-1][bcdx] + scoreMat[nuc_index[(unsigned char)seqa.seq[r]]][nuc_index[(unsigned char)seqb.seq[c]]];
      bcdx = KBandCol(rdx-1, c+1, k);
      gapAScore  = score[rdx-1][bcdx] + gap;
      bcdx = KBandCol(rdx, c+1-1, k);
      gapBScore  = score[rdx][bcdx] + gap;
      
      minScore = std::min(matchScore, std::min(gapAScore, gapBScore));
      score[rdx][cdx] = minScore;
      if (minScore == matchScore) 
				path[rdx][cdx] = BAND_MATCH;
      else if (minScore == gapAScore)
				path[rdx][cdx] = BAND_GAP_A;
      else 
				path[rdx][cdx] = BAND_GAP_B;
    }
  }
	if (locations != NULL) {
		GetBandedAlignedLocations(score,
															path, 
															rows, bandedCols, 
															locations,
															optScores, // Don't store scores
															k,
															rows-1, 
															KBandCol(rows-1, cols-1,k));
	}
  return score[rows-1][KBandCol(rows-1, cols-1,k)];
}



ssize_t LocalAlign(DNASequence &seqa, DNASequence &seqb, 
								 Score &scoreMat, ssize_t *&locations) {
  return LocalAlign(seqa, seqb, 0, scoreMat.gapOpen, 0, locations, scoreMat.scoreMat);
}

ssize_t AffineLocalAlign(DNASequence &seqa, DNASequence &seqb,
											 Score &scoreMat, ssize_t *&locations) {
  return AffineLocalAlign(seqa, seqb, 
													1000, 1000, 1000,
													scoreMat.gapOpen, scoreMat.gapExtend,
													locations, scoreMat.scoreMat);
}


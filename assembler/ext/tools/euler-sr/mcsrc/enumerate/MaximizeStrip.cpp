/***************************************************************************
 * Title:          MaximizeStrip.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "MaximizeStrip.h"
#include "alignutils.h"
#include "TupleLib.h"

#include <stdio.h>
#include <iostream>
#include <ios>
#include <iomanip>
#include <cmath>

ssize_t Compatible(ssize_t *perm, 
	       ssize_t *refStartPos, ssize_t *refEndPos,
	       ssize_t *qryStartPos, ssize_t *qryEndPos,
	       ssize_t index1, ssize_t index2) {
  ssize_t diff, val;
  diff= perm[index2] - perm[index1];
  // Cases for compatibilty: 
  //   1.  Two markers are in the same direction, and they follow eachother in the query strand
  //       (this means that a sequence is colinear)
  //   2.  The previous marker is in the forward direction, and the current marker is in
  //       the reverse direction, and the position in the query sequence of the current 
  //       marker is beyond that of the previous marker.
  //   3.  The pervious marker is in the reverse direction, and the current marker is in
  //       the forward direction, and the current marker is in front of the previous in the query 
  //       strand.

  if (Sign(perm[index1]) == 1 && Sign(perm[index2]) == 1 && 
      qryStartPos[index2] > qryEndPos[index1])
    return 1;
  
  if (Sign(perm[index1]) == 0 && Sign(perm[index2]) == 0 && 
      qryStartPos[index2] < qryEndPos[index1])
    return 1;
  
  if (Sign(perm[index1]) != Sign(perm[index2]) && 
      qryStartPos[index2] > qryEndPos[index1])
    return 1;

  return 0;
}

void GetPermutation(ssize_t *perm, ssize_t *permSizes, ssize_t length, ssize_t startPos,
		    ssize_t*maxPerm, ssize_t &maxPermLength, ssize_t &totalPermLength,
		    double **scores, ssize_t *locations) {

  // The locations array contains a list of back pointers.  Follow the 
  // back pointers (really just indices)
  // to the beginning of the list, storing the locations
  // of the pointers in an array.
  ssize_t row;
  ssize_t maxRow, index;

  // Count the number of elements in the permutation
  maxPermLength = 0;
  totalPermLength = 0;
  row = startPos; // The optimally ending row is stored in the last spot
  while (locations[row] >= 0) {
    assert(row >= 0 && row < length || printf("row: %d length: %d\n", row, length) == 0);
    maxPermLength++;
    totalPermLength +=permSizes[row];
    row = locations[row];
  }
  totalPermLength += permSizes[row];

  // Now fill the permutations array
  row = startPos;
  index = maxPermLength;
  while (locations[row] >= 0) {
    index--;
    maxPerm[index] = row;
    row = locations[row];
  }
  assert(index == 0 || printf("maxpermlength: %d\n", index) == 0);
}


double GetDist(ssize_t x1, ssize_t y1, ssize_t x2, ssize_t y2) {
	double dx = x2-x1;
	double dy = y2-y1;
  return std::sqrt( dx*dx + dy*dy );
}

double ScoreGap(ssize_t x1, ssize_t y1, ssize_t x2, ssize_t y2, double g, double i) {
  ssize_t dist, ins;
  dist = std::min(std::abs(y2 - y1), std::abs(x2 - x1));
  ins  = std::max(std::abs(y2 - y1), std::abs(x2 - x1)) - dist;
  return dist*g + ins*i;
}

double GetMinThreading(ssize_t *perm, ssize_t *permSize, 
		      ssize_t *refStartPos, ssize_t *refEndPos, 
		      ssize_t *qryStartPos, ssize_t *qryEndPos,
		      ssize_t length, 
		      double s, double g, double d,
		      ssize_t *maxPerm, ssize_t &maxPermLength, ssize_t &totalPermLength) {
  // Find minimum cost path from first marker to last.
  ssize_t *locations;
  double **score = CreateMatrix<double>(length, length);
  locations = new ssize_t[length];
  ssize_t row, col;
  ssize_t idx;
  score[0][0] = 0;
  locations[0]= -1;
  
  for (row = 1; row < length ; row++) {
    // row is index of last permutation
    // col is index of marker to include in LIS

    // score[row][col] holds the optimal score to get too the marker at pos 'row' from
    // the marker at pos [col]. 
    // Initialize it to only be that row.
    score[row][row] = permSize[row]*s + ScoreGap(refEndPos[0], qryEndPos[0], 
						 refStartPos[row], qryStartPos[row], 
						 g, d);
    locations[row]  = 0;
    for (col = 0; col < row ; col++) {
      // Look to see if the value at the current position can be 
      // added to the set of current values.
      if (Compatible(perm, refStartPos, refEndPos,
		 qryStartPos, qryEndPos,
		     col, row )) {
	score[row][col] = score[col][col] + permSize[col]*s+ 
	  ScoreGap(refEndPos[col], qryEndPos[col], 
		   refStartPos[row], qryStartPos[row],
		   g, d);
	if (score[row][row] < score[row][col]) {
	  score[row][row] = score[row][col];
	  locations[row] = col;
	}
      }
    }
  }
  GetPermutation(perm, permSize, length, length-1,
		 maxPerm, maxPermLength, totalPermLength,
		 score,  locations);

  std::cout << "scores of min threading: ";
  ssize_t pos = length-1;
  while (pos >= 0) {
    std::cout << "(" << pos << ")" << score[pos][pos] << "  ";
    pos = locations[pos];
  }
  double minDist = score[length-1][length-1];
  DeleteMatrix(score, length);
  delete []locations;
  return minDist;
}


double GetLongestIncreasingSubset(ssize_t *perm, ssize_t *permSize, 
				 ssize_t *refStartPos, ssize_t *refEndPos, 
				 ssize_t *qryStartPos, ssize_t *qryEndPos,
				 double s, double g, double d,
				 ssize_t length, ssize_t diffThreshold,
				 ssize_t *maxPerm, ssize_t &maxPermLength, ssize_t &totalPermLength) {

  ssize_t *locations;
  double **score = CreateMatrix<double>(length, length);
  locations = new ssize_t[length];
  ssize_t row, col;
  ssize_t sufficient = 0;
  score[0][0] = 0;
  locations[0]= -1;
  double maxScore = 0;
  ssize_t   maxRow   = 0;
  for (row = 1; row < length && ! sufficient ; row++) {
    // row is index of last permutation
    // col is index of marker to include in LIS

    // score[row][row] holds the maximum score for permutation at row
    // Initialize it to only be that row.
    score[row][row] = ScoreGap(refEndPos[0], refEndPos[0], qryStartPos[row], qryStartPos[row], 0, 0);
    locations[row] = -1;
    for (col = 0; col < row && ! sufficient ; col++) {
      // Look to see if the value at the current position can be 
      // added to the set of current values.
      score[row][col] = 
	score[col][col] +  // previous best score for this marker
	permSize[col]*s +  // score of this marker
	ScoreGap(refEndPos[col], qryEndPos[col],  // cost to reach this marker
		 refStartPos[row], qryStartPos[row],
		 g, d);
      if (score[row][row] < score[row][col]) {
	score[row][row] = score[row][col];
	locations[row] = col;
      }
      if (score[row][row] > maxScore) {
	maxScore = score[row][row];
	maxRow   = row;
      }
    }
  }
  GetPermutation(perm, permSize, length, maxRow,
		 maxPerm, maxPermLength, totalPermLength,
		 score,  locations);
  
  maxScore = score[maxRow][maxRow];
  DeleteMatrix(score, length);
  delete []locations;
  return maxScore;
}
    
ssize_t RemoveStrips(T_Strips::iterator &stripIt, T_Strips *&strips, 
		 ssize_t *permutation, ssize_t length, 
		 ssize_t *maxPermutations, ssize_t maxPermLength) {
  ssize_t lisPos = 0; 
  ssize_t end = 0;
  ssize_t numerased = 0;
  ssize_t winPos;
  ssize_t removedStrips = 0;
  for (winPos = 0; winPos < length && stripIt != strips->end() ; winPos++) {
    assert(stripIt != strips->end());
    assert(lisPos <= maxPermLength || printf("lispos: %d maxperml: %d\n", lisPos, maxPermLength)==0);
    if (lisPos < maxPermLength && 
	stripIt->startQryEnum  == permutation[maxPermutations[lisPos]]) {
      // Found an enumertation that is part of the LIS, leave it
      // in the enueration
      ++stripIt;
      lisPos++;
    } 
    else { 
      // The current permutation is not part of the lis, remove it.
      removedStrips++;
      stripIt = strips->erase(stripIt);
    }
  }
  return removedStrips;
}


void MinThreadingPath(T_Strips *&strips, ssize_t window) {
  if (window < 2) 
    return;
  if (strips->size() < 2)
    return;

  // Find the path through markers (k-mers) that travels
  // the least distance from the first to last.

  ssize_t *permutation        = new ssize_t[window];
  ssize_t *refStartLocations  = new ssize_t[window];
  ssize_t *refEndLocations    = new ssize_t[window];
  ssize_t *qryStartLocations  = new ssize_t[window];
  ssize_t *qryEndLocations    = new ssize_t[window];
  ssize_t *permSize           = new ssize_t[window];
  ssize_t *maxPermutations    = new ssize_t[window];
  ssize_t maxPermLength = 0;
  ssize_t totalPermLength = 0;
  T_Strips::iterator stripIt, it, permStartIt, permIt;
  ssize_t stripPos, i;
  it = strips->begin();
  stripIt = strips->begin();
  ssize_t lisPos;
  ssize_t removedStrips;
  //  std::cout << "beginning ms" << std::endl;
  ssize_t done = 0;
  std::cout << "MS: " << strips->size() << " start " << std::endl;
  while (! done && strips->size() > 1 && stripIt != strips->end()) {
    assert(stripIt != strips->end() 
	   || printf("strippos: %d size: %d  window: %d\n", 
		     stripPos, strips->size(), window) == 0);

    // Copy the permutation into temporary arrays that are
    // passed to the LIS function.
    it = stripIt;
    ssize_t minRef, minQry, maxRef, maxQry;
    ssize_t minEnum, maxEnum;

    for (i = 0; i < window && it != strips->end(); i++) {
      refStartLocations[i] = it->startRefPos;
      refEndLocations[i]   = it->startRefPos + it->size;
      qryStartLocations[i] = it->startQryPos;
      qryEndLocations[i]   = it->startQryPos + it->size;
      permutation[i]       = it->startQryEnum;
      qryStartLocations[i] = it->startQryPos;
      permSize[i]          = it->size;
      ++it;
    }

    if (it == strips->end()) {
      // This is the last iteration through the loop
      window = i; // just in case the window extends past the end of the list
      //   done = 1;
    }

    // Find the indices of the longset increasng subset within the
    // array permutation.
    double score;
    score = GetMinThreading(permutation, permSize, 
			    refStartLocations, refEndLocations,
			    qryStartLocations, qryEndLocations,
			    window, 
			    1, 0.0, -0.1, 
			    maxPermutations, 
			    maxPermLength, 
			    totalPermLength);

    // Advance the iterator to the first position in strips 
    // that corresponds to the first permutation in the LIS.
    // stripIt is at the first element in a longest increaseng subset.
    // Merge all elements of the longest increasng subset into stripIt.

    removedStrips = RemoveStrips(stripIt, strips,
				  permutation, window,
				  maxPermutations, maxPermLength);
  }
  delete [] permutation;
  delete [] permSize;
  delete [] maxPermutations;
  delete [] refStartLocations;
  delete [] refEndLocations;
  delete [] qryStartLocations;
  delete [] qryEndLocations;
}



void MaximizeStrips(T_Strips *strips, ssize_t window, ssize_t threshold) {
  if (window < 2) 
    return;
  if (strips->size() < 2)
    return;
  ssize_t *permutation    = new ssize_t[window];
  ssize_t *refStartLocations  = new ssize_t[window];
  ssize_t *refEndLocations    = new ssize_t[window];
  ssize_t *qryStartLocations  = new ssize_t[window];
  ssize_t *qryEndLocations    = new ssize_t[window];
  ssize_t *permSize           = new ssize_t[window];
  ssize_t *maxPermutations    = new ssize_t[window];
  ssize_t maxPermLength = 0;
  ssize_t totalPermLength = 0;
  T_Strips::iterator stripIt, it, permStartIt, permIt;
  ssize_t stripPos, i;
  it = strips->begin();
  stripIt = strips->begin();
  ssize_t lisPos;
  ssize_t removedStrips;
  //  std::cout << "beginning ms" << std::endl;
  ssize_t done = 0;
  std::cout << "MS: " << strips->size() << " start " << std::endl;
  while (! done && strips->size() > 1 && stripIt != strips->end()) {
    assert(stripIt != strips->end() 
	   || printf("strippos: %d size: %d  window: %d\n", 
		     stripPos, strips->size(), window) == 0);
    
    // Copy the permutation into temporary arrays that are
    // passed to the LIS function.
    it = stripIt;
    ssize_t minRef, minQry, maxRef, maxQry;
    ssize_t minEnum, maxEnum;
    
    for (i = 0; i < window && it != strips->end(); i++) {
      refStartLocations[i] = it->startRefPos;
      refEndLocations[i]   = it->endRefPos;
      qryStartLocations[i] = it->startQryPos;
      qryEndLocations[i]   = it->endQryPos;
      permutation[i]       = it->startQryEnum;
      qryStartLocations[i] = it->startQryPos;
      permSize[i]          = it->size;
      ++it;
    }

    if (it == strips->end()) {
      // This is the last iteration through the loop
      window = i; // just in case the window extends past the end of the list
      //   done = 1;
    }
  	
    // Find the indices of the longset increasng subset within the
    // array permutation.
    double score;
    score = GetLongestIncreasingSubset(permutation, permSize, 
				       refStartLocations, refEndLocations,
				       qryStartLocations, qryEndLocations,
				       1, -0.05, -0.05, // take a guess for now
				       window, threshold,
				       maxPermutations, 
				       maxPermLength, 
				       totalPermLength);

    // Advance the iterator to the first position in strips 
    // that corresponds to the first permutation in the LIS.
    // stripIt is at the first element in a longest increaseng subset.
    // Merge all elements of the longest increasng subset into stripIt.
    
    removedStrips = RemoveStrips(stripIt, strips,
				 permutation, window,
				 maxPermutations, maxPermLength);
  }
  delete [] permutation;
  delete [] permSize;
  delete [] maxPermutations;
  delete [] refStartLocations;
  delete [] refEndLocations;
  delete [] qryStartLocations;
  delete [] qryEndLocations;
}

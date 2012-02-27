/***************************************************************************
 * Title:          EmbossAlignment.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "EmbossAlignment.h"
#include "DNASequence.h"
#include "align/alignutils.h"


void EmbossAlignment::
CalculateAverageScore(ssize_t window, double *&scores, ssize_t & scoresLength,
		      FloatMatrix &scoreMat, double gapOpen, double gapExtend,
		      DNASequence &refSeq, DNASequence &qrySeq,
		      ssize_t *markers, ssize_t numMarkers) {
  ssize_t i, l, p, q, pl, s;
	s = 0;

  // l: index into the location array.  
  // r: position in the reference sequence
  double score;
  assert(locations != NULL);
  ssize_t qs, qe;
  ssize_t inGap;
  unsigned char refc, qryc;
  ssize_t localWindow;

  ssize_t markerPos = 0;
  // determine the alignment length;
  ssize_t alignLength = 1;
  for (i = 1;i < length; i++) {
    if (locations[i] == -1)
      alignLength++;
    else {
      //      std::cout << "adidng: " << locations[i] - locations[i-1] << std::endl;
      if (locations[i-1] != -1)
	alignLength+= locations[i] - locations[i-1];
      else 
	alignLength += 1;
    }
  }

  scoresLength = alignLength - window + 1;
  scores = new double[scoresLength+1];
  
  // Fill out the alignment arrays, either the character aligned, or -1 for gap.

  qs = locations[0];  
  for(i = length-1; i >= 0 && locations[i] == -1; i--) ;
  qe = locations[i];

  if (qe - qs - window <= 0) {
    scores = NULL;
    return;
  }

  // Score all aligned positions in  the reference sequence.
  l = 0;
  for (q = qs; q <= qe - window; ) {
    // score a window in the alignment
    score = 0;
    pl = l;
    inGap = 0;
    localWindow = window;

    for (p = 0; p < localWindow;) {
      // Is the current position in a gap? 
      if (locations[pl] == -1 or
	  p + q < locations[pl]) {
	// Score the gap
	/*
	std::cout << "gap at location " << pl << " (" << locations[pl] << ")" << std::endl;
	*/
	if (inGap) // are we already in a gap?
	  score += gapExtend;
	else {
	  // Gap is opening.
	  score += gapOpen;
	  inGap = 1;
	}
      }
      else if (p + q == locations[pl])  {
	refc = (unsigned char) refSeq.seq[refStart + pl]; //
	qryc = (unsigned char) qrySeq.seq[qryStart + locations[pl]];
	/*	
	  std::cout << "scoring : " << refc << "(" <<  refStart + pl << ") ";
	  std::cout << " : " << qryc << "(" << qryStart + locations[pl] << ") " << std::endl;
	*/
	score += scoreMat[nuc_index[refc]][nuc_index[qryc]];
	inGap = 0;
	// 
      }
      
      // Advance in the locations
      if (locations[pl] == -1) {
	pl++;
	localWindow--;
      }
      else if (locations[pl] == p + q) {
	pl++;
	p++;
      }
      else {
	p++;
      }
    }
    /*
      std::cout << " score for pos: " << s << "/" << scoresLength << " : " <<  score<< std::endl;
    */
    scores[s] = score / window;
    s++;
    ssize_t lt = l;

    if (locations[lt] == -1 or
	locations[lt] == q ) {
      l++;
      if (markerPos < numMarkers && 
	  markers[markerPos] == l + refStart) {
	std::cout << "marker: " << markers[markerPos] << " -> " << l << std::endl;
	markers[markerPos] = l;
	markerPos++;
      }
    }
    if (locations[lt] != -1) 
      q++;
  }
}

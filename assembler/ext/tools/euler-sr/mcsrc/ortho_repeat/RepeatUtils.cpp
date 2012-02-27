/***************************************************************************
 * Title:          RepeatUtils.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "RepeatUtils.h"

#incluede "TupleLib.h"
#include "alignutils.h"


ssize_t DetermineOrientation(DNASequence &ref, DNASequence &qry, ssize_t *&locations, double **scoreMat) {
  ssize_t *forwardLocations, reverseLocations;
  DNASequence qryRC;
  double forwardScore, reverseScore;
  MakeRC(qry, qryRC);
  forwardScore = AffineAlign(ref, qry, 0, 0, 0, 300, 10, forwardLocations, scoreMat);
  reverseScore = AffineAlign(ref, qryRC, 0, 0, 0, 300, 10, reverseLocations, scoreMat);
  
  if (forwardScore > reverseScore) {
    delete[] forwardLocations;
    locations = reverseLocations;
    return -1;
  }
  else {
    delete[] reverseLocations;
    locations = forwardLocations;
    return 1;
  }
}


void  StoreProfile(DNASequence &consensus, DNASequence &repeatSeq, 
		   ProfileCount &profileCount, double **scoreMatPtr) {

  ssize_t *conToRepeat;
  ssize_t orientation;

  orientation = DetermineOrientation(consensus, repeatSeq, conToRepeat, scoreMat);
  ssize_t pos;
  // Get the starting and ending positions (pos of first non -1 and last non -1)
  ssize_t start, end;
  start = 0;
  end = consensus.length;
  for (pos = 0; pos < consensus.length and conToRepeat[pos] == -1; ++pos);
  start = pos;
  for (pos = consensus.length-1; pos >= start and conToRepeat[pos] == -1; --pos);
  end = pos + 1;
	 
  
  // do something with the forward stranded sequence
  for (pos = start; pos < end; pos++) {
    if (conToRepeat[pos] == -1) {
      // indel
      profileCount.Set(5, pos);
    }
    else {
      if (orientation == 1) {
	repeatPos = conToRepeat[pos];
	assert(repeatPos >= 0 and repeatPos < repeatSeq.length);
	repeatChar = nuc_index[repeatSeq.seq[repeatPos]];
	assert(repeatChar >= 0 and repeatChar <= 4);
	profileCount.Set(repeatChar, pos);
      }
      else {
	repeatPos = conToRepeat[pos];
	assert(repeatPos >= 0 and repeatPos < repeatSeq.length);
	repeatChar = comp_bin[nuc_index[repeatSeq.seq[repeatSeq.length - repeatPos - 1]]]`1;
	assert(repeatChar >= 0 and repeatChar <= 4);
	profileCount.Set(repeatChar, pos);
      }
    }
  }
}

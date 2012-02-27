/***************************************************************************
 * Title:          ContigMap.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "ContigMap.h"


ssize_t ContigMap::FindInterval(ssize_t pos, ssize_t &intervalIndex, ssize_t start, ssize_t end) {
  assert(startPositions.size() == endPositions.size());

  if (startPositions.size() == 0) {
    intervalIndex = -1;
    return 1;
  }
  
  if (start == -1 or end == -1) {
    start = 0;
    end   = startPositions.size()-1;
  }
  
  if (pos < startPositions[start] or pos > endPositions[end]) {
    // Did not find the interval anywhere in boundaries
    intervalIndex = -1;
    return 0;
  }
  if (pos >= startPositions[start] and pos <= endPositions[end]) {
    // Found the interval, return it.
    intervalIndex = start;
		return 1;
  }
  else {
    // Binary search for the interval.
    ssize_t middle = (start + end) / 2;
    if (pos <= endPositions[middle]) {
      return FindInterval(pos, intervalIndex, start, middle);
    }
    else {
      return FindInterval(pos, intervalIndex, middle, end);
    }
  }

	assert(0); // should never reach here
}


//void ContigMap::CreateIntervals(SeqVector &sequences, int gap) {
//  int i;
//  int startPos, endPos;
//  startPos = 0;
//  for (i = 0; i < sequences.size(); i++) {
//    startPositions.push_back(startPos);
//    endPositions.push_back(endPos);
//    startPos = endPos + gap;
//    endPos   = startPos + sequences[i]->length;
//  }
//}
    
    

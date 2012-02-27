/***************************************************************************
 * Title:          ContigMap.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef CONTIG_MAP_H_
#define CONTIG_MAP_H_

#include <vector>

#include "DNASequence.h"
#include "mctypes.h"

//typedef std::vector<DNASequence*> SeqVector;

class ContigMap {
public:
  IntVector startPositions, endPositions;
  ssize_t size() { return startPositions.size(); }
  ssize_t FindInterval(ssize_t pos, ssize_t &intervalIndex, ssize_t start = -1, ssize_t end = -1);
  //  void CreateIntervals(SeqVector &sequences, int gap);
};


#endif

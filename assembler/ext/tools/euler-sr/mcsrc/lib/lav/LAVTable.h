/***************************************************************************
 * Title:          LAVTable.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef LAV_TABLE_H_
#define LAV_TABLE_H_

#include <string>
#include <vector>
#include "DNASequence.h"
#include "mctypes.h"

class LAVRow {
public:
  ssize_t tStart;
  ssize_t tEnd;
  ssize_t qStart;
  ssize_t qEnd;
  ssize_t strand;
  ssize_t chainId;
  double ComputeScore(DNASequence &t, DNASequence &q, 
		     ssize_t earlyStart, ssize_t earlyEnd,
		     FloatMatrix &scoreMat, ssize_t &lengthAligned, 
		     ssize_t &numMatches, ssize_t &numMismatches);
  LAVRow() {
    tStart = -1;
    tEnd   = -1;
    qStart = -1;
    qEnd   = -1;
    strand = -1;
    chainId= -1;
  }
};

double ScoreRows(DNASequence &t, DNASequence &q, 
		std::vector<LAVRow> &rows,
		ssize_t tStart, ssize_t tEnd, ssize_t &lengthAligned,
		ssize_t &numMatches, ssize_t &numMismatches,
		FloatMatrix &scoreMat);


double ScoreRows(DNASequence &t, DNASequence &q, 
		std::vector<LAVRow> &rows,
		ssize_t tStart, ssize_t tEnd,
		FloatMatrix &scoreMat, 
		ssize_t &lengthAligned, ssize_t &numMatches, ssize_t &numMismatches,
		double &percentIdentity);
  
#endif

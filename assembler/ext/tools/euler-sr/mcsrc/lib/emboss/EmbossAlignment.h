/***************************************************************************
 * Title:          EmbossAlignment.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef EMBOSS_ALIGNMENT_H_ 
#define EMBOSS_ALIGNMENT_H_

#include <stdlib.h>
#include "DNASequence.h"
#include "mctypes.h"

class EmbossAlignment {
public:
  ssize_t    length;
  ssize_t    identEq, identTot;
  double identity;
  ssize_t    simEq, simTot;
  double similarity;
  double pctIdentity;
  ssize_t    gapEq, gapTot;
  double gaps;
  double alignScore;
  
  ssize_t *locations;

  ssize_t refStart, qryStart;
  EmbossAlignment() {
    locations = NULL;
    refStart = -1;
    qryStart = -1;
  }

  void CalculateAverageScore(ssize_t window, double *&scores, ssize_t &scoreLength,
			     FloatMatrix &scoreMat, 
			     double gapOpen, double gapExtend,
			     DNASequence &refSeq, DNASequence &qrySeq,
			     ssize_t *markers, ssize_t numMarkers);
};
#endif


/***************************************************************************
 * Title:          InversionUtils.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "SeqUtils.h"
#include "InversionUtils.h"
#include "InversionBins.h"

ssize_t CheckForDeletion(ssize_t tStart, ssize_t tEnd, ssize_t qStart, ssize_t qEnd, double maxRat) {
  double tNumRat, qNumRat;

  if (tEnd < 0 or tStart < 0)
    return 0;
  tNumRat = (double(tEnd - tStart) ) / (qEnd - qStart);
  qNumRat = (double(qEnd - qStart) )  / (tEnd - tStart);
  if ( tNumRat > 10.0  or qNumRat > 10 ) {
    // This is likely a deletion. 
    // maybe for now I should just output the ratios to check them out.
    return 1;
  }
  else 
    return 0;
}


void DoLocalAlign(DNASequence &refSeq, DNASequence &qrySeq,
		  ssize_t invStart, ssize_t  invEnd, ssize_t qryStart, ssize_t qryEnd, 
		  ssize_t &type) {

  std::string alignCommand;
  DNASequence invRCSeq, // the forward sequence of the
    // inversion (reversed complement)
    invSeq,  // the inverted sequence
    qryGapSeq;  // the unaligned sequence in query.
		
  double forwardScore, reverseScore;
  invRCSeq._ascii  = 1;
  invSeq._ascii    = 1;
  qryGapSeq._ascii = 1;
  // Prepare the inversion sequence
  invRCSeq.seq = &refSeq.seq[invStart];
  invRCSeq.length = invEnd - invStart + 1;
  *invRCSeq.titlestream << "ref_strand ";
  MakeRC(invRCSeq, invSeq);
  *invSeq.titlestream << "inversion";
	      
  // Prepare the query gap sequence
  qryGapSeq.seq = &qrySeq.seq[qryStart];
  qryGapSeq.length = qryEnd - qryStart + 1;
  *qryGapSeq.titlestream << "query_gap";

  alignCommand = "~/projects/sw/bin/needle -aformat srspair -gapOpen 10.0 -gapExtend 0.5 ";

  EmbossAlignment alignment;
  reverseScore     = -99999;
  forwardScore = -99999;
  if (double(invRCSeq.length) * qryGapSeq.length >= 10000000.0) {
    type = DELETED; //foundWith = 0; foundWithout = 0; foundUnknown = 1;
  }
  else {
		  
    if (!EmbossAlign(alignCommand, invSeq, qryGapSeq, alignment) ) {
      // Todo: handle this error, maybe with unknown
      std::cout <<" error aligning inversion " << std::endl;
      exit(1);
    }
    else {
      reverseScore = alignment.alignScore;
    }
		  
    *qryGapSeq.titlestream << "query_gap";
    alignCommand = "~/projects/sw/bin/needle -aformat srspair -gapOpen 10.0 -gapExtend 0.5 ";

    if (!EmbossAlign(alignCommand, invRCSeq, qryGapSeq, alignment) ) {
      // Todo: handle this error, maybe with unknown
      std::cout << "error aligning inversion rc " << std::endl;
      exit(1);
    }
    else {
      forwardScore = alignment.alignScore;
    }
  }
  double randAlignSlope = 0.9;
  double randAlignYIsect = 55;
  double randStdDevSlope = 0.15;
  double randStdDevYIsect = 12;
  double randomScore;
  // Now how to determine if one alignment is much better than the other alignment? 
  double randomAlignScore = randAlignYIsect + randAlignSlope * invSeq.length;
  double randomAlignStdDev = randStdDevSlope * invSeq.length + randStdDevYIsect;
  randomScore = randomAlignScore +2* randomAlignStdDev ;
  //  std::cout<< "got scores: " << reverseScore << " " << forwardScore << " " << randomScore << std::endl;
  if ( reverseScore > randomScore) {
    // We can definitely say that there was an inversion here that we didn't detect 
    // with blastz.
    type = WITH_INVERSION; //foundWith = 1; foundWithout = 0; foundUnknown = 0;
  }
  else {
    if (forwardScore >= randomScore) {
      // Found a good alignment in the forward direction, so
      // there is no inversion here.
      type = WITHOUT_INVERSION; // foundWith = 0; foundWithout = 1; foundUnknown = 0;
    } 
    else {
      // There is no positive reverse or forward alignment for
      // this inversion save it as an unknown score.
      type = DELETED; // foundWith = 0; foundWithout = 0; foundUnknown = 1;
    }
  }
}

void ReadValidatedFile(std::string &fileName, StringVector &species, InversionList &invList) {
  std::cout << "this function has been deprecated " << std::endl;
  assert(0);
}



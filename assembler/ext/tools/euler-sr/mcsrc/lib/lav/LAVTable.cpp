/***************************************************************************
 * Title:          LAVTable.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "LAVTable.h"
#include "SeqUtils.h"
#include "align/alignutils.h"

double LAVRow::ComputeScore(DNASequence &tSeq, DNASequence &qSeq, 
			   ssize_t earlyStart, ssize_t earlyEnd,
			   FloatMatrix &scoreMat, ssize_t &lengthAligned, ssize_t &numMatches, ssize_t &numMismatches) { 
  // The alignment is contiguous, so just score each nucleotide
  //UNUSED+// ssize_t  i,dir;
  ssize_t tStartPos, qStartPos ;
  ssize_t tPos, qPos, tEndPos, qEndPos; // qEndPos doesn't need to be used..
  double score = 0;
  // {t,q}StartPos are half-open 0 based intervals, this 
  // is a bit of a pain, but the indexing starts at 1 in the blastz
  // output.
  // So everything has to be decreased by 1 here .
  numMatches = 0;
  numMismatches = 0;
  tStartPos = tStart + earlyStart - 1;
  qStartPos = qStart + earlyStart - 1;
  tEndPos   = tEnd - earlyEnd - 1;
  qEndPos   = qEnd - earlyEnd - 1;
  lengthAligned = tEndPos - tStartPos + 1;
  char t, q;
  //UNUSED// char nucs[] = {'G', 'A', 'C', 'T', 'N'};

  // Print the sequences 
  DNASequence tSeqFrag, qSeqFrag, qrc;
  tSeqFrag._ascii = 1;
  qSeqFrag._ascii = 1;
  /*
      tSeqFrag.seq = &tSeq.seq[tStartPos];
      tSeqFrag.length = tEndPos - tStartPos + 1;
      
      qSeqFrag.seq = &qSeq.seq[qStartPos];
      qSeqFrag.length = qEndPos - qStartPos + 1;
      tSeqFrag.PrintSeq(std::cout); std::cout << std::endl;
      if (strand == 0) {
      qSeqFrag.PrintSeq(std::cout); std::cout << std::endl;
      }
      else {
      qSeqFrag.seq = &qSeq.seq[qSeq.length - qEndPos - 1];
    qrc._ascii = 1;
    MakeRC(qSeqFrag, qrc);
    qrc.PrintSeq(std::cout); std::cout << std::endl;
  }
  */
  unsigned char tc, qc;
  if (strand == 0) {
    // Score according to the forward strand
    // tPos, qPos are already set
    for (tPos = tStartPos, qPos = qStartPos; 
	 tPos < tEndPos; tPos++, qPos++) {
      tc = toupper(tSeq.seq[tPos]);
      qc = toupper(qSeq.seq[qPos]);
      //      std::cout << tc << ", " << qc;
      t = nuc_index[tc];
      q = nuc_index[qc];
      //      std::cout << ", " << scoreMat[t][q];
      if (t == q)
	++numMatches;
      else
	++numMismatches;
      score += scoreMat[t][q];
    }
    //    std::cout << std::endl;
  }
  else { 
    qStartPos = qSeq.length - (qStartPos) - 1; 

    for (tPos = tStartPos,
	   qPos = qStartPos; tPos < (tEndPos - earlyEnd-1); tPos++, qPos--) {
      tc = toupper(tSeq.seq[tPos]);
      qc = comp_ascii[toupper(qSeq.seq[qPos])];
      t = nuc_index[tc];
      q = nuc_index[qc];
      if (t == q) 
	++numMatches;
      else
	++numMismatches;
      score += scoreMat[t][q];
    }
  }
  return score;
}


double ScoreRows(DNASequence &tSeq, DNASequence &qSeq, 
		std::vector<LAVRow> &rows,
		ssize_t tStart, ssize_t tEnd,
		FloatMatrix &scoreMat, 
		ssize_t &lengthAligned,  ssize_t &numMatches, ssize_t &numMismatches,
		double &percentIdentity) {

  
  // Score a series of rows, possibly starting in the middle of 
  // the first row, and before the end of the last row.
  ssize_t earlyStart;
  ssize_t earlyEnd;
  ssize_t gapOpen, gapExtend;
  double score, regionScore;
  gapOpen = 300;
  gapExtend = 10;
  score = 0;
  numMismatches = 0;
  numMatches = 0;
  if (rows.size() == 0) 
    return 0;

  if (tStart > 0 && tStart > rows[0].tStart) {
    earlyStart = tStart - rows[0].tStart + 1;
  }
  else {
    earlyStart = 0;
  }
  if (tStart >= 0 and tStart < rows[0].tStart)
    score += gapOpen + ((rows[0].tStart - tStart -1) * gapExtend);

  ssize_t lr = rows.size()- 1;
  if (tEnd > 0 && tEnd < rows[lr].tEnd) {
    earlyEnd = rows[lr].tEnd  - tEnd;
  }
  else {
    earlyEnd = 0;
  }

  if (tEnd > rows[lr].tEnd) 
    score += gapOpen + (tEnd - rows[lr].tEnd-1)*gapExtend;
  ssize_t numRowMatched, numRowMismatched = 0;
  ssize_t rowLength;
  lengthAligned = 0;
  if (rows.size() == 1) {
    score += rows[0].ComputeScore(tSeq, qSeq, earlyStart, earlyEnd, scoreMat, 
				  lengthAligned,  numRowMatched, numRowMismatched);
    percentIdentity = double(numRowMatched) / lengthAligned;
    numMatches = numRowMatched;
    numMismatches = numRowMismatched;
    return score;
  }
  numMatches = 0; 
  numMismatches = 0;
  // Else, there are multiple rows.
  score += rows[0].ComputeScore(tSeq, qSeq, earlyStart, 0, scoreMat, 
				rowLength, numRowMatched, numRowMismatched); 
  lengthAligned += rowLength;
  ssize_t i;
  ssize_t gapLen;
  for (i = 1; i < rows.size()-1; i++) {
    regionScore = rows[i].ComputeScore(tSeq, qSeq, 0, 0, scoreMat, 
				       rowLength, numRowMatched, numRowMismatched);
    lengthAligned += rowLength;
    score += regionScore;
    gapLen = std::max(rows[i].tStart - rows[i-1].tEnd -2, rows[i].qStart - rows[i-1].qEnd - 2);
    score += gapOpen + gapExtend*gapLen;
    numMatches += numRowMatched;
    numMismatches += numRowMismatched;
  }

  regionScore = rows[i].ComputeScore(tSeq, qSeq, 0, earlyEnd, scoreMat, 
				     rowLength, numRowMatched, numRowMismatched);
  lengthAligned += rowLength;
  numMatches += numRowMatched;
  numMismatches += numRowMismatched;
  score += regionScore;
  //  std::cout << "last: " << score << std::endl;
  gapLen = std::max(rows[i].tStart - rows[i-1].tEnd -2, rows[i].qStart - rows[i-1].qEnd - 2);
  score += gapOpen + gapExtend*gapLen;
  percentIdentity = double(numMatches) / rowLength;
  return score;
}

double ScoreRows(DNASequence &tSeq, DNASequence &qSeq, 
		std::vector<LAVRow> &rows,
		ssize_t tStart, ssize_t tEnd,
		ssize_t &lengthAligned,
		ssize_t &numMatches, ssize_t &numMismatches,
		FloatMatrix &scoreMat) {
  double percentIdentity;
  return ScoreRows(tSeq, qSeq, rows, tStart, tEnd, scoreMat, 
		   lengthAligned, numMatches, numMismatches, 
		   percentIdentity);
}

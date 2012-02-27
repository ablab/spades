/***************************************************************************
 * Title:          AlignmentPrinter.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  09/05/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "AlignmentPrinter.h"

#include "DNASequence.h"
#include "SeqUtils.h"
#include "align/alignutils.h"

char GetAlignChar(unsigned char refChar, unsigned char qryChar) {
  if (nuc_index[refChar] < 4 and
      nuc_index[qryChar] < 4 and
      nuc_index[refChar] == nuc_index[qryChar]) {
    return '|';
  }
  else {
    return ' ';
  }
}


void PrintAlignment(DNASequence &refSeq, 
		    DNASequence &qrySeq,
		    ssize_t refStartPos, ssize_t qryStartPos,
		    ssize_t *alignment, ssize_t alignmentLength,
		    std::ostream &out, ssize_t lineLength ) {
  std::string refAlignString, qryAlignString, alignString;
  ssize_t refAlignStart, qryAlignStart;
  MakeAlignmentString(refSeq, qrySeq, refStartPos, qryStartPos, 
		      alignment, alignmentLength,
		      refAlignString, qryAlignString, alignString, 
		      refAlignStart, qryAlignStart);

  PrintAlignmentString(refAlignString, qryAlignString, alignString,
		       refAlignStart, qryAlignStart, out, lineLength);
  out << std::endl;
}

void MakeAlignmentString(DNASequence &refSeq, 
			 DNASequence &qrySeq,
			 ssize_t refStartPos, ssize_t qryStartPos,
			 ssize_t *alignment, ssize_t alignmentLength,
			 std::string &refAlignString, 
			 std::string &qryAlignString, 
			 std::string &alignString,
			 ssize_t &refAlignStart, ssize_t &qryAlignStart) {

  // Find the first aligned position in locations.
  //UNUSED+// ssize_t  r,p;
  ssize_t  a;
  ssize_t alignStartIndex, alignEndIndex;
  for (a = 0; a < alignmentLength and alignment[a] == -1; a++) ;
  alignStartIndex = a;
  for (a = alignmentLength-1; a >= 0 and alignment[a] == -1; --a) ;
  alignEndIndex = a;

  std::string refStr, qryStr, alignStr;
  ssize_t refPos, qryPos;
  a = alignStartIndex;
  refPos = refStartPos + alignStartIndex;
  qryPos = qryStartPos + alignment[a];
  refAlignStart = a;
  qryAlignStart = qryPos;
  unsigned char refNuc, qryNuc;
  refAlignString = qryAlignString = alignString = "";
  while (a <= alignEndIndex) {
    refNuc = (unsigned char) refSeq.seq[refPos];
    qryNuc = (unsigned char) qrySeq.seq[qryPos];
    refAlignString += nuc_char[refNuc];
    qryAlignString += nuc_char[qryNuc];
    alignString    += GetAlignChar(refNuc, qryNuc);
    refPos++;
    qryPos++;
    a++;
    while (a < alignEndIndex and alignment[a] == -1) {
      refPos++;
      a++;
      refNuc = (unsigned char) refSeq.seq[refPos];
      refAlignString += nuc_char[refNuc];
      alignString += ' ';
      qryAlignString += '-';
    }
    if (a < (alignmentLength-1)) {
      ssize_t g;
      for (g = 0; g < (alignment[a+1] - alignment[a] - 1); ++g) {
				qryPos++;
				qryNuc = (unsigned char) qrySeq.seq[qryPos];
				refAlignString += '-';
				alignString += ' ';
				qryAlignString += nuc_char[qryNuc];
      }
    }
  }
}


void PrintAlignmentString(std::string refAlignString, 
			  std::string qryAlignString,
			  std::string alignString,
			  ssize_t refPos, ssize_t qryPos,
			  std::ostream &out,
			  ssize_t lineLength) {
  assert(refAlignString.size() == qryAlignString.size());
  assert(qryAlignString.size() == alignString.size());

  ssize_t linePos, refStringPos, qryStringPos, alignStringPos;
  //UNUSED// ssize_t i;
  refStringPos= qryStringPos = alignStringPos = 0;
  while (refStringPos < refAlignString.size()) {
    out << "ref: "; out.width(8); out << refPos << " ";
    for (linePos = 0; linePos < lineLength; linePos++) {
      if (refStringPos >= refAlignString.size()) 
	break;
      out << refAlignString[refStringPos];
      refStringPos++;
      ++refPos;
    }
    out << std::endl;
    out << "              ";
    for (linePos = 0; linePos < lineLength; linePos++) {
      if (alignStringPos >= alignString.size())
	break;
      out << alignString[alignStringPos];
      alignStringPos++;
    }
    out << std::endl;
    out << "qry: "; out.width(8); out << qryPos << " ";
    for (linePos = 0; linePos < lineLength; linePos++) {
      if (qryStringPos >= qryAlignString.size())
	break;
      out << qryAlignString[qryStringPos];
      qryStringPos++;
      ++qryPos;
    }
    out << std::endl << std::endl;

  }
}

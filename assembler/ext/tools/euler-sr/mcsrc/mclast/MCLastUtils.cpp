/***************************************************************************
 * Title:          MCLastUtils.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "MCLastUtils.h"
#include "MCHsp.h"
#include "align/alignutils.h"
#include "align/graphalign.h"
#include "AlignmentPrinter.h"
#include "utils.h"
#include "mctypes.h"

ssize_t ParseHashOptions(int argc, char**argv,
		     HashPattern &hashPattern, 
		     ssize_t &mask) {
  int argi = 0;
  std::string hashDescriptor = "";
  while (argi < argc) {
    if (strcmp(argv[argi], "-n") == 0) {
      ++argi;
      HashValue::indexSize = atoi(argv[argi]);
      HashValue::convertSize = HashValue::size - HashValue::indexSize;
    }
    if (strcmp(argv[argi], "-maskPattern") == 0) {
      ++argi;
      hashDescriptor = argv[argi];
    }
    if (strcmp(argv[argi], "-maskFile") == 0) {
      ++argi;
      hashDescriptor = argv[argi];
      ReadHashPattern(hashDescriptor, hashPattern);
    }
    if (strcmp(argv[argi], "-im") == 0) {
      mask = 0;
    }
    ++argi;
  }
  if (hashDescriptor == "")
    InitHashPattern(hashPattern);
  else 
    ReadHashPatternString(hashDescriptor, hashPattern);
  return 0;
}



void ReadHashPatternString(std::string &hashPatternStr,
			   HashPattern &hashPattern) {
  ssize_t i;
  const char *str = hashPatternStr.c_str();
  for (i = 0; i < hashPatternStr.size(); i++) {
    if (str[i] == '0')
      hashPattern.mask.push_back(0);
    else if (str[i] == '1')
      hashPattern.mask.push_back(1);
  }
  for (i = 0; i < hashPattern.mask.size(); i++ ) {
    if (hashPattern.mask[i])
      hashPattern.maskIndices.push_back(i);
  }
  hashPattern.SetNChars();
}

void ReadHashPattern(std::string &patFileName,
		     HashPattern &hashPattern) {
  std::ifstream patIn;
  openck(patFileName, patIn);
  char c;
  while (patIn) {
    c = patIn.get();
    if (c == '0') 
      hashPattern.mask.push_back(0);
    if (c == '1')
      hashPattern.mask.push_back(1);
  }
  hashPattern.SetNChars();
}

double ExtendMatch(DNASequence &refSeq, ssize_t refPos,
		  DNASequence &qrySeq, ssize_t qryPos,
		  HashPattern &hashPattern,
		  FloatMatrix &scoreMat, double gapOpen, double gapExtend,
		  MCHsp &hsp,
		  ssize_t *&alignment, ssize_t &alignmentLength) {

  ssize_t refStartPos, qryStartPos, refEndPos, qryEndPos;

  refStartPos = refPos; refEndPos = refPos + hashPattern.size();
  qryStartPos = qryPos; qryEndPos = qryPos + hashPattern.size();
  
  ssize_t refAlignStart, qryAlignStart;
  double alignScore;
  alignScore = AlignRegion(refSeq, refStartPos, refEndPos, 
			   qrySeq, qryStartPos, qryEndPos,
			   scoreMat, gapOpen, gapExtend,
			   refAlignStart, qryAlignStart,
			   alignment, alignmentLength);


  hsp.refStart = refAlignStart;
  hsp.qryStart = qryAlignStart;

  hsp.refEnd   = refAlignStart + alignmentLength;
  hsp.qryEnd   = alignment[alignmentLength-1];

  hsp.score = alignScore;
  return alignScore;
}
		 
ssize_t FindHsps(DNASequence &refSeq, DNASequence &qrySeq,
	     HashTable &refHashTable, HashTable &qryHashTable,
	     ssize_t window,
	     HashPattern &hashPattern,
	     FloatMatrix &scoreMat, double gapOpen, double gapExtend,
	     ssize_t maxScore, ssize_t maxDistance,
	     std::vector<MCHsp> &hsps) {
  ssize_t pos;
  HashValue hashValue;
  ssize_t lastPos = refSeq.length - hashPattern.size();
  Chain *chain;
  MCHsp hsp;
  ssize_t *alignment, alignmentLength;
  MCHspIntervalList hspIntervals;
  for (pos = 0; pos < lastPos; pos++) {
    if (GetHashValue(refSeq, pos, hashPattern, 1, hashValue)) {
      // Try to find the reference in query
      qryHashTable.Find(hashValue, chain);
      if (chain != NULL) {
	ssize_t c;
	// Extend each match.
	ssize_t forPos;
	std::cout<<" found hash match " << std::endl;
	for (c = 0; c < chain->numAlign(); c++) {
	  forPos = qrySeq.GetPos(chain->
													 DNASequence::FORWARD_STRAND);
	  std::cout << "extending " << pos << " " << chain->positions[c].pos << std::endl;
	  if (abs(pos - forPos ) < maxDistance) {
	    ExtendMatch(refSeq, pos, 
			qrySeq, chain->positions[c].pos, 
			hashPattern,
			scoreMat, gapOpen, gapExtend,
			hsp, alignment, alignmentLength);
	    if (hsp.score < maxScore) {
	      MCHsp *hspPtr = new MCHsp(hsp);
	      MCHsp *replaced;
	      replaced = hspIntervals.Insert(hspPtr, 1.2);
	      
	    }
	  }
	}
      }
    }
  }
  ssize_t i;
  for (i = 0; i < hspIntervals.hsps.size(); i++) {
    hsps.push_back(*hspIntervals.hsps[i]);
  }
}


ssize_t FindCommonK(DNASequence &refSeq, DNASequence &qrySeq,
		HashTable &refHashTable, HashTable &qryHashTable,
		HashPattern &hashPattern, std::vector<ssize_t> & refPos, std::vector<ssize_t> &qryPos,
		ssize_t maskRepeats) {
  ssize_t pos;
  HashValue hashValue;
  ssize_t lastPos = refSeq.length - hashPattern.size();
  Chain *chain;
  ssize_t *alignment, alignmentLength;
  for (pos = 0; pos < lastPos; pos++) {
    if (GetHashValue(refSeq, pos, hashPattern, maskRepeats, hashValue)) {
      // Try to find the reference in query
      qryHashTable.Find(hashValue, chain);
      //      if (chain != NULL and chain->positions.size() == 1) {
      if (chain != NULL) {
	ssize_t forPos;
	ssize_t index;
	for (index = 0; index < chain->positions.size(); index++) {
	  forPos = qrySeq.GetPos(chain->positions[index].pos,
				 DNASequence::FORWARD_STRAND);
	  refPos.push_back(pos);
	  qryPos.push_back(forPos);
	}
      }
    }
  }
}



void InitHashPattern(std::string fileName, HashPattern &hashPattern) {
  std::ifstream in;
  openck(fileName, in);
  char c;
  while ( (c=in.get()) and c != EOF and c != '\n') {
    if (c == '1') {
      hashPattern.mask.push_back(1);
    }
    else if (c == '0') {
      hashPattern.mask.push_back(0);
    }
  }
  hashPattern.SetNChars();
}


void InitHashPattern(HashPattern &hashPattern) {
  hashPattern.mask.push_back(1);
  hashPattern.mask.push_back(1);
  hashPattern.mask.push_back(0);
  hashPattern.mask.push_back(1);
  hashPattern.mask.push_back(1);
  hashPattern.mask.push_back(0);
  hashPattern.mask.push_back(1);
  hashPattern.mask.push_back(1);
  hashPattern.mask.push_back(0);
  hashPattern.mask.push_back(1);
  hashPattern.mask.push_back(1);
  hashPattern.mask.push_back(0);
  hashPattern.mask.push_back(1);
  hashPattern.mask.push_back(1);
  hashPattern.mask.push_back(0);
  hashPattern.mask.push_back(1);
  hashPattern.mask.push_back(1);
  hashPattern.mask.push_back(0);
  hashPattern.mask.push_back(1);
  hashPattern.SetNChars();
}

/***************************************************************************
 * Title:          MaskRepeats.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
// std includes
#include <ios>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <unistd.h>

// my stuff
#include "utils.h"
#include "lav/LAVUtils.h"
#include "lav/LAVReader.h"
#include "lav/LAVTable.h"
#include "lav/LAVFile.h"
#include "lav/LAVBlock.h"
#include "lav/LAVAlignedContig.h"
#include "DNASequence.h"
#include "SeqReader.h"
#include "align/alignutils.h"
#include "SeqUtils.h"


int main(int argc, char* argv[]) {

  if (argc < 5) { 
    std::cout << "usage: " << argv[0] 
	      << " lavfile refseq qryseq refOut qryOut [-t threshold][-w windowSize]"
	      << std::endl;
    exit(0);
  }
  std::string lavFileName = argv[1];
  std::string refSeqName  = argv[2];
  std::string qrySeqName  = argv[3];
  std::string refOutName  = argv[4];
  std::string qryOutName  = argv[5];

  std::string opt;
  ssize_t windowSize = 15;
  ssize_t threshold = 10;
  int argi = 5;
  while (argi < argc) {
    argi = 5;
    opt = argv[argi++];
    if (opt == "-w") {
      windowSize = atoi(argv[argi]);
      ++argi;
    }
    if (opt == "-t") {
      threshold = atoi(argv[argi]);
      ++argi;
    }
  }
  LAVFile lavFile;
  LAVReader::ReadLAVFile(lavFileName, lavFile);

  DNASequence refSeq, qrySeq, qrySeqRC, maskedRefSeq, maskedQrySeq;
  
  SeqReader::GetSeq(refSeqName, refSeq, SeqReader::noConvert);
  SeqReader::GetSeq(qrySeqName, qrySeq, SeqReader::noConvert);
  
  MakeRC(qrySeq, qrySeqRC);
  maskedRefSeq = refSeq;

  FloatVector refSeqScores;
  IntVector   refSeqMult, qrySeqMult;
  refSeqMult.resize(refSeq.length);
  qrySeqMult.resize(qrySeq.length);
  ssize_t i;
  for (i = 0; i < refSeq.length; i++)
    refSeqMult[i]   = 0;

  for (i = 0; i < qrySeq.length; i++)
    qrySeqMult[i] = 0;

  LAVAlignedContig *alignedContig;
  LAVBlock *block;
  std::vector<ssize_t> rcBegin, rcEnd;
  double forwardScore = 0.0;
  double reverseScore = 0.0;
  // Find out where all of the 
  FloatVector alignScore;
  ssize_t blockLength;
  
  ssize_t ac,b,bi;

  ssize_t pos;
  
  // Count multiplicities
  for (ac = 0; ac < lavFile.alignments.size(); ac++) {
    alignedContig = lavFile.alignments[ac];
    for (b = 0; b < alignedContig->size(); b++) {
      block = alignedContig->alignments[b];
      for (bi = 0; bi < block->size(); bi++) {
	for (pos = block->refALBegin[bi];
	     pos < block->refALEnd[bi]; pos++) {
	  refSeqMult[pos]++;
	}
	for (pos = block->qryALBegin[bi];
	     pos < block->qryALEnd[bi]; pos++) {
	  if (alignedContig->qryContig.strand == 0) 
	    qrySeqMult[pos]++;
	  else
	    qrySeqMult[qrySeq.length - pos + 1]++;
	}
      }
    }
  }

  ssize_t r;
  for (r = 0; r < refSeq.length; r++ ) 
    if (refSeqMult[r] > threshold) {
      refSeq.MaskPosition(r);
    }
  ssize_t q;
  for (q = 0; q < qrySeq.length; q++ ) 
    if (qrySeqMult[q] > threshold) {
      qrySeq.MaskPosition(q);
    }
  
  

  std::ofstream refOut, qryOut;
  openck(refOutName, refOut, std::ios::out);
  refSeq.PrintSeq(refOut);
  refOut << std::endl;
  refOut.close();


  openck(qryOutName, qryOut, std::ios::out);
  qrySeq.PrintSeq(qryOut);
  qryOut << std::endl;
  qryOut.close();

  return 0;
}

/***************************************************************************
 * Title:          InvRepFinder.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <vector>
#include <iostream>
#include <string>
#include <limits.h>
// defined in common
#include "DNASequence.h"
#include "TupleLib.h"
#include "SeqReader.h"
#include "utils.h"
#include "alignutils.h"
#include "SeqUtils.h"
#include "mctypes.h"

// defined here
#include "hash/HashUtils.h"
#include "MCLastUtils.h"
#include "MCHsp.h"

void InitHashPattern(HashPattern &hashPattern);


int main(int argc, char* argv[]) {

  if (argc < 2) {
    std::cout << "usage: invrepfind seq1 [-n nbins] [-im ] [-m maxScore] [-maskFile mask_file "
	      << "[-d distancethresh]" << std::endl;
    std::cout << "do hashing on 13 / 19 nucleotides " << std::endl;
    std::cout << "nbins are the number of nucleotides to use in a bin (8) " 
	      << std::endl;
    std::cout << "-im ignores mask " << std::endl;
    exit(0);
  }
  ssize_t nBins = 8;
  std::string refSeqFileName;
  refSeqFileName = argv[1];
  int argi = 2;
  ssize_t mask = 1;
  ssize_t maxScore = -1000;
  ssize_t maxDistance = 10000;
  std::string hashDescriptor = "";
  HashPattern hashPattern;
  while (argi < argc) {
    if (strcmp(argv[argi], "-n") == 0) {
      ++argi;
      HashValue::indexSize = atoi(argv[argi]);
      HashValue::convertSize = HashValue::size - HashValue::indexSize;
    }
    if (strcmp(argv[argi], "-maskFile") == 0) {
      ++argi;
      hashDescriptor = argv[argi];
      ReadHashPattern(hashDescriptor, hashPattern);
    }
    if (strcmp(argv[argi], "-im") == 0) {
      mask = 0;
    }
    if (strcmp(argv[argi], "-m") == 0){
      ++argi;
      maxScore = atoi(argv[argi]);
    }
    if (strcmp(argv[argi], "-d") == 0) {
      ++argi;
      maxDistance = atoi(argv[argi]);
    }
    ++argi;
  }

  std::ifstream refIn, qryIn;
  openck(refSeqFileName, refIn, std::ios::binary | std::ios::in);
  DNASequence refSeq, qrySeq;
  SeqReader::MaskRepeats();
  SeqReader::GetSeq(refIn, refSeq);

  MakeRC(refSeq, qrySeq);

    
  if (hashDescriptor == "")
    InitHashPattern(hashPattern);

  ssize_t nBuckets;
  HashTable refHashTable, qryHashTable;
  nBuckets = HashTable::CalcNBuckets(HashValue::indexSize);
  refHashTable.Init(nBuckets);
  qryHashTable.Init(nBuckets);

  HashGenome(refSeq, 1, 
	     hashPattern, mask,
	     refHashTable);

  HashGenome(qrySeq, 1, 
	     hashPattern, mask,
	     qryHashTable);

  std::cout << "done hashing genome " << std::endl;

  // Now find matches in ref/qry
  std::vector<MCHsp> hsps;
  
  double gapOpen, gapExtend;
  std::string scoreMatName;
  FloatMatrix scoreMat;
  char *home = getenv("HOME");
  scoreMatName = std::string(home) + "/projects/mcsrc/align/data/scoremat.txt";
  gapOpen = 400; gapExtend  = 50;
  ReadScoreMatFile(scoreMatName, scoreMat);


  //  QueryHash(refSeq, hashPattern, qryHashTable);
  FindHsps(refSeq, qrySeq, refHashTable, qryHashTable, 
	   0, hashPattern, scoreMat, gapOpen, gapExtend, 
	   maxScore, maxDistance, hsps);

  ssize_t h;
  for (h = 0; h < hsps.size(); h++ ) {
    std::cout << hsps[h].score << " "
	      << hsps[h].refStart << " " << hsps[h].refEnd  << " "
	      << qrySeq.length - hsps[h].qryEnd << " " 
	      << qrySeq.length - hsps[h].qryStart  << std::endl;
  }
  return 0;
}


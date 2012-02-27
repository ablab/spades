/***************************************************************************
 * Title:          DotPlot.cpp 
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

// defined here
#include "hash/HashUtils.h"
#include "MCLastUtils.h"
#include "MCHsp.h"

void InitHashPattern(HashPattern &hashPattern);


int main(int argc, char* argv[]) {

  if (argc < 4) {
    std::cout << "usage: dotplot seq1 seq2 out.dot "
	      << "[-n nbins] [-im ] [-maskFile mask_file] " 
	      << "[ -maskPattern 1011... ] " 
	      << "[ -maskRepeats ] " 
	      << "[-d distancethresh]" << std::endl;
    std::cout << "do hashing on 13 / 19 nucleotides " << std::endl;
    std::cout << "nbins are the number of nucleotides to use in a bin (8) " 
	      << std::endl;
    std::cout << " -maskRepeats masks lower-case nucleotides " << std::endl;
    std::cout << "-im ignores mask " << std::endl;
    exit(0);
  }
  ssize_t nBins = 8;
  std::string refSeqFileName, qrySeqFileName, dotFileName;
  refSeqFileName = argv[1];
  qrySeqFileName = argv[2];
  dotFileName    = argv[3];

  
  int argi = 2;
  ssize_t mask = 0;
  std::string hashDescriptor = "";
  HashPattern hashPattern;
  while (argi < argc) {
    if (strcmp(argv[argi], "-n") == 0) {
      ++argi;
      HashValue::indexSize = atoi(argv[argi]);
      HashValue::convertSize = HashValue::size - HashValue::indexSize;
    }
    if (strcmp(argv[argi], "-maskRepeats") == 0) {
      mask = 1;
    }
    if (strcmp(argv[argi], "-maskPattern") == 0) {
      ++argi;
      hashDescriptor = argv[argi];
      ReadHashPatternString(hashDescriptor, hashPattern);
      std::cout << "got hash pattern " << hashDescriptor << std::endl;
    }
    if (strcmp(argv[argi], "-maskFile") == 0) {
      ++argi;
      hashDescriptor = argv[argi];
      ReadHashPattern(hashDescriptor, hashPattern);
    }
    ++argi;
  }

  std::ifstream refIn, qryIn;
  openck(refSeqFileName, refIn, std::ios::binary | std::ios::in);
  openck(qrySeqFileName, qryIn, std::ios::binary | std::ios::in);

  std::ofstream dotFile;
  openck(dotFileName, dotFile, std::ios::out);

  DNASequence refSeq, qrySeq, qrySeqRC;
  if (mask) 
    SeqReader::MaskRepeats();
  SeqReader::GetSeq(refIn, refSeq, 0);
  SeqReader::GetSeq(qryIn, qrySeq, 0);
  if (hashDescriptor == "" ) {
    std::cout << "using default pattern " << std::endl;
    InitHashPattern(hashPattern);
  }

  ssize_t nBuckets;
  HashTable refHashTable, qryHashTable;
  nBuckets = HashTable::CalcNBuckets(HashValue::indexSize);
  refHashTable.Init(nBuckets);
  qryHashTable.Init(nBuckets);

  HashGenome(refSeq, 1, hashPattern, mask, refHashTable);
  HashGenome(qrySeq, 1, hashPattern, mask, qryHashTable);

  double gapOpen, gapExtend;
  std::string scoreMatName;

  std::vector<ssize_t> refPos, qryPos, qryPosRC;

  FindCommonK(refSeq, qrySeq, refHashTable, qryHashTable, 
	      hashPattern, refPos, qryPos, mask);

  ssize_t h;
  for (h = 0; h < refPos.size(); h++ ) {
    dotFile << refPos[h] << " " << refPos[h] + hashPattern.nChars << " "
	      << qryPos[h] << " " << qryPos[h] + hashPattern.nChars << std::endl;
  }

  MakeRC(qrySeq, qrySeqRC);

  HashTable qryRCHashTable;
  qryRCHashTable.Init(nBuckets);

  HashGenome(qrySeqRC, 1, hashPattern, mask, qryRCHashTable);
  refPos.clear();
  FindCommonK(refSeq, qrySeqRC, refHashTable, qryRCHashTable, 
	      hashPattern, refPos, qryPosRC, mask);

  
  std::cerr << "printing " << refPos.size() << " reverse " << hashPattern.nChars << std::endl;
  for (h = 0; h < refPos.size(); h++ ) {
    dotFile << refPos[h] << " " << refPos[h] + hashPattern.nChars << " "
	    << qryPosRC[h]  << " " 
	    << qryPosRC[h] - hashPattern.nChars << std::endl;
  }
  dotFile.close();
  return 0;
}


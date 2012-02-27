/***************************************************************************
 * Title:          printLavRCCoords.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <ostream>
#include <stdio.h>
#include <unistd.h>


#include "DNASequence.h"
#include "SeqReader.h"
#include "align/alignutils.h"
#include "lav/LAVReader.h"
#include "lav/LAVFile.h"
#include "lav/LAVAlignedContig.h"
#include "lav/LAVBlock.h"


void InitEnv(int argc, char* argv[], std::string &lavFileName, ssize_t &printRef);

void PrintUsage();

int main(int argc, char* argv[]) {
  
  std::string lavFileName;
  DNASequence refSeq, qrySeq; 
  ssize_t printRef;
  printRef = 0;
  InitEnv(argc, argv, lavFileName, printRef);
  // Get input.

  LAVFile lavFile;
  LAVReader::ReadLAVFile(lavFileName, lavFile);

    ssize_t a, b;
  LAVAlignedContig *alignedContig;
  
  DNASequence refSubseq, qrySubseq;

  ssize_t qrySeqLen; 
  for (a = 0; a < lavFile.alignments.size(); a++) {
    alignedContig = lavFile.alignments[a];
    if (alignedContig->qryContig.strand == 1) {
      qrySeqLen = alignedContig->qryContig.end - alignedContig->qryContig.start + 1;
      LAVBlock *block;
      for (b = 0; b < alignedContig->alignments.size(); b++) {
	block = alignedContig->alignments[b];
	if (printRef)
	  std::cout << block->refBegin  << " " << block->refEnd << std::endl;
	else
	  std::cout << qrySeqLen - block->qryEnd << " " << qrySeqLen - block->qryBegin << std::endl;
      }
    }
  }
  return 0;
}



void InitEnv(int argc, char* argv[], 
	     std::string &lavFileName,
	     ssize_t &printRef) {

  ssize_t copt;
  ssize_t i;
  while ( (copt=getopt(argc, argv, "r")) != EOF){
    switch(copt) {
    case 'r':
      printRef = 1;
      continue;
    default:
      PrintUsage();
      exit(1);
    }
  }
  i = optind;
  if (i >= argc) {
    std::cout << "You must specify a lav file. " << std::endl;
    PrintUsage();
    exit(1);
  }
  lavFileName = argv[i];
}


void PrintUsage() {
  std::cout << "plrcc. Print LAVfile reverse-complement coords  . " << std::endl;
  std::cout << "usage: plrcc [-r] lavFile " << std::endl;
  std::cout << "-r  print reference coords " << std::endl;
}

/***************************************************************************
 * Title:          CombineMasking.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  10/24/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <string>
#include "SeqReader.h"
#include "DNASequence.h"
#include "utils.h"

int main(int argc, char *argv[]) {

  if (argc != 4) {
    std::cout << "usage: combmsk seq1 seq2 outfile " << std::endl;
    exit(1);
  }
  std::string outFileName;
  std::string inA, inB;
  inA = argv[1];
  inB = argv[2];
  outFileName= argv[3];

  std::ifstream inAFile, inBFile;
  std::ofstream out;
  
  openck(inA,inAFile);
  openck(inB,inBFile);
  openck(outFileName, out, std::ios::out);

  DNASequence seqA, seqB, joined;

  SeqReader::MaskRepeats();
  SeqReader::GetSeq(inAFile, seqA, SeqReader::noConvert);
  SeqReader::GetSeq(inBFile, seqB, SeqReader::noConvert);

  if (seqA.length != seqB.length) {
    std::cout << seqA.length << "  != " << seqB.length << std::endl;
    exit(0);
  }
  ssize_t i;
  char c;
  for (i = 0; i < seqA.length; i++) {
    c = seqB.seq[i];
    if (c >= 'a' and c <= 'z') {
      seqA.seq[i] = c;
    }
  }
  seqA.PrintSeq(out);
}
  





  
  
    

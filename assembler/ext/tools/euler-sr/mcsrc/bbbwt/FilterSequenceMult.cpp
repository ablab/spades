/***************************************************************************
 * Title:          FilterSequenceMult.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/28/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <time.h>
#include <string>

#include "bbbwt/BBBWTQuery.h"
#include "DNASequence.h"
#include "SeqReader.h"
#include "utils.h"

int main(int argc, char* argv[]) {

  BBBWT csa;
  std::string seqFileName, csaFileName, seqOutFile;
	//UNUSED//  int wordlen;
  if (argc < 5) {
    std::cout << "usage: querySequenceMult reads csa tupleSize minMult seqOut" << std::endl;
    std::cout << "	 Search the multiplicity of each substring of length 'tuplesize' in each"<<std::endl
	      << "   read in 'reads' using the compressed suffix array 'csa'."<<std::endl;
    std::cout << "	 Write to 'seqOut' all sequences that are all above 'minMult' multiplicity" 
	      << std::endl;
    return 1;
  }
  int tupleSize = 0;
  ssize_t minMult   = 0;
  seqFileName = argv[1];
  csaFileName = argv[2];
  tupleSize   = atoi(argv[3]);
  minMult     = atoi(argv[4]);
  seqOutFile  = argv[5];
	
  std::ofstream seqOut;
  openck(seqOutFile, seqOut, std::ios::out);
  BW::Read(csaFileName, csa);
  DNASequence seq, tmp;
  clock_t startTime, curTime;
  startTime = clock();
  double elapsedTime;
  std::ifstream seqIn;
  openck(seqFileName, seqIn, std::ios::in);
	//UNUSED//  int qryVal;
	ssize_t low, high;
  ssize_t seqNumber = 0;
  ssize_t mult;
  ssize_t good;
  std::vector<DNASequence*> goodSequences;
  while(SeqReader::GetSeq(seqIn, seq, SeqReader::noConvert)) {
    ssize_t p;
    if (seqNumber % 1000 == 0) {
      curTime = clock();
      elapsedTime = (double(curTime) - double(startTime)) / CLOCKS_PER_SEC;
      std::cout << seqNumber << " " << elapsedTime << std::endl;
    }
    good = 1;
    for (p = 0; p < seq.length - tupleSize + 1; p++ ) {
      BW::Query(seq, p, tupleSize, csa, low, high);
      mult = high - low;
      if (mult < minMult){
	good = 0;
	break;
      }
    }
    if (good) {
      goodSequences.push_back(new DNASequence(seq));
      //			seq.PrintSeq(seqOut);
      //			seqOut << std::endl;
    }
    ++seqNumber;
  }
  
  ssize_t g;
  for (g =0 ; g < goodSequences.size(); g++) {
    goodSequences[g]->PrintSeq(seqOut);
    seqOut << std::endl;
  }

  return 0;
}

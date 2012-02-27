/***************************************************************************
 * Title:          MatchSequences.cpp 
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
	ssize_t maxQueryLength = -1;
  if (argc < 2) {
    std::cout << "usage: matchSequences reads csa " << std::endl;
    std::cout << "	 print sequenes found in a csa " << std::endl;
    return 1;
  }
	//UNUSED//  int verbose = 0;
  seqFileName = argv[1];
  csaFileName = argv[2];

  BW::Read(csaFileName, csa);
  ssize_t low, high;
	//UNUSED//  int numMatched = 0;
  NamedSequenceList sequences;
	//  ReadNamedSequences(seqFileName, sequences);
	//UNUSED//  int i;
	ssize_t mult;
	std::ifstream seqIn;
	openck(seqFileName, seqIn, std::ios::in);
	DNASequence seq;
	ssize_t seqIndex = 0;
	while(SeqReader::GetSeq(seqIn, seq, SeqReader::noConvert)) {
		if (maxQueryLength != -1 and maxQueryLength < seq.length)
			continue;
    BW::Query(seq, csa, low, high);
    mult = high - low;
    if (mult > 0) {
			seq.PrintlnSeq(std::cout);
    }
		seqIndex++;
		if (seqIndex % 10000 == 0) 
			std::cerr << seqIndex << std::endl;
  }
  
	//  std::cout << numMatched << " " << sequences.size() << std::endl;
  return 0;
}

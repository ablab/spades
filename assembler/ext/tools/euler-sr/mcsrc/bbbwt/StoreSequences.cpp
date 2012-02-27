/***************************************************************************
 * Title:          StoreSequences.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
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
#include "SeqUtils.h"
#include "utils.h"

int main(int argc, char* argv[]) {

  BBBWT csa;
  std::string seqFileName, csaFileName;
  //UNUSED// ssize_t wordlen;
  if (argc < 3) {
    std::cout << "usage: storeReads reads csa [--storeRC]" << std::endl;
		std::cout << "load all reads in 'reads' into csa 'csa'. " << std::endl;
    return 1;
  }
  seqFileName = argv[1];
  csaFileName = argv[2];
  int argi = 3;
  ssize_t storeRC = 0;
  while (argi < argc) {
    if (strcmp(argv[argi], "--storeRC") == 0) {
      storeRC = 1;
    }
    else {
      std::cout << "didn't understand option: " << argv[argi] << std::endl;
      exit(1);
    }
    ++argi;
  }
  
  DNASequence seq, tmp;
  clock_t startTime, curTime;
  startTime = clock();
  double elapsedTime;
  std::ifstream seqIn;
  openck(seqFileName, seqIn, std::ios::in);
  //UNUSED// ssize_t qryVal, low, high;
  ssize_t seqNumber = 0;
  //UNUSED// ssize_t size = 0;
  while(SeqReader::GetSeq(seqIn, seq, SeqReader::noConvert)) {
    //UNUSED// ssize_t p;
    if (seqNumber % 1000 == 0) {
      curTime = clock();
      elapsedTime = (double(curTime) - double(startTime)) / CLOCKS_PER_SEC;
			std::cout << seqNumber << " " << elapsedTime << std::endl;
    }
    BW::Store(seq, csa);
    DNASequence rc;
    MakeRC(seq, rc);
    BW::Store(rc, csa);
    ++seqNumber;
  }
  BW::Write(csa, csaFileName);


  
  return 0;
}

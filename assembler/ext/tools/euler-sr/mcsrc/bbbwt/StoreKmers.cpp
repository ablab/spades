/***************************************************************************
 * Title:          StoreKmers.cpp 
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
#include "utils.h"

int main(int argc, char* argv[]) {

  BBBWT csa;
  std::string seqFileName, csaFileName;
  ssize_t wordlen;
  if (argc < 4) {
    std::cout << "usage: storeKmers infile k outfile" << std::endl;
    return 1;
  }
  seqFileName = argv[1];
  wordlen = atoi(argv[2]);
  csaFileName = argv[3];
  
  DNASequence seq, tmp;
  clock_t startTime, curTime;
  startTime = clock();
  double elapsedTime;
  std::ifstream seqIn;
  openck(seqFileName, seqIn, std::ios::in);
  ssize_t qryVal, low, high;
  ssize_t seqNumber = 0;
  ssize_t size = 0;
  while(SeqReader::GetSeq(seqIn, seq, SeqReader::noConvert)) {
    ssize_t p;
    if (seqNumber % 100 == 0) {
      curTime = clock();
      elapsedTime = (double(curTime) - double(startTime)) / CLOCKS_PER_SEC;
    }
    for (p = 0; p < seq.length - wordlen + 1; p++ ) {
      qryVal = BW::Query(seq, p, wordlen, csa, low, high);
      if (qryVal <= 0) {
				BW::Store(seq, csa, p, wordlen);
				size++;
				std::cout << "size: " << size << std::endl;
      }
      else {
				tmp._ascii = seq._ascii;
				tmp.seq = &(seq.seq[p]);
				tmp.length = wordlen;
				std::cout << "skipped : ";
				tmp.PrintSeq(std::cout);
				std::cout << std::endl;
      }
    }
    BW::Store(seq, csa);
    ++seqNumber;
  }
  std::cout << "stored " << size << std::endl;
  BW::Write(csa, csaFileName);
  // Store the count of the number of tuples stored.  This will help later.
  std::string sizeFile = csaFileName + ".count";
  std::ofstream out;
  openck(csaFileName, out, std::ios::out);
  out << size << std::endl;
  out.close();
  
  return 0;
}

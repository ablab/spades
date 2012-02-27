/***************************************************************************
 * Title:          QuerySequenceMult.cpp 
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
#include "SeqUtils.h"
#include "utils.h"

int main(int argc, char* argv[]) {

  BBBWT csa;
  std::string seqFileName, csaFileName;
	//UNUSED//  int wordlen;
  if (argc < 2) {
    std::cout << "usage: querySequenceMult reads csa k" << std::endl;
		std::cout << "   Search the multiplicity of each substring of length 'k' in each"<<std::endl
							<< "   read in 'reads' using the compressed suffix array 'csa'."<<std::endl;
    return 1;
  }
	//UNUSED//	int tupleSize = 0;
  seqFileName = argv[1];
  csaFileName = argv[2];
	int argi = 3;
	ssize_t printMissing = 0;
	ssize_t printFound = 0;
	std::string missingFileName, foundFileName;
	ssize_t compareToRef = 0;
	std::string refSeqName = "";
	while (argi < argc) {
		if (strcmp(argv[argi], "-printMissing") == 0){ 
			missingFileName = argv[++argi];
			printMissing = 1;
		}
		else if (strcmp(argv[argi], "-printFound") ==0) {
			foundFileName = 1;
			printFound = 1;
		}
		else if (strcmp(argv[argi], "-refSeq") == 0) {
			compareToRef = 1;
			refSeqName = argv[++argi];
		}
		++argi;
	}
	DNASequence refSeq;
	if (compareToRef == 1) {
		SeqReader::GetSeq(refSeqName, refSeq, SeqReader::noConvert);
	}

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
	std::ofstream foundOut, missingOut;
	if (printFound) {
		openck(foundFileName, foundOut, std::ios::out);
	}
	if (printMissing){
		openck(missingFileName, missingOut, std::ios::out);
	}
	ssize_t totalChars = 0;
	ssize_t totalMismatch = 0;
	ssize_t numMisMatch = 0;
  while(SeqReader::GetSeq(seqIn, seq, SeqReader::noConvert)) {
		//UNUSED//    int p;
    if (seqNumber % 10000 == 0) {
      curTime = clock();
      elapsedTime = (double(curTime) - double(startTime)) / CLOCKS_PER_SEC;
			std::cout << seqNumber << " " << elapsedTime << std::endl;
			if (compareToRef) {
				std::cout << totalChars << " " << totalMismatch << " tm: " << double(totalMismatch) / double(totalChars) << std::endl;
			}
    }
		BW::Query(seq, csa, low, high);
		mult = high - low;
		if (printMissing and mult <= 0) {
			seq.PrintlnSeq(missingOut);
		}
		if (printFound and mult >= 1) {
			seq.PrintlnSeq(foundOut);
		}
		
		if (compareToRef && mult <= 0) {
			ssize_t r, p;
			DNASequence seqRC;
			MakeRC(seq, seqRC);
			ssize_t refMatch, refRCMatch;
			refMatch = seq.length;
			refRCMatch = seq.length;
			for (r = 0; r < refSeq.length - seq.length; r++) {
				ssize_t nm, nrcm;
				nm = 0; nrcm =0;
				for (p = 0; p < seq.length; p++) {
					if (seq.seq[p] != refSeq.seq[r + p]) nm++;
					if (seqRC.seq[p] != refSeq.seq[r + p]) nrcm++;
				}
				if (nm < refMatch)
					refMatch = nm;
				if (nrcm < refRCMatch)
					refRCMatch = nrcm;
			}
			if (refMatch < refRCMatch)
				totalMismatch += refMatch;
			else
				totalMismatch += refRCMatch;
			numMisMatch++;
		}
		totalChars += seq.length;
    ++seqNumber;
  }
	if (compareToRef) {
		std::cout << numMisMatch << " / " << seqNumber << " reads had mismatches." << std::endl;
		std::cout << totalChars << " " << totalMismatch << " tm: " << double(totalMismatch) / double(totalChars) << std::endl;
	}
  return 0;

}

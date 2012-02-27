/***************************************************************************
 * Title:          MatchSeq.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "DNASequence.cpp"
#include "SeqReader.h"
#include "SeqUtils.h"
#include "utils.h"

#include <string>
#include <iostream>
#include <fstream>
#include <vector>

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::cout << "usage:" << argv[0] 
	      << " ref query" << std::endl;
    exit(1);
  }
  std::string refSeqName = argv[1];
  std::string qrySeqName = argv[2];

  DNASequence refSeq, refRC;
  SeqReader::MaskRepeats();
  SeqReader::GetSeq(refSeqName, refSeq, SeqReader::noConvert);
	MakeRC(refSeq, refRC);

  std::vector<DNASequence*> qrySequences;
  DNASequence *qrySeq;
  std::ifstream qryIn;
  openck(qrySeqName, qryIn, std::ios::in);
  while(SeqReader::GetSeq(qryIn, qrySeq, SeqReader::noConvert)) {
    qrySequences.push_back(qrySeq);
  }

  std::vector<std::vector<ssize_t > > qryPositions;
  qryPositions.resize(qrySequences.size());
  ssize_t i;
  ssize_t q;
	for (q = 0; q < qrySequences.size(); q++ ) {
		std::cout << q;
		for (i = 0; i < refSeq.length; i++ ) {
      if (i + qrySequences[q]->length < refSeq.length) {
				if (memcmp(&(refSeq.seq[i]), qrySequences[q]->seq, qrySequences[q]->length) == 0) {
					qryPositions[q].push_back(i);
					std::cout << " " << i;
				}
			}
    }
		std::cout << " rc:";
		for (i = 0; i < refRC.length; i++ ) {
      if (i + qrySequences[q]->length < refSeq.length) {
				if (memcmp(&(refRC.seq[i]), qrySequences[q]->seq, qrySequences[q]->length) == 0) {
					qryPositions[q].push_back(i);
					std::cout << " " << refRC.length - i - 1;
				}
			}
    }

		std::cout << std::endl;
  }
  
  ssize_t s;
  for (q = 0; q < qrySequences.size(); q++ ) {
    std::cout << q << "\t";
    for (s = 0; s < qryPositions[q].size(); s++ ) {
      std::cout << qryPositions[q][s] << "\t";
    }
    std::cout << std::endl;
  }

  return 0;
}

/***************************************************************************
 * Title:          BuildProfile.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>

#include "OrthoRepeatReader.h"
#include "SeqReader.h"
#include "DNASequence.h"
#include "utils.h"
#include "alignutils.h"

int main(int argc, char* argv[]) {


  std::string orthoMapFileName;
  std::string refChromFileName;
  std::string consensusFileName;
  if (argc != 4) {
    std::cout << "usage: buildprofile orhtomapfile seqfile consensus " << std::endl;
    std::cout << " seqfile corresponds to the first sequence in orhtomapfile " << std::endl;
    exit(1);
  }

  orthoMapFileName = argv[1];
  refChromFileName = argv[2];
  consensusFileName = argv[3];

  double scoreMat[5][5];
  double **scoreMatPtr;
  char *home = getenv("HOME");
  scoreMatName = std::string(home) + "/projects/mcsrc/align/data/scoremat.txt";

  ReadScoreMatFile(scoreMatName, scoreMat);
  AssignScoreMatPtr(scoreMatPtr, scoreMat);


  DNASequence consensus;
  SeqReader::GetSeq(consensusFileName, consensus);

  std::ifstream orthoMapFile;
  openck(orthoMapFileName, orthoMapFile);
  
  OrthoRepeatReader orthoReader(orthoMapFile);
  
  std::vector<OrthoRepeat*> orthoRepeats;
  OrthoRepeat *newRepeat;
  while (orthoReader.GetNextRepeat(newRepeat)) {
    orthoRepeats.push_back(newRepeat);
    std::cout << "got repeat " << newRepeat->refStart << " "
	      << newRepeat->refEnd << " " << newRepeat->refType
	      << " " << newRepeat->qryStart 
	      << " " << newRepeat->qryEnd  
	      << " " << newRepeat->qryType << std::endl;
  }

  DNASequence refSeq;
  SeqReader::GetSeq(refChromFileName, refSeq);

  ssize_t r;
  DNASequence repeatSeq;
  OrhtoRepeat *orthoRep;
  
  ProfileCount profileCount(5, consensus.length);

  for (r = 0; r < orthoRepeats.size(); r++) {
    orthoRep = orthoRepeats[r];
    repeatSeq.seq = &refSeq.seq[orthoRep->refStart];
    repeatSeq.length = orthoRep->refEnd - orthoRep->refStart + 1;
    StoreProfile(consensus, repeatSeq, profileCount, scoreMatPtr);
  }
  
}

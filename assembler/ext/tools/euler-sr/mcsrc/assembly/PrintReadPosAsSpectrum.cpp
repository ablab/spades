/***************************************************************************
 * Title:          PrintReadPosAsSpectrum.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "SimpleSequence.h"
#include "SeqReader.h"
#include "DeBruijnGraph.h"
#include "IntegralTupleStatic.h"


int main(int argc, char* argv[]) {

  std::string readsFileName, readPosFileName, tuplesOutFileName;
  ssize_t tupleLength;

  if (argc != 5) {
    std::cout << "usage: printReadPosAsSpectrum readsFile readPosFile tupleSize tupleOut" 
	      << std::endl;
    return 1;
  }
  readsFileName = argv[1];
  readPosFileName = argv[2];
  tupleLength = atoi(argv[3]);
  tuplesOutFileName = argv[4];

  std::ifstream readPosIn;
  std::ofstream tuplesOut;
  openck(readPosFileName, readPosIn, std::ios::in);
  openck(tuplesOutFileName, tuplesOut, std::ios::out);

  // get the reads + reverse complements
  SimpleSequenceList reads;
  ReadSimpleSequences(readsFileName, reads);
  AppendReverseComplements(reads);

  ssize_t numReadPos;
  readPosIn >> numReadPos;
  
  std::vector<ReadPos> readPosList;
  readPosList.resize(numReadPos);
  ssize_t i;
  for (i = 0; i < numReadPos; i++) {
    readPosIn >> readPosList[i];
  }
  readPosIn.close();

  tuplesOut << numReadPos << std::endl;
  unsigned char *tuple = new unsigned char[tupleLength+1];
  tuple[tupleLength] = 0;
  for (i = 0; i < numReadPos; i++ ) {
    memcpy(tuple, &(reads[readPosList[i].read].seq[readPosList[i].pos]), tupleLength);
    tuplesOut << tuple << std::endl;
  }
  return 0;
}

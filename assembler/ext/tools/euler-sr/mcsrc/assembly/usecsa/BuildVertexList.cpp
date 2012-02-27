/***************************************************************************
 * Title:          BuildVertexList.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/28/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <string>
#include "DNASequence.h"
#include "SeqReader.h"
#include "bbbwt/BBBWTQuery.h"
#include "../DeBruijnGraph.h"
#include "PagedList.h"

// program that creates the list of vertices given reads.
// the vertex list is given as a list of read positions

typedef PagedList<ReadPos, 1000> ReadPosList;
ssize_t ReadPos::hashLength = 0;
SimpleSequenceList* ReadPos::sequences = NULL;

void PrintUsage() {
  std::cout << "usage: buildVertexList readsFile outputFile [-vertexSize size] " << std::endl;
  exit(0);
}

int main(int argc, char* argv[]) {
  
  std::string readFileName;
  std::string vertexListFileName;
  ssize_t wordLength = 20;
  if (argc < 3) {
    PrintUsage();
    exit(1);
  }
  int argi = 1;
  readFileName       = argv[argi++];
  vertexListFileName = argv[argi++];

  while (argi < argc) { 
    if (strcmp(argv[argi], "-vertexSize") == 0) {
      wordLength = atoi(argv[++argi]);
    }
    ++argi;
  }
  BBBWT wordCSA;

  SimpleSequenceList sequences;
  ReadSimpleSequences(readFileName, sequences);
  std::cout << "read all sequences " << std::endl;
  AppendReverseComplements(sequences);
  std::cout << "and appended complements " << std::endl;
  ssize_t seq;
  ReadPosList readPositions; // TODO: Fix compiler "warning: 'readPositions.PagedList<ReadPos, 1000l>::curPtr' may be used uninitialized in this function"
  ssize_t low, high;
  ReadPos rp;
  DNASequence tmpSeq;
  std::cout << "storing of length " << wordLength << std::endl;
  for (seq = 0; seq < sequences.size(); seq++) {
    ssize_t p;
    ssize_t searchVal;
    if (seq % 100 == 0) {
      std::cout << "processed " << seq << " reads " << std::endl;
    }
    for (p = 0; p < sequences[seq].length - wordLength + 1; p++ ) {
      if (p % 100 == 0)
				std::cout << "querying " << p << std::endl;
      searchVal = BW::Query(sequences[seq], p, wordLength, wordCSA, low, high);
      if (searchVal <= 0) {
				if (p % 100 == 0)
					std::cout << "storing " << p << std::endl;
				/*
					tmpSeq.seq = &(sequences[seq].seq[p]);
					tmpSeq.length = wordLength;
					tmpSeq.PrintSeq(std::cout);
					std::cout << " " << seq << " " << p;
					std::cout << std::endl;
				*/
				BW::Store( sequences[seq], wordCSA, p, wordLength);
				rp.read = seq;
				rp.pos  = p;
				readPositions.Append(rp);
      }
    }
  }
  std::cout << "stored " << readPositions.size() << " read positions " << std::endl;
  std::ofstream vertexListOut;
  openck(vertexListFileName, vertexListOut, std::ios::out);

  vertexListOut << readPositions;
  vertexListOut.close();
  return 0;
}

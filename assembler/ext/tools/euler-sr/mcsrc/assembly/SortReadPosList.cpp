/***************************************************************************
 * Title:          SortReadPosList.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <vector>
#include "utils.h"
#include "DeBruijnGraph.h"
#include "DNASequence.h"
#include "SeqReader.h"
#include "SeqUtils.h"

void PrintUsage() {
  std::cout << "usage: sortreadpos readPos.in sequences tupleLength readPos.out " << std::endl;
  std::cout << "readPos.in -  the list of read positions.  The first line has the number of " << std::endl
	    << "              number of read positions, and each line contains one read position, " << std::endl
	    << "               (read, pos). " << std::endl
	    << " sequences -  the list of sequences (without their reverse complement " << std::endl
	    << " tupleLength- the length of the tuple to sort on " << std::endl
	    << "readPos.out-  Sort, and store in this file. " << std::endl;
}

int main(int argc, char* argv[]) {

  if (argc < 5) {
    PrintUsage();
    exit(1);
  }
  std::string inFileName, outFileName, readFileName;
  int tupleSize;
  inFileName = argv[1];
  readFileName = argv[2];
  tupleSize = atoi(argv[3]);
  outFileName = argv[4];
  ReadPositions readPositions;
  std::ifstream in;
  
  openck(inFileName, in, std::ios::in);

  // Read in the read positions
  in >> readPositions;

  // Read in all reads
  SimpleSequenceList sequences;
  ReadSimpleSequences(readFileName, sequences);
  AppendReverseComplements(sequences);

  // Set up the functor to compare the tuples
  CompareTuples comp;
  comp.sequencesPtr = &sequences;
  comp.sequencesPtr = &sequences;
  comp.length = tupleSize;

  // Do quicksort on the read positions
  std::sort(readPositions.begin(), readPositions.end(), comp);

  // Done, print them
  std::ofstream out;
  openck(outFileName, out, std::ios::out);
  out << readPositions;
  out.close();

  return 0;
}

/***************************************************************************
 * Title:          CountSequenceMatches.cpp 
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
    std::cout << "usage: countSequenceMatches reads csa [-verbose] [-maxLength l]" << std::endl;
    std::cout << "	 count the number of times a sequence is in a csa, reporting the number " << std::endl
	      << "       sequences that have at least one match. " << std::endl;
    return 1;
  }
  ssize_t verbose = 0;
  seqFileName = argv[1];
  csaFileName = argv[2];
  int argi = 3;
  while (argi < argc) {
    if (strcmp(argv[argi], "-verbose") == 0) {
      verbose = 1;
    }
		else if (strcmp(argv[argi], "-maxLength") == 0) {
			maxQueryLength = atoi(argv[++argi]);
		}
    ++argi;
  }

  BW::Read(csaFileName, csa);
  ssize_t low, high;
  ssize_t numMatched = 0;
  NamedSequenceList sequences;
  ReadNamedSequences(seqFileName, sequences);
  ssize_t i, mult;
  for (i = 0; i < sequences.size();i++) {
		if (maxQueryLength != -1 and maxQueryLength < sequences[i].length)
			continue;
    BW::Query(sequences[i], csa, low, high);
    mult = high - low;
    if (mult > 0) {
      numMatched++;
      if (verbose) {
				std::cout << mult << " ";
      }
    }
    else {
      if (verbose) {
				std::cout <<"0 ";
      }
    }
    if (verbose) {
      std::cout << sequences[i].namestr << " " 
								<< sequences[i].length << std::endl;
    }
  }
  
  std::cout << numMatched << " " << sequences.size() << std::endl;
  return 0;
}

/***************************************************************************
 * Title:          TruncateHomopolymericRegions.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  10/24/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "DNASequence.h"
#include "SeqReader.h"
#include "utils.h"
#include <iostream>


void PrintUsage() {
  std::cout << "usage: truncateHPR readsIn readsOut [-length l] " << std::endl;
  std::cout << "    truncate homopolymeric regions to length 'l' " << std::endl;
}

int main(int argc, char* argv[]) {

  std::string readsInName, readsOutName;
  ssize_t minLength = 4;
  
  if (argc < 3) {
    PrintUsage();
    exit(1);
  }
  int argi = 0;
  readsInName = argv[++argi];
  readsOutName = argv[++argi];	
  std::cout << readsInName << " " << readsOutName << std::endl;
  argi++;
  while (argi< argc) { 
    if (strcmp(argv[argi], "-length") == 0) {
      minLength = atoi(argv[++argi]);
    }
    else {
      PrintUsage();
      return 1;
    }
    ++argi;
  }
  DNASequenceList reads;
  ReadDNASequences(readsInName, reads);

  ssize_t r;
  ssize_t totalDeleted = 0;
  ssize_t numReadsModified = 0;
  for (r = 0; r < reads.size(); r++ ) {

    ssize_t p;
    ssize_t hpStart;
    ssize_t nDeleted = 0;
    hpStart = 0;
    ssize_t readModified = 0;
    for (p = 1; p < reads[r].length; p++) { 
      /*
	std::cout << "old seq ";
	reads[r].PrintSeq(std::cout);
	std::cout << std::endl;
      */
      if (reads[r].seq[hpStart] != reads[r].seq[p]) {
	if (p - hpStart > minLength) {
	  readModified = 1;
	  nDeleted = (p - hpStart - minLength);
	  ssize_t p2;
	  for (p2 = hpStart + minLength; p2 < reads[r].length-nDeleted; p2++) {
	    reads[r].seq[p2] = reads[r].seq[p2+nDeleted];
	  }
	  reads[r].length -= nDeleted;
	}
	hpStart = p;
      }
    }
    /*
      if (readModified) {
      std::cout << "new seq " << std::endl;
      reads[r].PrintSeq(std::cout);
      std::cout << std::endl;
      }
    */

    numReadsModified += readModified;
    totalDeleted += nDeleted;
  }
  
  std::ofstream readsOut;
  openck(readsOutName, readsOut, std::ios::out);
  for (r = 0; r < reads.size(); r++ ) {
    reads[r].PrintSeq(readsOut);
    readsOut << std::endl;
  }
  readsOut.close();
  
  std::cout << "modified " << numReadsModified << " deleted " << totalDeleted << std::endl;
  return 0;
}
  

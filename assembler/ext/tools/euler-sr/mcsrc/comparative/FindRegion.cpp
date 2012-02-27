/***************************************************************************
 * Title:          FindRegion.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "lav/LAVUtils.h"
#include "lav/LAVFile.h"
#include "lav/LAVReader.h"

#include "DNASequence.h"
#include "SeqReader.h"
#include "utils.h"
#include "FragmentUtils.h"


#include <sys/types.h>
#include <unistd.h>
#include <stdlib.h>
#include <sstream>

typedef ssize_t Pair[2];


int main(int argc, char* argv[]) {
  //
  // locatefragments
  // 
  // The purpose of this program is to take a list of fragments 
  // of one sequence, a pairwise alignment of that sequence and another
  // and to try and find them in the second sequence.
  // 
  std::string fragmentFileName;
  std::string positionFileName;
  std::string alignmentFileName;
  std::string subjectFileName;
  std::string blastdb;

  // configure defaults
  blastdb    = "";
  alignmentFileName = "";
  // process command line options
  if (argc < 4) {
    std::cout << "usage: locator fragmentFile positionFile alignmentFile subjectFileName " << std::endl;
    exit(1);
  }
  int argi = 1;
  alignmentFileName= argv[argi++];

  std::vector<ssize_t> positions;
  while (argi < argc) {
    positions.push_back(atoi(argv[argi++]));
  }

  // Read in the alignment that specifies the mapping of 
  // the reference sequence to the query
  LAVFile alignment;
  LAVReader::ReadLAVFile(alignmentFileName, alignment);

  std::vector<LAVBlock*> refOrderBlocks;
  BlockReferenceOrder refOrder;
  StoreBlockArray(alignment, refOrderBlocks, refOrder);

  //
  // Sanity check.
  //


  ssize_t pos;
  // attempt to locate each fragment in the query sequence
  for (pos = 0; pos < positions.size(); pos++) {
    // Find the positions surrounding the query
    ssize_t sbjctBeginPos, sbjctEndPos;
    if (FindSurroundingRegion(positions[pos], positions[pos],
			      refOrderBlocks,
			      sbjctBeginPos, sbjctEndPos,
			      -1)) {
      std::cout << positions[pos] << " is between "
		<< sbjctBeginPos << " and " 
		<< sbjctEndPos << std::endl;
    }
  } 
}

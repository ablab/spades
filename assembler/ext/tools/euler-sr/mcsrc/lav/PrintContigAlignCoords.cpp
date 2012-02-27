/***************************************************************************
 * Title:          PrintContigAlignCoords.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <string>
#include <unistd.h>


#include "utils.h"
#include "blocks/BlockDB.h"
#include "lav/LAVReader.h"
#include "lav/LAVFile.h"
#include "lav/LAVUtils.h"
#include "lav/LAVPrinter.h"
#include "net/NetDB.h"


int main(int argc, char *argv[]) {

  std::string lavFileName, refSpecies, qrySpecies, seqName, database;
  std::string outFileName;
  std::string seqTable;
  ssize_t qrySeqId, refSeqId;
  refSpecies = "";
  qrySpecies = "";
  seqName    = "";
  database   = "alignments";

  if (argc < 3) { 
		std::cout << "printContigAlignCoords" << std::endl;
		std::cout << "Given an alignment file of a reference sequence" << std::endl;
		std::cout << "and many short subsequences, print the coordinates. " << std::endl;
		std::cout << "of the highest scoring block. " << std::endl;
		std::cout << "This assumes that there are no large gaps in the subsequenes." 
							<< std::endl << std::endl;
    std::cout << "Usage: printContigAlignCoords lavFile outfile " << std::endl;
		std::cout << std::endl;


    exit(1);
  }

  lavFileName = argv[1];
  outFileName = argv[2];

	LAVFile lavFile;
  LAVReader::ReadLAVFile(lavFileName, lavFile);

	ssize_t a;
	LAVAlignedContig *alignedContig;
	for (a = 0; a < lavFile.alignments.size(); a++ ){ 
		alignedContig = lavFile.alignments[a];
		ssize_t qryStart = alignedContig->qryContig.start;
		ssize_t qryEnd   = alignedContig->qryContig.end;
		std::cout << alignedContig->qryContig.sequenceName << " " << qryStart << " " << qryEnd 
							<< " with: " << alignedContig->alignments.size() << std::endl;
	}

	return 0;
}

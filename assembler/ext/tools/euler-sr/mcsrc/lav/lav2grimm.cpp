/***************************************************************************
 * Title:          lav2grimm.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
// std library
#include <iostream>
#include <fstream>
#include <string>
#include <exception>
#include <stdio.h>
#include <sstream>

// mine
#include "utils.h"
#include "lav/LAVTable.h"
#include "lav/LAVFile.h"
#include "lav/LAVAlignedContig.h"
#include "lav/LAVBlock.h"
#include "lav/LAVReader.h"
#include "lav/LAVUtils.h"

int main(int argc, char* argv[]) {

  if (argc < 3) {
    std::cout << "usage: lav2grimm inlav outpoints " << std::endl;
    exit(0);
  }
  std::string inLavName = argv[1];
  std::string outBlockName = argv[2];
  std::ofstream blockOut;
  openck(outBlockName, blockOut, std::ios::out);
  LAVFile lavFile;
  LAVReader::ReadLAVFile(inLavName, lavFile);

  //UNUSED+// ssize_t bl;
  ssize_t a, b,  index;
  LAVAlignedContig *alignedContig;
  LAVBlock *block;
  index = 0;
  for (a = 0; a < lavFile.size(); a++) {
    alignedContig = lavFile.alignments[a];
    for (b = 0; b < alignedContig->size(); b++) {
      block = alignedContig->alignments[b];
      blockOut << index << " 1 " 
	       << block->refBegin << " " << block->refEnd - block->refBegin
	       << " +1 "  // always +1 on reference strand
	       << " 1 "; // chromosome
      if (alignedContig->qryContig.strand == 0) {
	blockOut << block->qryBegin 
		 << " " << block->qryEnd - block->qryBegin;
	blockOut << " +1 ";
      }
      else {
	blockOut << (alignedContig->qryContig.end - 
		     alignedContig->qryContig.start + 1) - block->qryEnd 
		 << " " << block->qryEnd - block->qryBegin ;
	blockOut << " -1 ";
      }
      

      blockOut << std::endl;
    ++index;
    }
  }
  blockOut.close();

}



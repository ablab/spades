/***************************************************************************
 * Title:          LavToBlocks.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <unistd.h>
#include "lav/LAVFile.h"
#include <iostream>
#include <map>
#include <vector>
#include "lav/LAVReader.h";
#include "blocks/BlockGraph.h"
#include "blocks/Block.h"


void InitEnv(int argc, char *argv[], 
	     std::string &lavFile,
	     std::string &blocksFile);

void PrintUsage();
  
int main(int argc, char* argv[]) {
  
  std::string lavFileName, blocksFileName;
  ssize_t a, b, c;
  BlockGraph graph;
  LAVFile lavFile;
  ssize_t seqIndex;
  ssize_t refIndex, qryIndex;
  seqIndex = 0;
  InitEnv(argc, argv, lavFileName, blocksFileName);
  
  LAVReader::ReadLAVFile(lavFileName, lavFile);
  LAVAlignedContig *alignedContig;
  LAVBlock *lavBlock;
  // Create the vertices that represent the sequences
  for (a = 0; a < lavFile.alignments.size(); a++) {
    alignedContig = lavFile.alignments[a];
    if (alignedContig->refContig.strand == 0 && 
	alignedContig->qryContig.strand == 0) {
      refIndex = graph.AddVertex(alignedContig->refContig.start, 
				 alignedContig->refContig.end);
      qryIndex = graph.AddVertex(alignedContig->qryContig.start, 
				 alignedContig->qryContig.end);

      // Process the alignments that are stored for these contigs
      for (b = 0; b < alignedContig->alignments.size(); b++) {
	// Each alignment is a collection of contiguous blocks
	// join them
	lavBlock = alignedContig->alignments[b];
	for (c = 0; c < lavBlock->size(); c++) {
	  graph.GlueSequences(refIndex, lavBlock->refALBegin[c], 
			      lavBlock->refALEnd[c],
			      qryIndex, lavBlock->qryALBegin[c], 
			      lavBlock->qryALEnd[c]);
	}
	ssize_t m;
      }
    }
  }
}


void InitEnv(int argc, char *argv[], 
	     std::string &lavFile,
	     std::string &blocksFile) {
  ssize_t copt;
  lavFile = "";
  blocksFile = "";
  // no options for now
  while ( (copt=getopt(argc, argv, "")) != EOF) {
    
  }
  ssize_t ind = optind;

  if (ind >= argc) {
    PrintUsage();
    exit(0);
  }
  lavFile = argv[ind];
  ind++;
  if (ind >= argc) {
    PrintUsage();
    exit(0);
  }
  blocksFile = argv[ind];
}
void PrintUsage()
 {
  std::cout << "lavtoblocks: a file to take a number of lav pairwise alignments " << std::endl
	    << " and create a multiple alignment from the blocks.  For now  " << std::endl
	    << " this considers the first alignment in the lav files to be the reference " << std::endl
	    << " alignment, and all others are aligned according to that.  It will take " << std::endl
	    << " some thought for how to have this in the general case. " << std::endl;
  std::cout << "usage: " << std::endl;
  std::cout << "   lavToBlocks lavFile blocksFile" << std::endl;
}

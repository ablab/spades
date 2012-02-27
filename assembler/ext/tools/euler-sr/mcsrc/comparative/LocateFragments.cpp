/***************************************************************************
 * Title:          LocateFragments.cpp 
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
  fragmentFileName = argv[argi++];
  positionFileName = argv[argi++];
  alignmentFileName= argv[argi++];
  subjectFileName  = argv[argi++];
  
  // Read in fragments to search for
  std::ifstream fragIn;
  openck(fragmentFileName, fragIn);
  std::vector<DNASequence*> fragments;
  DNASequence *fragPtr;
  while(SeqReader::GetSeq(fragIn, fragPtr, SeqReader::noConvert)) {
    fragments.push_back(fragPtr);
  }
  

  // Read in the locations of those fragments
  std::vector<ssize_t*> fragPositions;
  std::ifstream posIn;
  openck(positionFileName, posIn, std::ios::in);
  ssize_t p = 0;
  while (posIn) {
    fragPositions.push_back(new ssize_t[2]);
    posIn >> fragPositions[p][0] >> fragPositions[p][1];
    p++;
  }


  // Read in the alignment that specifies the mapping of 
  // the reference sequence to the query
  LAVFile alignment;
  LAVReader::ReadLAVFile(alignmentFileName, alignment);

  std::vector<LAVBlock*> refOrderBlocks;
  BlockReferenceOrder refOrder;
  StoreBlockArray(alignment, refOrderBlocks, refOrder);

  // Read in the subject where the fragments will be searched.
  DNASequence sequence;
  SeqReader::GetSeq(subjectFileName, sequence);

  //
  // Sanity check.
  //

  // Prepare the name of the coordinates of frames
  std::string queryName, sbjctName;
  MakeTempName(queryName, ".locator.qry.fasta");
  MakeTempName(sbjctName, ".locator.ref.fasta");

  ssize_t frag;
  // attempt to locate each fragment in the query sequence
  for (frag = 0; frag < fragments.size(); frag++ ) {
    // Find the positions surrounding the query
    ssize_t sbjctBeginPos, sbjctEndPos;
    FindSurroundingRegion(fragPositions[frag][0],
			  fragPositions[frag][1],
			  refOrderBlocks,
			  sbjctBeginPos, sbjctEndPos,
			  sequence.length);

    std::cout << "searching subject from " << sbjctBeginPos << " to " <<  sbjctEndPos << std::endl;
    
    assert(sbjctEndPos > sbjctBeginPos);

    // Output the subject sequence to a file
    PrintTempSeq(sequence, sbjctBeginPos, sbjctEndPos, sbjctName);
    PrintTempSeq(*fragments[frag], 0, fragments[frag]->length-1, queryName);

    Pos position;
    // Locate query in subjec using bl2seq

    ssize_t blockIndex;
    ssize_t beginBlock, endBlock;
    ssize_t qryStart, qryEnd;

    if (RunBlast(queryName, sbjctName, position) ) {
      std::cout << "blasting fragment resulted in: "
		<< position.rBegin << " " << position.rEnd << " " 
		<< position.qBegin << " " << position.qEnd << std::endl;
    }

    // Clean up
    system(std::string("rm " + sbjctName).c_str());
    system(std::string("rm " + queryName).c_str());
  }
  // Look for alignments containing each position
}





	
      /*
	if (FindBlock(positions[p].qBegin, blocks, beginBlock, contained)) {
	if (contained) {
	qryStart = FindPosition(positions[p].qBegin, *blocks[beginBlock]);
	std::cout << "position: " << positions[p].qBegin
	<< " is contained in block: " << beginBlock
	<< " " << blocks[beginBlock]->refBegin 
	<< " " << blocks[beginBlock]->refEnd << " "
	<< " " << blocks[beginBlock]->qryBegin 
	<< " " << blocks[beginBlock]->qryEnd << " at: "
	<< qryStart << std::endl;
		    
	}
	else {
	std::cout << "position: " << positions[p].qBegin
	<< " is after block: " << beginBlock
	<< " " << blocks[beginBlock]->refEnd
	<< " " << blocks[beginBlock+1]->refBegin << " "
	<< " " << blocks[beginBlock]->qryEnd
	<< " " << blocks[beginBlock+1]->qryBegin << std::endl;
	}
	}
	if (FindBlock(positions[p].qEnd, blocks, beginBlock, contained)) {
	if (contained) {
	qryEnd = FindPosition(positions[p].qEnd, *blocks[beginBlock]);
	std::cout << "position: " << positions[p].qEnd
	<< " is contained in block: " << beginBlock
	<< " " << blocks[beginBlock]->refBegin 
	<< " " << blocks[beginBlock]->refEnd << " "
	<< " " << blocks[beginBlock]->qryBegin 
	<< " " << blocks[beginBlock]->qryEnd 
	<< " at: " << qryEnd << std::endl;
	}
	else {
	std::cout << "position: " << positions[p].qEnd
	<< " is after block: " << beginBlock
	<< " " << blocks[beginBlock]->refEnd
	<< " " << blocks[beginBlock+1]->refBegin << " "
	<< " " << blocks[beginBlock]->qryBegin 
	<< " " << blocks[beginBlock+1]->qryEnd << std::endl;
	}
	}
      */

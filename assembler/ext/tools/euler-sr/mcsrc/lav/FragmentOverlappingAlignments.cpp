/***************************************************************************
 * Title:          FragmentOverlappingAlignments.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <set>

#include "utils.h"
#include "lav/LAVTable.h"
#include "lav/LAVFile.h"
#include "lav/LAVAlignedContig.h"
#include "lav/LAVBlock.h"
#include "lav/LAVReader.h"
#include "lav/LAVUtils.h"
#include "lav/LAVPrinter.h"

#include "align/alignutils.h"
#include "DNASequence.h"
#include "SeqReader.h"

const ssize_t REFERENCE = 0;
const ssize_t QUERY     = 1;

ssize_t SplitBlock(ssize_t sequence, ssize_t position, 
	       LAVBlock *block, 
	       LAVBlock *&startBlock, LAVBlock *&endBlock,
	       ssize_t strand);

ssize_t SplitBlockRev(ssize_t sequence, ssize_t position, 
		  LAVBlock *block, 
		  LAVBlock *&startBlock, LAVBlock *&endBlock,
		  ssize_t strand);

void StoreCoordinates(LAVFile &lavFile, std::set<ssize_t> &coords, ssize_t sequence);

void SplitBlocks(LAVFile &lavFile, std::set<ssize_t> coords, ssize_t sequence);

void AdvanceIt(std::set<ssize_t>::iterator &setIt, ssize_t strand) {
  if (strand == 0)
    ++setIt;
  else
    --setIt;
}
    

int main(int argc, char* argv[]) {

  if (argc < 3) {
    std::cout << "usage: " << argv[0] << " inLavFile outlavFile " << std::endl;
    exit(0);
  }
  std::string inLavName  = argv[1];
  std::string outLavName = argv[2];
  LAVFile lavFile;
  DNASequence refSeq, qrySeq;
  LAVReader::ReadLAVFile(inLavName, lavFile);
  std::map<char, double> keywordOptions;

  //UNUSED// ssize_t refBlockStart, refBlockEnd;
  //UNUSED// ssize_t qryBlockStart, qryBlockEnd;
  std::set<ssize_t> refCoords, qryCoords;
  //UNUSED// ssize_t  a, b;
  //UNUSED// LAVAlignedContig *alignedContig;
  //UNUSED// LAVBlock *block;
  //UNUSED// ssize_t refLength, qryLength;

  // change UCSC 1-based coordinate system.
  BeginAtZero(lavFile);
  StoreCoordinates(lavFile, refCoords, REFERENCE);
  
  // Now split all other blocks into smaller blocks 
  // according to the reference sequence
  
  SplitBlocks(lavFile, refCoords, REFERENCE);

  // Find the coordinates of where the query intersects
  // with the reference.
  StoreCoordinates(lavFile, qryCoords, QUERY);

  // Split along the query
  SplitBlocks(lavFile, qryCoords, QUERY);

  refCoords.clear();
  qryCoords.clear();

  // One more time. Add comment later explaining why... inversions
  // may cause the reference to intersect the query that intersects
  // the reference.  
  StoreCoordinates(lavFile, refCoords, REFERENCE);
  SplitBlocks(lavFile, refCoords, REFERENCE);

  // And yet again for the second splitting of the query.
  StoreCoordinates(lavFile, qryCoords, QUERY);
  SplitBlocks(lavFile, qryCoords, QUERY);

  refCoords.clear();
  qryCoords.clear();
  // One more time. Add comment later explaining why... inversions
  // may cause the reference to intersect the query that intersects
  // the reference.  
  StoreCoordinates(lavFile, refCoords, REFERENCE);
  SplitBlocks(lavFile, refCoords, REFERENCE);

  // And yet again for the second splitting of the query.
  StoreCoordinates(lavFile, qryCoords, QUERY);
  SplitBlocks(lavFile, qryCoords, QUERY);

  BeginAtOne(lavFile);
  std::ofstream outLAVFile;
  openck(outLavName, outLAVFile, std::ios::out);
  LAVPrinter::PrintLAVFile(lavFile, outLAVFile);
}

ssize_t SplitBlockRev(ssize_t sequence, ssize_t position, 
		  LAVBlock *block, 
		  LAVBlock *&startBlock, LAVBlock *&endBlock,
		  ssize_t strand) {
  //UNUSED+// ssize_t splitInQryInterval;
  ssize_t splitInInterval ;
  splitInInterval = 0;

  ssize_t alnBegin;
  ssize_t alnEnd;
  ssize_t alnNextBegin;
  ssize_t i;
  ssize_t splitIndex;
  
  splitIndex = -1;
  for (i = 0; i < block->size()-1; i++) {
    alnBegin = block->IntvBegin(i, sequence);
    alnEnd   = block->IntvEnd(i, sequence);
    alnNextBegin = block->IntvBegin(i+1, sequence);;

    if (position >= alnBegin and
	position < alnNextBegin) {
      splitIndex = i;
      if (position < alnEnd) {
	splitInInterval = 1;
      }
      break;
    }
  }
  if (splitIndex == -1) {
    // check the last block (or only one if there is just one)
    assert(position >= block->IntvBegin(i,sequence) and
	   position < block->IntvEnd(i,sequence));
    if ( position >= block->IntvBegin(i,sequence) and
	 position < block->IntvEnd(i,sequence) ) {
      splitIndex         = i;
      splitInInterval    = 1;
    }
  }
    
  startBlock = new LAVBlock;
  endBlock   = new LAVBlock;
  startBlock->strand = block->strand;
  endBlock->strand = block->strand;
  startBlock->refBegin = block->refBegin;
  startBlock->qryBegin = block->qryBegin;
    
  if (splitIndex >= 0) {
    // split the block into two halves
    startBlock->Resize(splitIndex + 1);
    startBlock->CopyIntervals(0, *block, 0, splitIndex + 1 - splitInInterval);
    
    if (splitInInterval) {
      if (sequence == REFERENCE) {
	startBlock->refEnd = position;
	startBlock->refALEnd[splitIndex] = startBlock->refEnd;
	startBlock->qryEnd = block->qryALBegin[splitIndex] + 
	  (position - block->refALBegin[splitIndex]);
	startBlock->qryALEnd[splitIndex] = startBlock->qryEnd;
      }
      else {
	startBlock->refEnd = block->refALBegin[splitIndex] + 
	  (position - block->qryALBegin[splitIndex]);
	startBlock->refALEnd[splitIndex] = startBlock->refEnd;
	startBlock->qryEnd = position;
	startBlock->qryALEnd[splitIndex] = startBlock->qryEnd;
      }
    }
    else {
      startBlock->refEnd = block->refALEnd[splitIndex];
      startBlock->refALEnd[splitIndex] = startBlock->refEnd;
      startBlock->qryEnd = block->qryALEnd[splitIndex];
      startBlock->qryALEnd[splitIndex] = startBlock->qryEnd;
    }      
    startBlock->refALBegin[splitIndex] = block->refALBegin[splitIndex];
    startBlock->qryALBegin[splitIndex] = block->qryALBegin[splitIndex];
    startBlock->alIdentity[splitIndex] = block->alIdentity[splitIndex];
  }
  else {
    // ended in a gap
    startBlock->refEnd = block->refALEnd[splitIndex];
    startBlock->qryEnd = block->qryALEnd[splitIndex];
  }

  // the second half
  if (splitIndex <= block->size() -1 ) {
    // Allocate for the remaining intervals + one 
    // if the split happened inside an interval.
    endBlock->Resize(block->size() - splitIndex - 1 + splitInInterval);
    // If the split is in an interval, create a new one with 
    // the second half of the interval
    if (splitInInterval) {
      if (sequence == REFERENCE) {
	endBlock->refBegin      = position + 1;
	endBlock->refALBegin[0] = endBlock->refBegin;
	
	endBlock->qryBegin      = block->qryALBegin[splitIndex] + 
	  (position - block->refALBegin[splitIndex]) + 1;
	endBlock->qryALBegin[0] = endBlock->qryBegin;
      }
      else {
	endBlock->refBegin      = block->refALBegin[splitIndex] + 
	  (position - block->qryALBegin[splitIndex]) + 1;
	endBlock->refALBegin[0] = endBlock->refBegin;
	
	endBlock->qryBegin      = position + 1;
	endBlock->qryALBegin[0] = endBlock->qryBegin;
      }
      // Common to both cases
      endBlock->refALEnd[0]   = block->refALEnd[splitIndex];
      endBlock->qryALEnd[0]   = block->qryALEnd[splitIndex];
      endBlock->alIdentity[0] = block->alIdentity[splitIndex];
    }
    else {
      assert(block->refALBegin.size() > 1);
      endBlock->refBegin = block->refALBegin[splitIndex+1];
      endBlock->qryBegin = block->qryALBegin[splitIndex+1];
    }

    endBlock->CopyIntervals(splitInInterval, *block, splitIndex+1, block->size());

    endBlock->refEnd = block->refEnd;
    endBlock->qryEnd = block->qryEnd;
  }
   
  assert(endBlock->size() > 0);
  assert(startBlock->size() > 0);
  /*
  std::cout << "from: (" <<  block->refBegin << ", " << block->refEnd
	    << ") (" << block->qryBegin << ", " << block->qryEnd << ") " << std::endl;
  std::cout << "split reverse strand on: " << position << " " << sequence << " " 
	    << strand << std::endl;         
  std::cout << "two brand new intervals: (" << startBlock->refBegin << ", " << startBlock->refEnd
	    << ") (" << startBlock->qryBegin << ", " << startBlock->qryEnd << ") " << std::endl;
  std::cout << "                       : (" << endBlock->refBegin << ", " << endBlock->refEnd
	    << ") (" <<  endBlock->qryBegin << ", " << endBlock->qryEnd << ") " << std::endl;
  */
  // All sorts of debugging macros that catch various bugs I've had in the past
  
  assert(startBlock->refEnd >= startBlock->refBegin);
  assert(startBlock->qryEnd >= startBlock->qryBegin);
  assert(endBlock->refEnd >= endBlock->refBegin);
  assert(endBlock->qryEnd >= endBlock->qryBegin);
  assert(endBlock->refEnd == endBlock->refALEnd[endBlock->size() - 1]);
  assert(endBlock->qryEnd == endBlock->qryALEnd[endBlock->size() - 1]);
  
  assert(startBlock->refALBegin[0] <= startBlock->refALEnd[0]);
  assert(startBlock->qryALBegin[0] <= startBlock->qryALEnd[0]);
  assert(endBlock->refALBegin[0] <= endBlock->refALEnd[0]);
  assert(endBlock->qryALBegin[0] <= endBlock->qryALEnd[0]);

  assert(startBlock->refALBegin[startBlock->size() - 1] <= 
	 startBlock->refALEnd[startBlock->size() - 1]);
  assert(startBlock->qryALBegin[startBlock->size() - 1] <= 
	 startBlock->qryALEnd[startBlock->size() - 1]);

  assert(endBlock->refALBegin[endBlock->size() - 1] <= 
	 endBlock->refALEnd[endBlock->size() - 1]);
  assert(endBlock->qryALBegin[endBlock->size() - 1] <= 
	 endBlock->qryALEnd[endBlock->size() - 1]);
  return 1;
}


ssize_t SplitBlock(ssize_t sequence, ssize_t position, 
	       LAVBlock *block, 
	       LAVBlock *&startBlock, LAVBlock *&endBlock,
	       ssize_t strand) {

  ssize_t i;
  ssize_t splitIndex = -1;
  ssize_t splitInRefInterval, splitInQryInterval, splitInInterval;
  splitInRefInterval = 0;
  splitInQryInterval = 0;
  splitInInterval = 0;

  if (sequence == REFERENCE) {
    // sanity check on splitting forward strand
    if (strand == 0)
      assert(position >= block->refBegin and 
	     position <= block->refEnd);
    /*
    std::cout << "splitting block " << block->refBegin << " ... " 
	      << block->refEnd << " at " << position << std::endl;
    */
    for (i = 0; i < block->size()-1; i++) {
      if (position > block->refALBegin[i] and 
	  position <= block->refALBegin[i+1]) {
	splitIndex = i;
	if (position <= block->refALEnd[i]) {
	  splitInRefInterval = 1;
	  splitInInterval = 1;
	}
	break;
      }
    }
    if (splitIndex == -1) {
      // check the last block (or only one if there is just one)
      assert(position > block->refALBegin[i] and
	     position <= block->refALEnd[i]);
      if (position >= block->refALBegin[i] and 
	  position <= block->refALEnd[i]) {
	splitIndex         = i;
	splitInInterval    = 1;
	splitInRefInterval = 1;
      }
    }
    if (splitIndex == -1) {
      assert(0);
      return 0;
    }
  }

  if (sequence == QUERY) {
    /*
    std::cout << "splitting block on query " << block->qryBegin << " ... " 
	      << block->qryEnd << " at " << position << std::endl;
    */
    ssize_t qryAlnBegin;
    ssize_t qryAlnEnd;
    ssize_t qryAlnNextBegin;
      
    for (i = 0; i < block->size()-1; i++) {
      if (strand == 0){
	qryAlnBegin = block->qryALBegin[i];
	qryAlnEnd   = block->qryALEnd[i];
	qryAlnNextBegin = block->qryALBegin[i+1];
      }
      else {
	qryAlnBegin = block->qryALBegin[i];
	qryAlnEnd   = block->qryALEnd[i];
	qryAlnNextBegin = block->qryALBegin[i+1];
      }

      if (position > qryAlnBegin and
	  position <= qryAlnNextBegin) {
	splitIndex = i;
	if (position <= qryAlnEnd) {
	  splitInQryInterval = 1;
	  splitInInterval = 1;
	}
	break;
      }
    }
    if (splitIndex == -1) {
      // check the last block (or only one if there is just one)
      assert(position > block->qryALBegin[i] and
	     position <= block->qryALEnd[i]);
      if ((strand == 0 and 
	   position >= block->qryALBegin[i] and
	   position <= block->qryALEnd[i])) {
	splitIndex         = i;
	splitInInterval    = 1;
	splitInQryInterval = 1;
      }
      if ((strand == 1 and 
	   position >= block->qryALBegin[i] and
	   position <= block->qryALEnd[i])) {
	splitIndex         = i;
	splitInInterval    = 1;
	splitInQryInterval = 1;
      }
    }
    if (splitIndex == -1) {
      std::cout << "bad query split " << std::endl;
      assert(0);
      return 0;
    }
  }

  startBlock = new LAVBlock;
  endBlock   = new LAVBlock;

  startBlock->strand = block->strand;
  endBlock->strand = block->strand;
  startBlock->refBegin = block->refBegin;
  startBlock->qryBegin = block->qryBegin;

  if (splitIndex >= 0) {
    // split the block into two halves
    startBlock->refALBegin.resize(splitIndex + 1);
    startBlock->refALEnd.resize(splitIndex + 1);
    startBlock->qryALBegin.resize(splitIndex + 1);
    startBlock->qryALEnd.resize(splitIndex + 1);
    startBlock->alIdentity.resize(splitIndex + 1);
    for (i = 0; i < splitIndex + (1 - splitInInterval); i++) {
      startBlock->refALBegin[i] = block->refALBegin[i];
      startBlock->refALEnd[i]   = block->refALEnd[i];
      startBlock->qryALBegin[i] = block->qryALBegin[i];
      startBlock->qryALEnd[i]   = block->qryALEnd[i];
      startBlock->alIdentity[i] = block->alIdentity[i];
    }
    
    if (splitInInterval) {
      if (splitInRefInterval) {
	/*
      std::cout << "split inside block, start's end was at " 
		<< block->refALEnd[splitIndex]
		<< " now at " << position << std::endl;
	*/
	startBlock->refEnd = position-1;
	startBlock->refALEnd[splitIndex] = position-1;
	startBlock->qryEnd = block->qryALBegin[splitIndex] + 
	  (position - block->refALBegin[splitIndex]) - 1;
	startBlock->qryALEnd[splitIndex] = block->qryALBegin[splitIndex] + 
	  (position - block->refALBegin[splitIndex]) - 1;
      }
      else {
	startBlock->refEnd = block->refALBegin[splitIndex] + 
	  (position - block->qryALBegin[splitIndex]) - 1;
	startBlock->refALEnd[splitIndex] = block->refALBegin[splitIndex] + 
	  (position - block->qryALBegin[splitIndex]) - 1;
	startBlock->qryEnd = position - 1;
	startBlock->qryALEnd[splitIndex] = position - 1;
      }
      
      startBlock->refALBegin[splitIndex] = block->refALBegin[splitIndex];
      startBlock->qryALBegin[splitIndex] = block->qryALBegin[splitIndex];
      startBlock->alIdentity[splitIndex] = block->alIdentity[splitIndex];
    }
    else {
      // ended in a gap
      startBlock->refEnd = block->refALEnd[splitIndex];
      startBlock->qryEnd = block->qryALEnd[splitIndex];
    }

    // the second half
 
   if (splitIndex <= block->size() -1 ) {
      // Allocate for the remaining intervals + one 
      // if the split happened inside an interval.
      endBlock->refALBegin.resize(block->size() - splitIndex - 1 + splitInInterval);
      endBlock->qryALBegin.resize(block->size() - splitIndex - 1 + splitInInterval);
      endBlock->refALEnd.resize(block->size()   - splitIndex - 1 + splitInInterval);
      endBlock->qryALEnd.resize(block->size()   - splitIndex - 1 + splitInInterval);
      endBlock->alIdentity.resize(block->size() - splitIndex - 1 + splitInInterval);
      // If the split is in an interval, create a new one with 
      // the second half of the interval
      if (splitInInterval) {
	if (splitInRefInterval) {
	  /*
	  std::cout << "split in ref interval " << std::endl;
	  */
	  endBlock->refBegin = position;
	  endBlock->refALBegin[0] = position;
	  endBlock->qryBegin = block->qryALBegin[splitIndex] + 
	    (position  - block->refALBegin[splitIndex]);
	  endBlock->qryALBegin[0] = block->qryALBegin[splitIndex] + 
	    (position  - block->refALBegin[splitIndex]);
	}
	else {
	  /*
	  std::cout << "split in qry interval " << std::endl;
	  */
	  // split in the query interval
	  endBlock->refBegin = block->refALBegin[splitIndex] + 
	    (position - block->qryALBegin[splitIndex]);
	  endBlock->refALBegin[0] = block->refALBegin[splitIndex] + 
	    (position - block->qryALBegin[splitIndex]);

	  if (strand == 0) {
	    endBlock->qryBegin      = position;
	    endBlock->qryALBegin[0] = position;
	    endBlock->qryALEnd[0]   = block->qryALEnd[splitIndex];
	  }
	  else {
	    endBlock->qryBegin      = position;
	    endBlock->qryALBegin[0] = position;
	    endBlock->qryALEnd[0]   = block->qryALEnd[splitIndex];
	  }
	}
	endBlock->refALEnd[0]   = block->refALEnd[splitIndex];
	endBlock->qryALEnd[0]   = block->qryALEnd[splitIndex];
	endBlock->alIdentity[0] = block->alIdentity[splitIndex];
      }
      else {
	assert(block->refALBegin.size() > 1);
	endBlock->refBegin = block->refALBegin[splitIndex+1];
	endBlock->qryBegin = block->qryALBegin[splitIndex+1];
      }
      
      for (i = splitIndex + 1; i < block->size(); i++ ) { 
	endBlock->refALBegin[i  - 1 - splitIndex + splitInInterval] = block->refALBegin[i];
	endBlock->refALEnd[i    - 1 - splitIndex + splitInInterval] = block->refALEnd[i];
	endBlock->qryALBegin[i  - 1 - splitIndex + splitInInterval] = block->qryALBegin[i];
	endBlock->qryALEnd[i    - 1 - splitIndex + splitInInterval] = block->qryALEnd[i];
	endBlock->alIdentity[i  - 1 - splitIndex + splitInInterval] = block->alIdentity[i];
      }

      endBlock->refEnd = block->refEnd;
      endBlock->qryEnd = block->qryEnd;
      
    }
   
   assert(endBlock->size() > 0);
   assert(startBlock->size() > 0);
   /*
   std::cout << "from: (" <<  block->refBegin << ", " << block->refEnd
	     << ") (" << block->qryBegin << ", " << block->qryEnd << ") " << std::endl;
   std::cout << "split on: " << position << " " << sequence << " " << strand << std::endl;         
   std::cout << "two brand new intervals: (" << startBlock->refBegin << ", " << startBlock->refEnd
	     << ") (" << startBlock->qryBegin << ", " << startBlock->qryEnd << ") " << std::endl;
   std::cout << "                       : (" << endBlock->refBegin << ", " << endBlock->refEnd
	     << ") (" <<  endBlock->qryBegin << ", " << endBlock->qryEnd << ") " << std::endl;
   */
   // tidy up
   assert(startBlock->refEnd >= startBlock->refBegin);
   assert(startBlock->qryEnd >= startBlock->qryBegin);
   assert(endBlock->refEnd >= endBlock->refBegin);
   assert(endBlock->qryEnd >= endBlock->qryBegin);
   assert(endBlock->refEnd == endBlock->refALEnd[endBlock->size() - 1]);
   assert(endBlock->qryEnd == endBlock->qryALEnd[endBlock->size() - 1]);

   assert(startBlock->refALBegin[startBlock->size() - 1] <= 
	  startBlock->refALEnd[startBlock->size() - 1]);
   assert(startBlock->qryALBegin[startBlock->size() - 1] <= 
	  startBlock->qryALEnd[startBlock->size() - 1]);

   assert(endBlock->refALBegin[endBlock->size() - 1] <= 
	  endBlock->refALEnd[endBlock->size() - 1]);
   assert(endBlock->qryALBegin[endBlock->size() - 1] <= 
	  endBlock->qryALEnd[endBlock->size() - 1]);


    return 1;
  }
  assert(0);
  return 0;
}


void StoreCoordinates(LAVFile &lavFile, std::set<ssize_t> &coords, ssize_t sequence) {
  ssize_t a, b;
  LAVAlignedContig *alignedContig;
  LAVBlock *block;
  ssize_t refLength, qryLength;
  //UNUSED// ssize_t refBlockStart, refBlockEnd, qryBlockStart, qryBlockEnd;
  for (a = 0; a < lavFile.size(); a++) {
    alignedContig = lavFile.alignments[a];
    refLength = alignedContig->refContig.end - 
      alignedContig->refContig.start + 1;
    qryLength = alignedContig->qryContig.end - 
      alignedContig->qryContig.start + 1;

    for (b = 0; b < alignedContig->size(); b++) {
      block = alignedContig->alignments[b];
      if (sequence == REFERENCE) {
	// Store the coordinates of the block to chop other blocks later
	coords.insert(block->refBegin);
	coords.insert(block->refEnd + 1);
      }
      else {
	//UNUSED// ssize_t qryBlockStart, qryBlockEnd;
	// Cache the coordinates of this block
	if (alignedContig->qryContig.strand == 0) {
	  coords.insert(block->qryBegin);
	  coords.insert(block->qryEnd+1);
	}
	else {
	  coords.insert(qryLength - block->qryEnd - 1);
	  coords.insert(qryLength - block->qryBegin - 1 + 1 );
	}
      }
    }
  }
}


void SplitBlocks(LAVFile &lavFile, std::set<ssize_t> coords, ssize_t sequence) {
  ssize_t a, b;
  std::set<ssize_t>::iterator lower, upper, splitIterator;
  LAVBlock *startBlock, *endBlock;
  std::vector<LAVBlock*> newBlocks;
  ssize_t nContigs=  lavFile.size();
  ssize_t strand, length;
  ssize_t refLength, qryLength;
  LAVAlignedContig *alignedContig;
  LAVBlock* block;
  ssize_t refBlockStart, refBlockEnd;
  ssize_t qryBlockStart, qryBlockEnd;

  for (a = 0; a < nContigs; a++) {
    alignedContig = lavFile.alignments[a];
    refLength = alignedContig->refContig.end - 
      alignedContig->refContig.start + 1;
    qryLength = alignedContig->qryContig.end - 
      alignedContig->qryContig.start + 1;

    
    for (b = 0; b < alignedContig->size(); b++) {
      block = alignedContig->alignments[b];
     
      if (sequence == REFERENCE) {
	// Cache the coordinates of this block
	refBlockStart = block->refBegin;
	refBlockEnd   = block->refEnd;
	
	// Store the coordinates of the block to chop other blocks later
	strand = alignedContig->refContig.strand;
	//	if (strand == 0) {
	  lower = coords.upper_bound(refBlockStart);
	  upper = coords.lower_bound(refBlockEnd);
	  if (*lower > *upper) 
	    lower = coords.end();
	length = refLength;
      }
      else {
	// Determine the intervals in the query contig that overlap

	// Cache the coordinates of this block
	if (alignedContig->qryContig.strand == 0) {
	  qryBlockStart = block->qryBegin;
	  qryBlockEnd   = block->qryEnd;
	}
	else {
	  qryBlockStart = qryLength - block->qryEnd - 1;
	  qryBlockEnd   = qryLength - block->qryBegin - 1;
	}
	
	// Store the coordinates of the block to chop other blocks later
	if (alignedContig->qryContig.strand == 0) {
	  lower = coords.upper_bound(qryBlockStart);
	  upper = coords.lower_bound(qryBlockEnd);
	  if (*lower > *upper) 
	    lower = coords.end();
	}
	else {
	  // Need to search for bounds in opposite direction
	  lower = coords.lower_bound(qryBlockEnd);
	  --lower;
	  upper = coords.lower_bound(qryBlockStart);
	  if (*upper > *lower)
	    lower = coords.end();
	}
	strand = alignedContig->qryContig.strand;
	length = qryLength;
      }


      //UNUSED// ssize_t nits = 0;
      if (upper != coords.end() and
	  lower != coords.end() and 
	  *lower != *upper) {
	for (splitIterator = lower; 
	     splitIterator != upper; 
	     AdvanceIt(splitIterator, strand)) {

	  if (sequence == REFERENCE) {
	    if (*splitIterator > block->refBegin){
	      startBlock = NULL;
	      endBlock   = NULL;
	      //	      if (block->strand == 0) {
		if (SplitBlock(sequence, *splitIterator, 
			       block, startBlock, endBlock,
			       strand)) {
		  newBlocks.push_back(startBlock);
		  block = endBlock;
		}
		//	      }
	      /*	      else {
		if (SplitBlockRev(sequence, *splitIterator,
				  block, startBlock, endBlock,
				  strand)) {
		  newBlocks.push_back(startBlock);
		  block = endBlock;
		}
	      }*/
	    }
	  }
	  else {
	    if (((alignedContig->qryContig.strand == 0) and
		 (*splitIterator > block->qryBegin)) or
		((alignedContig->qryContig.strand == 1) and
		 (*splitIterator <= (qryLength - block->qryBegin - 1)))) {
	      // split according to position on forward strand
	      ssize_t splitResult;
	      if (alignedContig->qryContig.strand == 0)
		splitResult = SplitBlock(QUERY, *splitIterator,
					 block, startBlock, endBlock,
					 alignedContig->qryContig.strand);
	      else
		splitResult = SplitBlockRev(QUERY, qryLength - *splitIterator - 1,
					    block, startBlock, endBlock,
					    alignedContig->qryContig.strand);
	      if (splitResult) {
		/*
		std::cout << " split block into "
			  << startBlock << " and " << endBlock << std::endl;
		std::cout << "pushing back " << block << std::endl;
		*/
		newBlocks.push_back(startBlock);
		block = endBlock;
	      }
	    }
	  }
	}
	       
	assert(endBlock != NULL);
	assert(alignedContig->alignments[b] != block);
	alignedContig->alignments.erase(alignedContig->alignments.begin() + b);
	b--;
	newBlocks.push_back(block);
      }
    }
    // Copy all the new blocks

    ssize_t i;
    for (i = 0; i < newBlocks.size(); i++) {
      alignedContig->alignments.push_back(newBlocks[i]);
    }
    newBlocks.clear();
  }
}

/***************************************************************************
 * Title:          AxtLib.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  09/06/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "AxtLib.h"
#include <ostream>

void AxtToLAV(AxtEntries &axtList, LAVFile &lavFile, ssize_t refSeqLen, ssize_t qrySeqLen) { 

  ssize_t i, j;
  AxtEntry *axt;
  LAVBlock *lavBlock;
  LAVAlignedContig *lavAlignedContig;
  ssize_t contigExists;
  for (i = 0; i < axtList.size(); i++) {
    axt = axtList[i];
    contigExists = 0;
    for (j = 0; j < lavFile.alignments.size() && ! contigExists; j++) {
      if (lavFile.alignments[j]->qryContig.strand == axt->qryStrand) {
	lavAlignedContig = lavFile.alignments[j];
	contigExists = 1;
      }
    }
    if (contigExists == 0) {
      lavAlignedContig = new LAVAlignedContig;
      lavAlignedContig->refContig.start  = 1;
      lavAlignedContig->refContig.end    = refSeqLen;
      lavAlignedContig->refContig.strand = 0;
      lavAlignedContig->refContig.contig = 0;
      lavAlignedContig->qryContig.start  = 1; 
      lavAlignedContig->qryContig.end    = qrySeqLen;
      lavAlignedContig->qryContig.strand = axt->qryStrand;
      lavAlignedContig->qryContig.contig = 0;
      lavAlignedContig->refContigName = "\"" + axt->refTitle + "\"";
      lavAlignedContig->qryContigName = "\"" + axt->qryTitle + "\"";
	 
      lavFile.alignments.push_back(lavAlignedContig);
    }								      


    lavBlock = new LAVBlock;
    lavAlignedContig->alignments.push_back(lavBlock);

    // Create lav ungapped blocks from gap positions in axt 
    // file.
    ssize_t r, q;
    ssize_t refPos, qryPos;
    ssize_t matchLength = -1, gapLength = -1;
    lavBlock->refBegin = axt->refStart;
    lavBlock->refEnd   = axt->refEnd;
    lavBlock->qryBegin = axt->qryStart;
    lavBlock->qryEnd   = axt->qryEnd;

    assert(axt->refGapLocations.size() >= 2);
    assert(axt->qryGapLocations.size() >= 2);

    r = 0; q = 0;
    refPos = axt->refGapLocations[r]; 
    qryPos = axt->qryGapLocations[q];
    lavBlock->refALBegin.push_back(refPos + axt->refStart);
    lavBlock->qryALBegin.push_back(qryPos + axt->qryStart);
    //    r++; q++;
    ssize_t gappedQry = 0;
    ssize_t gappedRef = 0;
    ssize_t alignIndex = refPos;

    while (r < axt->refGapLocations.size() &&
	   q < axt->qryGapLocations.size()) {

      // Loop invariant: 
      //   At the beginning of each loop, refPos and qryPos point 
      //   to the beginning of an ungapped alignment.
      // also: 
      //   r and q index the beginning of a gap. 
      
      if (r < axt->refGapLocations.size()-1 && 
	  axt->refGapLocations[r+1] < axt->qryGapLocations[q+1]) {
	// The reference sequence is the next one to have a gap.
	matchLength = axt->refGapLocations[r+1] - alignIndex;
	//	gapLength   = axt->refGapLocations[r + 1] - axt->refGapLocations[r];
	r+=2;  // Advance r to the beginning of the next gap 
	gappedQry = 0;
	gappedRef = 1;
      }
      else if (q < axt->qryGapLocations.size() - 1) {
	// The query sequence is the next one to have a gap.
	matchLength = axt->qryGapLocations[q+1] - alignIndex;
	//	gapLength = axt->qryGapLocations[q+1] - axt->qryGapLocations[q];
	q += 2; // advance to the beginning of the next gap
	gappedQry = 1;
	gappedRef = 0;
      }

      lavBlock->refALEnd.push_back(refPos + matchLength + axt->refStart);
      lavBlock->qryALEnd.push_back(qryPos + matchLength + axt->qryStart);
      lavBlock->alIdentity.push_back(100); // Actually compute score later?
      if (gappedQry && q < axt->qryGapLocations.size()) {
	gapLength = axt->qryGapLocations[q] - axt->qryGapLocations[q-1];
	refPos += (matchLength + gapLength);
	qryPos += (matchLength + 1);        
      }
      else if (gappedRef &&  r < axt->refGapLocations.size() ) {
	refPos += (matchLength + 1);	 
	gapLength = axt->refGapLocations[r] - axt->refGapLocations[r-1];
	qryPos += (matchLength + gapLength );

      }
      alignIndex += matchLength + gapLength;
      if (r < axt->refGapLocations.size() - 1 &&
	  q < axt->qryGapLocations.size() - 1) {
	lavBlock->refALBegin.push_back(refPos + axt->refStart);
	lavBlock->qryALBegin.push_back(qryPos + axt->qryStart);
      }
    }
  }
}


void PrintAxt(AxtEntries &axtFile, std::ostream &out) {
  //UNUSED+// ssize_t j;
  ssize_t i ;
  for (i = 0; i < axtFile.size(); i++) {
    axtFile[i]->Print(out);
  }
}

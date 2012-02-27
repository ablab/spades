/***************************************************************************
 * Title:          LAVUtils.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef _INV_UTILS_
#define _INV_UTILS_

#include "lav/LAVTable.h"
#include "lav/LAVFile.h"
#include "lav/LAVBlock.h"
#include "lav/LAVAlignedContig.h"


// functors for sorting blocks

class BlockReferenceOrder {
public:
  ssize_t operator()(LAVBlock* a, LAVBlock *b) {
    if (a->refBegin < b->refBegin)
      return 1;
    // sort everything that is equal on ending position
    else if (a->refBegin == b->refBegin and a->refEnd < b->refEnd) 
      return 1;
    else
      return 0;
  }
};

class BlockQueryOrder {
public:
  ssize_t length;
  ssize_t operator()(LAVBlock* a, LAVBlock *b) {
    ssize_t aStart, bStart;
    if (a->strand == 1) {
      aStart = length - a->qryBegin + 1;
    }
    else {
      aStart = a->qryBegin;
    }
    if (b->strand == 1) {
      bStart = length - b->qryBegin + 1;
    }
    else {
      bStart = b->qryBegin;
    }
    return aStart < bStart;
  }
};



ssize_t ScoreBlock(LAVBlock &block,
		 ssize_t strand,
		 DNASequence &refSeq,
		 DNASequence &qrySeq,
		 IntMatrix &scoreMat,
		 ssize_t gapOpen,
		 ssize_t gapExtend);

ssize_t DetermineOrientation(LAVFile *lavFile);

ssize_t FindComplementaryAlignment(LAVAlignedContig *alignment, 
			       LAVFile &lavFile, 
			       LAVAlignedContig *&compAlignment);

ssize_t DetermineOrientation(LAVAlignedContig *forwardContig,
			 LAVAlignedContig *reverseContig);


void ComputeAverageAlignmentScore(LAVBlock *block,
				  DNASequence &refSeq,
				  DNASequence &qrySeq,
				  ssize_t windowSize,
				  IntMatrix &scoreMat,
				  ssize_t gapOpen,
				  ssize_t gapExtend,
				  IntVector &avgAlignScore);


void ComputePointAlignmentScore(LAVBlock *block,
				DNASequence &refSeq,
				DNASequence &qrySeq,
				IntMatrix &scoreMat,
				ssize_t gapOpen,
				ssize_t gapExtend,
				IntVector &alignScore);



void ComputeAverageAlignmentScore(LAVFile &lavFile,
				  DNASequence &refSeq,
				  DNASequence &qrySeq,
				  DNASequence &qrySeqRC,
				  ssize_t windowSize,
				  IntMatrix &scoreMat,
				  ssize_t gapOpen,
				  ssize_t gapExtend,
				  IntVector &avgAlignScores);


void ComputeTotalAlignmentScore(LAVFile &lavFile,
				DNASequence &refSeq,
				DNASequence &qrySeq,
				DNASequence &qrySeqRC,
				IntMatrix &scoreMat, 
				ssize_t gapOpen,
				ssize_t gapExtend,
				IntVector &alignScore);

void GetBlockStatistics(LAVBlock &block,
			ssize_t strand,
			DNASequence &refSeq, DNASequence &qrySeq,
			ssize_t &numMatch, ssize_t &numMismatch, ssize_t &numGap);

void BeginAtZero(LAVFile &lavFile);
void BeginAtOne(LAVFile &lavFile);
void ChangeStart(LAVFile &lavFile, ssize_t offset);

void FlipReverseCoordinates(LAVFile &lavFile);

template<typename Order_T>
void StoreBlockArray(LAVFile &lavFile, 
		     std::vector<LAVBlock*> &blocks, 
		     Order_T &order) {
  LAVAlignedContig *alignedContig;
  LAVBlock *block;
  ssize_t ac, b, i;
  ssize_t nBlocks = 0;
  for (ac = 0; ac < lavFile.alignments.size(); ac++) {
    alignedContig = lavFile.alignments[ac];
    nBlocks += alignedContig->size();
  }
  blocks.resize(nBlocks);
  ssize_t bnum = 0;
  for (ac = 0; ac < lavFile.alignments.size(); ac++) {
    alignedContig = lavFile.alignments[ac];
    for (b = 0; b < alignedContig->size(); b++ )
      blocks[bnum++] = alignedContig->alignments[b];
  }
  sort(blocks.begin(), blocks.end(), order);
}

ssize_t FindPosition(ssize_t position, LAVBlock &block);

ssize_t FindBlock(ssize_t &position,
	      std::vector<LAVBlock *> &blocks,
	      ssize_t &curIndex,
	      ssize_t sequence,
	      ssize_t &contained);

#endif

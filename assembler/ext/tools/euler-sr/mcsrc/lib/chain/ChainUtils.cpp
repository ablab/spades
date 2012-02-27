/***************************************************************************
 * Title:          ChainUtils.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <vector>
#include "ChainUtils.h"

#include "lav/LAVAlignedContig.h"

void ChainsToAlignFile(std::vector<Chain*> &chains, LAVFile &lavFile) {
  // first add all chains in forward direction
  ssize_t c;
  ssize_t strand;
  LAVAlignedContig *alignedContig;
  lavFile.blastzOpts = LAVFile::standardBlastzOpts;
  for (strand = 1; strand >= 0; strand--) {
    alignedContig = new LAVAlignedContig;
    lavFile.alignments.push_back(alignedContig);
    alignedContig->refContig.strand = 0;
    alignedContig->qryContig.strand = strand;
    for (c = 0; c < chains.size(); c++ ) {
      if (chains[c]->header.qStrand == strand) {
	ChainToAlignment(*chains[c], *alignedContig);
      }
    }
  }
}

void ChainToAlignment(Chain &chain, LAVAlignedContig &alignment) {
  ssize_t c;
  ssize_t refPos = chain.header.tStart + 1;
  ssize_t qryPos = chain.header.qStart + 1;

  ssize_t increment;
  if (chain.header.qStrand == 0) {
    increment = 1;
  }
  else {
    increment = -1;
  }
  //UNUSED+// ssize_t blockRefEnd,blockQryEnd;
  ssize_t blockRefStart,  blockQryStart ;
  LAVBlock *block;
  for (c = 0; c < chain.numAlign(); c++ ) {
    block = new LAVBlock;
    alignment.alignments.push_back(block);
    blockRefStart = refPos;
    blockQryStart = qryPos;
    block->refBegin = blockRefStart;
    block->qryBegin = blockQryStart;
    while (c < chain.numAlign() and 
	   (chain.dt[c] == 0 or chain.dq[c] == 0)) {
      block->AddInterval(refPos, refPos + chain.size[c] - 1 ,
			 qryPos, qryPos + chain.size[c] - 1,
			 100);
      refPos += (chain.size[c] + chain.dt[c]);
      qryPos += (chain.size[c] + chain.dq[c]);
      c++;
    }
    if (c < chain.numAlign()) {
      // Add the last interval that is succeeded by a 
      // gap in both sequences.
      if (chain.size[c] > 0) {
	block->AddInterval(refPos, refPos + chain.size[c] - 1,
			   qryPos, qryPos + chain.size[c] - 1,
			   100);
      }
      // this block includes an interval
      block->refEnd = refPos + chain.size[c] - 1;
      block->qryEnd = qryPos + chain.size[c] - 1;
      // gap in both sequences, so next block starts after
      // the gap
      refPos += chain.size[c] + chain.dt[c];
      qryPos += chain.size[c] + chain.dq[c];
    }
    else {
      // the chain is over, end the block here as well.
      if (chain.numAlign() == 0) {
	block->refEnd = refPos;
	block->qryEnd = qryPos;
      }
      else {
	block->refEnd = refPos - 1;
	block->qryEnd = qryPos - 1;
      }
    }
  }
}

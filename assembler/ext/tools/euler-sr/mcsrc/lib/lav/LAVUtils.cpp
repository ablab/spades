/***************************************************************************
 * Title:          LAVUtils.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "LAVUtils.h"
#include "SeqUtils.h"
#include "align/alignutils.h"

ssize_t FindComplementaryAlignment(LAVAlignedContig *alignment, LAVFile &lavFile,
															 LAVAlignedContig *&compAlign) {
  // Look through the alignments for one that is of opposite
  // orientation, and has the same name except with (orwithout)
  // (reverse complement)

  std::string targetName;
  std::string revCompStr = " (reverse complement)";
  targetName = alignment->qryContigName;
  if (alignment->qryContig.strand == 0) {
    targetName += " (reverse complement)";
  }
  else {
    ssize_t pos;
    pos = targetName.find(revCompStr);
    if (pos != targetName.npos) {
      // current string is reverse complement, look for the forward
      // strand
      targetName = targetName.replace(pos, revCompStr.size(), "");
    }
  }
  ssize_t i;
  
  // Determined the target string name, now find it 
  // in all of the alignments.
  for (i = 0; i < lavFile.size(); i++ ) {
    if (lavFile.alignments[i]->qryContigName == targetName) {
      compAlign = lavFile.alignments[i];
      return i;
    }
  }
  compAlign = NULL;
  return -1;
}

ssize_t DetermineOrientation(LAVAlignedContig *forwardAlign, 
												 LAVAlignedContig *reverseAlign) {
  //UNUSED// ssize_t i;
  //UNUSED+// ssize_t ac,rcb;
  ssize_t a  ;
  //UNUSED// LAVBlock *block;
  std::vector<ssize_t> rcBegin, rcEnd;
  ssize_t forwardScore = 0;
  ssize_t reverseScore = 0;
  // Find out where all of the 
  
  if (reverseAlign != NULL)
    for (a = 0; a < reverseAlign->alignments.size(); a++ ) {
      reverseScore += reverseAlign->alignments[a]->score;
    }

  if (forwardAlign != NULL )
    for (a = 0; a < forwardAlign->alignments.size(); a++ ) {
      forwardScore += forwardAlign->alignments[a]->score;
    }

  if (forwardScore > reverseScore) return 0;
  return 1;

}

		      

ssize_t DetermineOrientation(LAVFile *lavFile) {

  //UNUSED// ssize_t i;
  //UNUSED+// ssize_t rcb;
  ssize_t a, ac ;
  LAVAlignedContig *alignedContig;
  //UNUSED// LAVBlock *block;
  std::vector<ssize_t> rcBegin, rcEnd;
  ssize_t forwardScore = 0;
  ssize_t reverseScore = 0;
  // Find out where all of the 
  for (ac = 0; ac < lavFile->alignments.size(); ac++) {
    alignedContig = lavFile->alignments[ac];
    if (alignedContig->qryContig.strand == 1) {
      for (a = 0; a < alignedContig->alignments.size(); a++ ) {
				reverseScore += alignedContig->alignments[a]->score;
      }
    }
    else {
      for (a = 0; a < alignedContig->alignments.size(); a++ ) {
				forwardScore += alignedContig->alignments[a]->score;
      }
    }
  }

  if (forwardScore > reverseScore) return 0;
  return 1;
}

#if 0
// Disabled - function appears to be incomplete
ssize_t ComputeBlockScore(LAVBlock &block,
											ssize_t blockPos,
											DNASequence &refSeq,
											DNASequence &qrySeq,
											IntMatrix &scoreMat) {

  ssize_t pos;
  ssize_t score;
  ssize_t refPos, qryPos;
  ssize_t blockLength;
  blockLength = block.refALEnd[blockPos] 
    - block.refALBegin[blockPos] + 1;
  for (pos = 0; pos < blockLength; pos++) {
    refPos=  block.refALBegin[refPos] + pos;
    qryPos = block.refALEnd[refPos] + pos;
  }
}
#endif


void ComputePointAlignmentScore(LAVBlock *block,
																DNASequence &refSeq,
																DNASequence &qrySeq,
																IntMatrix &scoreMat,
																ssize_t gapOpen,
																ssize_t gapExtend,
																IntVector &alignScore) {
  // Store the score at every position in the sequence. 
  // That uses some extra space, but it makes computing the mean score
  // a cinch.
  ssize_t pos;
  ssize_t curBlock = 0;
  ssize_t indexInBlock = 0;
  ssize_t refPos, qryPos;
  //UNUSED// ssize_t prevInGap = 0;
  for (pos = block->refBegin; pos < block->refEnd ; pos++) {
    // Look to see if the position is in a gap
    // Quick sanity check - if pos is past the end of the 
    // current block, we expect it to be 
    assert(pos <= block->refALEnd[curBlock] or curBlock < block->size());

    if (pos > block->refALEnd[curBlock] and
				pos < block->refALBegin[curBlock + 1]) {
      /*      alignScore[pos - block->refBegin] += gapOpen + 
							gapExtend * (pos - block->refALEnd[curBlock]); */
      alignScore[pos] += gapOpen + 
				gapExtend * (pos - block->refALEnd[curBlock]);

    }
    else {
      if (curBlock < block->size()-1 and 
					pos == block->refALBegin[curBlock+1]) {
				curBlock++;
				indexInBlock = 0;
      }
      // The blocks are not gapped, so we can just use the same offset 
      // into the block for each block.
      refPos = block->refALBegin[curBlock]  + indexInBlock - 1;
      qryPos = block->qryALBegin[curBlock]  + indexInBlock - 1;
      ssize_t refNuc, qryNuc;
      refNuc = nuc_index[refSeq[refPos]];
      qryNuc = nuc_index[qrySeq[qryPos]];
      /*
				alignScore[pos - block->refBegin] += scoreMat[refNuc][qryNuc];
      */
      alignScore[pos] += scoreMat[refNuc][qryNuc];
      indexInBlock++;
    }
  }
}


void ComputeTotalAlignmentScore(LAVFile &lavFile,
																DNASequence &refSeq,
																DNASequence &qrySeq,
																DNASequence &qrySeqRC,
																IntMatrix &scoreMat,
																ssize_t gapOpen,
																ssize_t gapExtend,
																IntVector &alignScore) {
  
  LAVAlignedContig *alignedContig;
  LAVBlock *block;
  alignScore.resize(refSeq.length);
  ssize_t i;

  ssize_t ac, b;
  for (i = 0; i < refSeq.length; i++) 
    alignScore[i] = 0;
  for (ac = 0; ac < lavFile.alignments.size(); ac++) {
    alignedContig = lavFile.alignments[ac];
    for (b = 0; b < alignedContig->size(); b++) {
      block = alignedContig->alignments[b];
      if (alignedContig->qryContig.strand == 0)
				ComputePointAlignmentScore(block, 
																	 refSeq, qrySeq,     
																	 scoreMat,
																	 gapOpen, gapExtend,alignScore);
      else
				ComputePointAlignmentScore(block, 
																	 refSeq, qrySeqRC,
																	 scoreMat,
																	 gapOpen, gapExtend,alignScore);
	
    }
  }
}

void ComputeAverageAlignmentScore(LAVFile &lavFile,
																	DNASequence &refSeq,
																	DNASequence &qrySeq,
																	DNASequence &qrySeqRC,
																	ssize_t windowSize,
																	IntMatrix &scoreMat,
																	ssize_t gapOpen,
																	ssize_t gapExtend,
																	IntVector &avgAlignScore) {


  IntVector totalAlignScore;
  ComputeTotalAlignmentScore(lavFile, refSeq, qrySeq, qrySeqRC,
														 scoreMat, gapOpen, gapExtend, 
														 totalAlignScore);
  ssize_t windowScore = 0;
  ssize_t i;
  // some error checking
  if (windowSize > refSeq.length)
    return;
  // window size
  for (i = 0; i < windowSize; i++) {
    windowScore += totalAlignScore[i];
  }
  ssize_t averageLength = refSeq.length - windowSize;
  for (i = 0; i < averageLength; i++ ) {
    avgAlignScore[i] = windowScore / windowSize;
    windowScore -= totalAlignScore[i];
    windowScore += totalAlignScore[i+windowSize];
  }
}


void ComputeAverageAlignmentScore(LAVBlock *block,
																	DNASequence &refSeq,
																	DNASequence &qrySeq,
																	ssize_t windowSize,
																	IntMatrix &scoreMat,
																	ssize_t gapOpen,
																	ssize_t gapExtend,
																	IntVector &avgAlignScores) {
  //UNUSED// ssize_t b;
  //UNUSED// ssize_t alignPos;
  //UNUSED// ssize_t start, end;

  ssize_t blockLength = block->refEnd - block->refBegin + 1;
  avgAlignScores.resize(blockLength);
  ssize_t pos;
  ssize_t curBlock = 0;
  ssize_t indexInBlock;
  ssize_t firstHalf, secondHalf;
  firstHalf = windowSize / 2;
  secondHalf = windowSize - firstHalf;
  ssize_t windowStart, windowEnd;

  IntVector refSeqScores;
  // Create a cached version of the scores of each nucleotide
  curBlock = 0;
  indexInBlock = 0;
  refSeqScores.resize(blockLength);
  avgAlignScores.resize(blockLength);
  ssize_t i;
  for (i = 0; i < blockLength; i++) {
    refSeqScores[i] = 0;
    avgAlignScores[i] = 0;
  }
    
  ComputePointAlignmentScore(block, refSeq, qrySeq,
														 scoreMat, gapOpen, gapExtend, 
														 refSeqScores);
  ssize_t totalScore = 0;
  for (pos = block->refBegin; pos < block->refEnd; pos++) {
    // Compute the starting and ending position of the average
    if (pos - firstHalf < block->refBegin) {
      windowStart = block->refBegin;
      windowEnd = windowSize - (pos - windowStart) + 1;
    }
    else {
      windowStart = pos - firstHalf;
      windowEnd   = pos + secondHalf;
    }
    
    if (pos + secondHalf > block->refEnd or 
				windowEnd > block->refEnd) {
      windowEnd = block->refEnd + 1;
    }
    else {
      windowEnd  = pos + secondHalf + 1;
    }
    totalScore = 0;
    ssize_t w;
    for (w = windowStart; w < windowEnd; w++) {
      totalScore += refSeqScores[w - block->refBegin];
    }    
    avgAlignScores[pos - block->refBegin] = 
      totalScore / (windowEnd - windowStart);
    // compute the average score at this position taking 
    // into account gaps
  }
}


void GetBlockStatistics(LAVBlock &block,
												ssize_t strand,
												DNASequence &refSeq, DNASequence &qrySeq,
												ssize_t &numMatch, ssize_t &numMismatch, ssize_t &numGap) {
  ssize_t r, q, i, fq;
  numMatch    = 0;
  numMismatch = 0;
  numGap      = 0;
  ssize_t ri, qi;
  ssize_t maxGap;
  for (i = 0; i < block.size(); i++) {
    q = block.qryALBegin[i];
    for (r = block.refALBegin[i]; r <= block.refALEnd[i]; r++, q++) {
      if (strand == 0) {
				ri = nuc_index[refSeq[r-1]];
				qi = nuc_index[qrySeq[q-1]];
      }
      else {
				fq = qrySeq.length - q + - 1;
				ri = nuc_index[refSeq[r-1]];
				qi = comp_bin[(unsigned char) nuc_index[qrySeq[fq-1]]];
      }
      if (ri == qi) numMatch++;
      else numMismatch++;
    }
    if (i < block.size()- 1) {
			// TODO: check this code
			//   Let DR = block.refALBegin[i+1] - block.refALEnd[i]
			//       DQ = block.qryALBegin[i+1] - block.qryALEnd[i]
			//   Code doesn't seem to handle case DR==DQ
			//   If DR==DQ should be put in, then code simplifies to maxGap = max(max(DR,DQ)-1,0)
			maxGap = 0;
      if ((block.refALBegin[i+1] - block.refALEnd[i] > 1)
					and (block.refALBegin[i+1] - block.refALEnd[i] > 
							 block.qryALBegin[i+1] - block.qryALEnd[i]))
				maxGap = block.refALBegin[i+1] - block.refALEnd[i]-1;
      if ((block.qryALBegin[i+1] - block.qryALEnd[i] > 1)
					and (block.qryALBegin[i+1] - block.qryALEnd[i] > 
							 block.refALBegin[i+1] - block.refALEnd[i]))
				maxGap = block.qryALBegin[i+1] - block.qryALEnd[i]-1;
      numGap+= maxGap;
    }
  }
}

ssize_t ScoreBlock(LAVBlock &block,
							 ssize_t strand,
							 DNASequence &refSeq,
							 DNASequence &qrySeq,
							 IntMatrix &scoreMat,
							 ssize_t gapOpen,
							 ssize_t gapExtend) {

  ssize_t r, q, i, fq;
  ssize_t maxGap;
  ssize_t score = 0;
  ssize_t ri, qi;
  ssize_t sc;

 
  for (i = 0; i < block.size(); i++) {
    q = block.qryALBegin[i];
    for (r = block.refALBegin[i]; r <= block.refALEnd[i]; r++, q++) {
      if (strand == 0) {
				ri = nuc_index[refSeq[r-1]];
				qi = nuc_index[qrySeq[q-1]];
				sc = scoreMat[ri][qi];
				score += sc;
      }
      else {
				fq = qrySeq.length - q;
				ri = nuc_index[refSeq[r-1]];
				qi = comp_bin[(unsigned char) nuc_index[qrySeq[fq]]];
				sc = scoreMat[ri][qi];
				score += sc;
      }
    }
    if (i < block.size()- 1) {
      maxGap = 0;
      if ((block.refALBegin[i+1] - block.refALEnd[i] > 1)
					and (block.refALBegin[i+1] - block.refALEnd[i] > 
							 block.qryALBegin[i+1] - block.qryALEnd[i]))
				maxGap = block.refALBegin[i+1] - block.refALEnd[i]-1;
      if ((block.qryALBegin[i+1] - block.qryALEnd[i] > 1)
					and (block.qryALBegin[i+1] - block.qryALEnd[i] > 
							 block.refALBegin[i+1] - block.refALEnd[i]))
				maxGap = block.qryALBegin[i+1] - block.qryALEnd[i]-1;

      score -= maxGap*gapExtend + (maxGap > 0 ? 1 : 0) * gapOpen;
    }
  }
  return score;
}

void ChangeStart(LAVFile &lavFile, ssize_t offset) {
  LAVAlignedContig *alignedContig;
  LAVBlock *block;
  ssize_t ac, b, i;
  for (ac = 0; ac < lavFile.alignments.size(); ac++) {
    alignedContig = lavFile.alignments[ac];
    if (alignedContig->refContig.start == 1) {
      for (b = 0; b < alignedContig->size(); b++) {
				block = alignedContig->alignments[b];
				block->refBegin += offset;
				block->refEnd   += offset;
				for (i = 0; i < block->size(); i++) {
					block->refALBegin[i] += offset;
					block->refALEnd[i] += offset;
				}
      }
    }
    if (alignedContig->qryContig.start == 1) {
      for (b = 0; b < alignedContig->size(); b++) {
				block = alignedContig->alignments[b];
				block->qryBegin += offset;
				block->qryEnd   += offset;
				for (i = 0; i < block->size(); i++) {
					block->qryALBegin[i] += offset;
					block->qryALEnd[i] += offset;
				}
      }
    }      
  }
}

void FlipReverseCoordinates(LAVFile &lavFile) {
  LAVAlignedContig *alignedContig;
  LAVBlock *block;
  ssize_t ac, b, i;
  ssize_t qryLength;
  for (ac = 0; ac < lavFile.alignments.size(); ac++) {
    alignedContig = lavFile.alignments[ac];
    if (alignedContig->qryContig.strand == 1) {
      qryLength = alignedContig->qryContig.end - alignedContig->qryContig.start + 1;
      for (b = 0; b < alignedContig->size(); b++) {
				block = alignedContig->alignments[b];
				block->qryBegin = qryLength - block->qryEnd - 1;
				block->qryEnd   = qryLength - block->qryBegin - 1;
				ssize_t last;
				last = block->size() / 2 + block->size() % 2;
				ssize_t bt, et;
				for (i = 0; i < last; i++) {
					bt = block->qryALBegin[i];
					et = block->qryALEnd[i];
					block->qryALBegin[i] = qryLength - block->qryALEnd[block->size() - i - 1] - 1;
					block->qryALEnd[i]   = qryLength - block->qryALBegin[block->size() - i - 1] - 1;
					block->qryALBegin[block->size() - i - 1] = qryLength - et - 1;
					block->qryALEnd[block->size() - i - 1]   = qryLength - bt - 1;
				}
      }
    }
  }
}

void BeginAtOne(LAVFile &lavFile) {
  ChangeStart(lavFile, 1);
}

void BeginAtZero(LAVFile &lavFile) {
  ChangeStart(lavFile, -1);
}

ssize_t FindBlock(ssize_t &position,
							std::vector<LAVBlock *> &blocks,
							ssize_t &curIndex,
							ssize_t sequence,
							ssize_t &contained) {
  // Make sure we can find the block at all (it is within bounds
  // of the alignment.
  curIndex = -1;
  if (blocks.size() == 0)
    return 0;
  if (position < blocks[0]->begin(sequence)) {
    curIndex = -1;
    return 0;
  }
  if (position > blocks[blocks.size()-1]->end(sequence)) {
    curIndex = blocks.size();
    return 0;
  }
  
  // do a binary search on the list of blocks to try and 
  // find either the block that contains 'position'
  // or is before 'position'.
  // some sanity checks
  ssize_t startIndex;
  ssize_t endIndex;

  startIndex = 0;
  endIndex   = blocks.size();
  curIndex   = endIndex / 2;

  assert(curIndex >= startIndex);
  assert(curIndex < endIndex);
  assert(curIndex < blocks.size());

  // We may have found the block.
  while (1) {
    if (blocks[curIndex]->begin(sequence) <= position and
				blocks[curIndex]->end(sequence)   >= position) {
      contained =1;
      return 1;
    }
    // gone too far
    if (curIndex > blocks.size() - 1)
      return 0;

    if ((curIndex < (blocks.size() - 1)) and 
				(position > blocks[curIndex]->end(sequence)) and
        (position < blocks[curIndex + 1]->begin(sequence))) {
      contained = 0;
      return 1;
    }
  
    if (position < blocks[curIndex]->begin(sequence)) {
      endIndex = curIndex;
      curIndex = (endIndex + startIndex) / 2;
    }
    else if (curIndex < blocks.size() - 1 and 
						 position >= blocks[curIndex+1]->begin(sequence)) {
      startIndex = curIndex;
      curIndex   = (endIndex + startIndex) / 2;
      assert (curIndex != startIndex);
    }
  }
  return 0;
}

ssize_t FindPosition(ssize_t position, LAVBlock &block) {
  //UNUSED// ssize_t i;
  // slow search, maybe replace with 
  ssize_t b, e, c;
  b = 0; e = block.size(); 
  c = e/2;
  while (1) {
    if ((position >= block.refALBegin[c] and
				 position <= block.refALEnd[c]) or
				((c < block.size() - 1) and
				 (block.refALBegin[c] <= position) and
				 (block.refALBegin[c+1] > position))) {
      break;
    }
    else if (position < block.refALBegin[c]) {
      e = c;
      c = (e + b) / 2;
    }
    else if (position >= block.refALEnd[c]) {
      b = c;
      c = (b + e) / 2;
    }
    else {
      std::cout << "reached a spot in binary search that shouldn't " 
								<< std::endl;
      assert(0);
    }
  }
  assert(c >= 0);
  assert(c < block.size());
  
  if (position <= block.refALEnd[c])
    // match in ungapped alignment
    // return exact position
    return block.qryALBegin[c] + (position - block.refALBegin[c]);
  else
    // match in gap in ungapped alignment
    // return closest position
    return block.qryALBegin[c];
}


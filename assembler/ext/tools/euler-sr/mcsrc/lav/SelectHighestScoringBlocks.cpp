/***************************************************************************
 * Title:          SelectHighestScoringBlocks.cpp 
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

ssize_t REFERENCE = 0;
ssize_t QUERY = 1;

double GetIdentityScore(DNASequence &seq, 
		       LAVBlock &block, Score &scoreMat,
		       ssize_t sequence);

void FilterSmallBlocks(LAVFile &lavFile, ssize_t size);

void FilterRemovedBlocks(LAVFile &lavFile, 
			 std::set<LAVBlock*> &removed);


    
int main(int argc, char* argv[]) {
  
  if (argc < 5) {
    std::cout << "usage: " << argv[0] 
	      << " inLavFile refseq qryseq outlavFile [scorethreshold] " 
	      << std::endl;
    exit(0);
  }

  std::vector<LAVBlock*> blocks;
  std::string inLavName  = argv[1];
  std::string refSeqName = argv[2];
  std::string qrySeqName = argv[3];
  std::string outLavName = argv[4];
  
  int argi = 5;
  ssize_t verbose = 0;
  while (argi < argc) {
    if (strcmp(argv[argi], "-verbose") == 0) {
      verbose = 1;
    }
    ++argi;
  }

  LAVFile lavFile;
  LAVReader::ReadLAVFile(inLavName, lavFile);

  DNASequence refSeq, qrySeq;
  SeqReader::GetSeq(refSeqName, refSeq, SeqReader::noConvert);
  SeqReader::GetSeq(qrySeqName, qrySeq, SeqReader::noConvert);
  char *src;
  src = getenv("MCSRC");
  std::string scorematName = 
    std::string(src) + "/align/data/blastzscoremat.txt"; 
  //  Score score(scorematName, 400,50);

  Score score; // use default scoring matrix
  std::map<char, double> keywordOptions;
  lavFile.ParseBlastzOpts(lavFile.blastzOpts, score.scoreMat, keywordOptions);

  // Compute the average score for estimating the 
  // re-scored divergence.
  double totalPenalty = 0.0;
  ssize_t i,j;
  for (i = 0; i < 4; i++) 
    for (j = 0; j < 4; j++) 
      if (i != j) totalPenalty+= (score.scoreMat[i][j] - score.scoreMat[i][i]);
  
  double averagePenalty = totalPenalty / 12.0;

  // Hard-wire the percent identity for now.
  double pctDivergence = 0.30; 
  
  // Also hard-wire the correction factor.
  double scoreEstimateErrorFactor = 1.0;

  // Remove blocks less than some size (10 here)
  FilterSmallBlocks(lavFile, 15);

  // count the number of blocks in the entire file
  ssize_t a, b, bi;
  ssize_t nBlocks = 0;
  for (a = 0; a < lavFile.size(); a++) {
    nBlocks+= lavFile.alignments[a]->size();
  }


  // Store references to all blocks in an array.
  LAVAlignedContig *alignment;
  blocks.resize(nBlocks);
  bi = 0;
  for (a = 0; a < lavFile.size(); a++) {
    alignment = lavFile.alignments[a];
    for (b = 0; b < alignment->size(); b++) {
      blocks[bi] = alignment->alignments[b];
      blocks[bi]->strand = alignment->qryContig.strand;
      bi++;
    }
  }

  BlockReferenceOrder referenceOrder;

  sort(blocks.begin(), blocks.end(), referenceOrder);
  std::vector<ssize_t> scores, lengths, validBlocks;
  ssize_t repeatStart;
  std::vector<double> blockScores, identityScores, expectedPenalties, measuredPenalties;
  
  double identityScore;
  ssize_t blockLength;
  ssize_t overlapStart, overlapQueryStart, repeatStartIndex;
  ssize_t overlapEnd;
  std::set<LAVBlock*> toRemove;

  ssize_t bo;
  for (b = 0; b < blocks.size()-1; b++) {
    bo = b+1;
    overlapStart = -1;
    while ( bo < blocks.size() and blocks[bo]->refBegin <= blocks[b]->refEnd){
      if (overlapStart == -1) {
	overlapQueryStart = blocks[b]->qryBegin;
	overlapStart  = blocks[b]->refBegin;
	overlapEnd    = blocks[b]->refEnd;
	repeatStartIndex = b;
	scores.clear();
	lengths.clear();
	validBlocks.clear();
	expectedPenalties.clear();
	measuredPenalties.clear();
	identityScores.clear();
	blockScores.clear();
	scores.push_back(blocks[b]->score);
	blockLength = blocks[b]->refEnd - blocks[b]->refBegin + 1;
	lengths.push_back(blockLength);
	identityScore = GetIdentityScore(refSeq, *blocks[b], score, REFERENCE);
	identityScores.push_back(identityScore);
	blockScores.push_back(blocks[b]->score);
	expectedPenalties.push_back(averagePenalty* blockLength*pctDivergence);
	measuredPenalties.push_back(scoreEstimateErrorFactor*
				    (blocks[b]->score - identityScore));
	if ((scoreEstimateErrorFactor * ( blocks[b]->score - identityScore )) >  
	    (averagePenalty* blockLength*pctDivergence)) {
	  validBlocks.push_back(b);
	}
      }
      scores.push_back(blocks[bo]->score);
      blockLength   = blocks[bo]->refEnd - blocks[bo]->refBegin + 1;
      identityScore = GetIdentityScore(refSeq, *blocks[bo], score, REFERENCE);
      identityScores.push_back(identityScore);
      blockScores.push_back(blocks[bo]->score);
      lengths.push_back(blockLength);
      expectedPenalties.push_back(averagePenalty* blockLength*pctDivergence);
      measuredPenalties.push_back(scoreEstimateErrorFactor * 
				  (blocks[bo]->score - identityScore )); 
      if ((scoreEstimateErrorFactor * (blocks[bo]->score - identityScore )) >   
	  (averagePenalty* blockLength*pctDivergence)) {
	validBlocks.push_back(bo);
      }
      if (verbose) {
	std::cout << "found reference overlap: " 
		  << "(" << b << ", " << blocks[b]->refBegin << ", " 
		  << blocks[b]->qryBegin << ") "
		  << " ... (" << blocks[b]->refEnd << ", "
		  << blocks[b]->qryEnd << ") " << blocks[b]->score << std::endl;
	std::cout << " with                    "
		  << "(" << bo << ", "<< blocks[bo]->refBegin << ", "
		  << blocks[bo]->qryBegin << ") "
		  << " ... (" << blocks[bo]->refEnd << ", "
		  << blocks[bo]->qryEnd << ")  " 
		  << blocks[bo]->score << " "
		  << std::endl;
      }
      /*
      ssize_t curMatch, curMismatch, curGap,
	nextMatch, nextMismatch, nextGap;
      GetBlockStatistics(*(blocks[b]), blocks[b]->strand, refSeq, qrySeq,
			 curMatch, curMismatch, curGap);
      GetBlockStatistics(*(blocks[bo]), blocks[bo]->strand, refSeq, qrySeq,
			 nextMatch, nextMismatch, nextGap);

      std::cout << "statistics: " << curMatch << " " 
		<< curMismatch << " " << curGap << std::endl;
      std::cout << "statistics: " << nextMatch << " " 
		<< nextMismatch << " " << nextGap << std::endl;
      */
      bo++;
    } // end while loop
    if (overlapStart != -1 and 
	scores.size() > 1) {
      ssize_t s;
      if (verbose) {
	std::cout << "position : " << overlapStart << " " << overlapEnd << " scores:";
	for (s = 0; s < scores.size(); s++) {
	  std::cout << " " << scores[s];
	}
	std::cout << std::endl;
	std::cout << "           " << overlapStart << " " << overlapEnd << " lengths:";
	for (s = 0; s < lengths.size(); s++) {
	  std::cout << " " << lengths[s];
	}
	std::cout << std::endl;
	std::cout << "           " << overlapStart << " " << overlapEnd << " exp penal:";
	for (s = 0; s < expectedPenalties.size(); s++) {
	  std::cout << " " << expectedPenalties[s];
	}
	std::cout << std::endl;
	std::cout << "           " << overlapStart << " " << overlapEnd << " meas penal:";
	for (s = 0; s < measuredPenalties.size(); s++) {
	  std::cout << " " << measuredPenalties[s];
	}

	std::cout << std::endl;
	std::cout << "           " << overlapStart << " " << overlapEnd << " identity :";
	for (s = 0; s < identityScores.size(); s++) {
	  std::cout << " " << identityScores[s];
	}

	std::cout << std::endl;
	std::cout << "           " << overlapStart << " " << overlapEnd << " block sc :";
	for (s = 0; s < blockScores.size(); s++) {
	  std::cout << " " << blockScores[s];
	}

	std::cout << std::endl;
	std::cout << "valid blocks: ";
	for (s = 0; s < validBlocks.size(); s++)
	  std::cout << "           " << validBlocks[s] << " ";
	std::cout << std::endl;
      }
      // Remove invalid blocks
      if (validBlocks.size() == 1) {
	ssize_t b1;
	/*
	std::cout << "checking " << b - repeatStartIndex + 1 << " blocks and keeping " 
		  << validBlocks[0] << std::endl;
	*/
	for (b1 = repeatStartIndex; b1 < bo; b1++) {
	  if (b1 != validBlocks[0]) {
	    toRemove.insert(blocks[b1]);
	  }
	}
      }
    }
    overlapStart = -1;
    
  }



  blocks.clear();
  nBlocks = 0;
  // Store references to all blocks in an array.
  for (a = 0; a < lavFile.size(); a++) {
    nBlocks+= lavFile.alignments[a]->size();
  }

  blocks.resize(nBlocks);
  bi = 0;
  for (a = 0; a < lavFile.size(); a++) {
    alignment = lavFile.alignments[a];
    for (b = 0; b < alignment->size(); b++) {
      blocks[bi] = alignment->alignments[b];
      blocks[bi]->strand = alignment->qryContig.strand;
      bi++;
    }
  }

  BlockQueryOrder queryOrder;
  queryOrder.length = qrySeq.length;
  sort(blocks.begin(), blocks.end(), queryOrder);
  ssize_t prevQryEnd, curQryBegin;
  for (b = 0; b < blocks.size()-1; b++) {
    bo = b+1;
    // determine the ending coordinate of overlap begin
    if (blocks[b]->strand == 0) {
      prevQryEnd = blocks[b]->qryEnd;
    }
    else {
      prevQryEnd = qrySeq.length - blocks[b]->qryBegin + 1;
    }
    while (1) {    
      // determine the beginning coordinate of the overlap end
      if (bo >= blocks.size())
	break;

      if (blocks[bo]->strand == 0) {
	curQryBegin = blocks[bo]->qryBegin;
      }
      else {
	curQryBegin = qrySeq.length - blocks[bo]->qryEnd + 1;
      }
      
      // Not overlapping, stop this loop
      if (curQryBegin > prevQryEnd) 
	break;
      
      if (overlapStart == -1) {
	overlapQueryStart = blocks[b]->qryBegin;
	overlapStart  = blocks[b]->qryBegin;
	overlapEnd    = blocks[b]->qryEnd;
	repeatStartIndex = b;
	scores.clear();
	lengths.clear();
	validBlocks.clear();
	expectedPenalties.clear();
	measuredPenalties.clear();
	identityScores.clear();
	blockScores.clear();
	scores.push_back(blocks[b]->score);
	blockLength = blocks[b]->qryEnd - blocks[b]->qryBegin + 1;
	lengths.push_back(blockLength);
	identityScore = GetIdentityScore(qrySeq, *blocks[b], score, QUERY);
	identityScores.push_back(identityScore);
	blockScores.push_back(blocks[b]->score);
	expectedPenalties.push_back(averagePenalty* blockLength*pctDivergence);
	measuredPenalties.push_back(scoreEstimateErrorFactor*
				    (blocks[b]->score - identityScore));
	if ((scoreEstimateErrorFactor * ( blocks[b]->score - identityScore)) >
	    (averagePenalty* blockLength*pctDivergence)) {
	  validBlocks.push_back(b);
	}
	if (verbose) { 
	  std::cout << "found query overlap: " 
		    << "(" << b << ", " << blocks[b]->refBegin << ", " 
		    << blocks[b]->qryBegin << ") "
		    << " ... (" << blocks[b]->refEnd << ", "
		    << blocks[b]->qryEnd << ") " << blocks[b]->score << std::endl;
	  std::cout << " with                    "
		    << "(" << bo << ", "<< blocks[bo]->refBegin << ", "
		    << blocks[bo]->qryBegin << ") "
		    << " ... (" << blocks[bo]->refEnd << ", "
		    << blocks[bo]->qryEnd << ")  " 
		    << blocks[bo]->score << " "
		    << std::endl;
	}
      }
      scores.push_back(blocks[bo]->score);
      blockLength  = blocks[bo]->qryEnd - blocks[bo]->qryBegin + 1;
      identityScore = GetIdentityScore(qrySeq, *blocks[bo], score, QUERY);
      identityScores.push_back(identityScore);
      blockScores.push_back(blocks[bo]->score);
      lengths.push_back(blockLength);
      expectedPenalties.push_back(averagePenalty* blockLength*pctDivergence);
      measuredPenalties.push_back(scoreEstimateErrorFactor * 
				  (blocks[bo]->score - identityScore )); 
      if ((scoreEstimateErrorFactor * (blocks[bo]->score - identityScore )) > 
	  (averagePenalty* blockLength*pctDivergence)) {
	validBlocks.push_back(bo);
      }

      /*
	ssize_t curMatch, curMismatch, curGap,
	nextMatch, nextMismatch, nextGap;
      GetBlockStatistics(*(blocks[b]), blocks[b]->strand, refSeq, qrySeq,
			 curMatch, curMismatch, curGap);
      GetBlockStatistics(*(blocks[bo]), blocks[bo]->strand, refSeq, qrySeq,
			 nextMatch, nextMismatch, nextGap);

      std::cout << "found reference overlap: " 
		<< "(" << b << ", " << blocks[b]->refBegin << ", " 
		<< blocks[b]->qryBegin << ") "
		<< " ... (" << blocks[b]->qryEnd << ", "
		<< blocks[b]->qryEnd << ") " << blocks[b]->score << std::endl;
      std::cout << " with                    "
		<< "(" << bo << ", "<< blocks[bo]->refBegin << ", "
		<< blocks[bo]->qryBegin << ") "
		<< " ... (" << blocks[bo]->qryEnd << ", "
		<< blocks[bo]->qryEnd << ")  " 
		<< blocks[bo]->score << " "
		<< std::endl;
      std::cout << "statistics: " << curMatch << " " 
		<< curMismatch << " " << curGap << std::endl;
      std::cout << "statistics: " << nextMatch << " " 
		<< nextMismatch << " " << nextGap << std::endl;
      */
      bo++;
    } // end while loop
    if (overlapStart != -1 and 
	scores.size() > 1) {
      ssize_t s;
      if (verbose) {
	std::cout << "position : " << overlapStart << " " << overlapEnd << " scores:";
	for (s = 0; s < scores.size(); s++) {
	  std::cout << " " << scores[s];
	}
	std::cout << std::endl;
	std::cout << "           " << overlapStart << " " << overlapEnd << " lengths:";
	for (s = 0; s < lengths.size(); s++) {
	  std::cout << " " << lengths[s];
	}
	std::cout << std::endl;
	std::cout << "           " << overlapStart << " " << overlapEnd << " exp penal:";
	for (s = 0; s < expectedPenalties.size(); s++) {
	  std::cout << " " << expectedPenalties[s];
	}
	std::cout << std::endl;
	std::cout << "           " << overlapStart << " " << overlapEnd << " meas penal:";
	for (s = 0; s < measuredPenalties.size(); s++) {
	  std::cout << " " << measuredPenalties[s];
	}

	std::cout << std::endl;
	std::cout << "           " << overlapStart << " " << overlapEnd << " identity :";
	for (s = 0; s < identityScores.size(); s++) {
	  std::cout << " " << identityScores[s];
	}

	std::cout << std::endl;
	std::cout << "           " << overlapStart << " " << overlapEnd << " block sc :";
	for (s = 0; s < blockScores.size(); s++) {
	  std::cout << " " << blockScores[s];
	}

	std::cout << std::endl;
	std::cout << "valid blocks: ";
	for (s = 0; s < validBlocks.size(); s++)
	  std::cout << "           " << validBlocks[s] << " ";
	std::cout << std::endl;
      }
      // Remove invalid blocks
      if (validBlocks.size() == 1) {
	ssize_t b1;
	/*
	std::cout << "checking " << b - repeatStartIndex + 1 << " blocks and keeping " 
		  << validBlocks[0] << std::endl;
	*/
	for (b1 = repeatStartIndex; b1 < bo; b1++) {
	  if (b1 != validBlocks[0]) {
	    toRemove.insert(blocks[b1]);
	  }
	}
      }
    }
    overlapStart = -1;
  }


  // Remove overlapping blocks from the alignment
  FilterRemovedBlocks(lavFile, toRemove);

  // Output the result
  std::ofstream lavOut;
  openck(outLavName, lavOut, std::ios::out);
  LAVPrinter::PrintLAVFile(lavFile, lavOut);
  return 0;
}


double GetIdentityScore(DNASequence &seq, 
		       LAVBlock &block,
		       Score &scoreMat, ssize_t sequence) {

  ssize_t p;
  double score = 0;
  ssize_t b;
  ssize_t rp;
  for (b = 0; b < block.size(); b++) {
    if (sequence == REFERENCE) {
      for (p = block.refALBegin[b]; p < block.refALEnd[b]; p++)
	score += scoreMat.scoreMat[nuc_index[seq.seq[p]]][nuc_index[seq.seq[p]]];
    }
    else {
      for (p = block.qryALBegin[b]; p < block.qryALEnd[b]; p++)
	if ( block.strand == 0 ) {
	  score += scoreMat.scoreMat[nuc_index[seq.seq[p]]][nuc_index[seq.seq[p]]];
	}
	else {
	  rp = seq.length - p + 1;
	  score += scoreMat.scoreMat[nuc_index[seq.seq[rp]]][nuc_index[seq.seq[rp]]];
	}
    }
  }
  return score;
}

void FilterSmallBlocks(LAVFile &lavFile, ssize_t size) {
  ssize_t a, b;
  LAVAlignedContig *alignment;
  LAVBlock *block;
  for (a = 0; a < lavFile.size(); a++) {
    alignment = lavFile.alignments[a];
    for (b = 0; b < alignment->size(); ) {
      block = alignment->alignments[b];
      if (block->refEnd - block->refBegin + 1 < size) {
	alignment->alignments.erase(alignment->alignments.begin() + b);
      }
      else {
	b++;
      }
    }
  }
}

void FilterRemovedBlocks(LAVFile &lavFile, 
			 std::set<LAVBlock*> &removed) { 
  LAVAlignedContig *alignment;
  LAVBlock *block;
  ssize_t a, b;
  for (a = 0; a < lavFile.size(); a++) {
    alignment = lavFile.alignments[a];
    for (b = 0; b < alignment->size(); ) {
      block = alignment->alignments[b];
      if (removed.find(block) != removed.end()) {
	/*
	std::cout << "removing block (" << block->refBegin << ", " << block->qryBegin
		  << ") - (" << block->refEnd << " , " << block->qryEnd << " ) " 
		  << std::endl;
	*/
	alignment->alignments.erase(alignment->alignments.begin() + b);
      }
      else {
	b++;
      }
    }
  }
}

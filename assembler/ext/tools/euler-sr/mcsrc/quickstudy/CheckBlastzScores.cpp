/***************************************************************************
 * Title:          CheckBlastzScores.cpp 
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


int main(int argc, char* argv[]) {

  if (argc < 4) {
    std::cout << "usage: " << argv[0] << " inLavFile refSeq qrySeq  " << std::endl;
    exit(0);
  }
  std::string inLavName  = argv[1];
  std::string refSeqName = argv[2];
  std::string qrySeqName = argv[3];

  LAVFile lavFile;
  DNASequence refSeq, qrySeq;
  LAVReader::ReadLAVFile(inLavName, lavFile);
  SeqReader::GetSeq(refSeqName, refSeq);
  SeqReader::GetSeq(qrySeqName, qrySeq);
  std::map<char, ssize_t> keywordOptions;
  Score score; // use default scoring matrix
  std::cout << "parsing blastz options " << std::endl;
  lavFile.ParseBlastzOpts(lavFile.blastzOpts, score.scoreMat, keywordOptions);


  //UNUSED// int refBlockStart, refBlockEnd;
  //UNUSED// int qryBlockStart, qryBlockEnd;
  std::set<ssize_t> refCoords, qryCoords;
  ssize_t  a, b;
  LAVAlignedContig *alignedContig;
  LAVBlock *block;
  ssize_t refLength, qryLength;
  
  double blockScore;
  for (a = 0; a < lavFile.size(); a++) {
    alignedContig = lavFile.alignments[a];
    refLength = alignedContig->refContig.end - 
      alignedContig->refContig.start + 1;
    qryLength = alignedContig->qryContig.end - 
      alignedContig->qryContig.start + 1;

    for (b = 0; b < alignedContig->size(); b++) {
      block = alignedContig->alignments[b];

      blockScore = ScoreBlock(*block, alignedContig->qryContig.strand,
			 refSeq, qrySeq,
			 score.scoreMat, score.gapOpen, score.gapExtend);

      std::cout << blockScore << "\t" << block->score << std::endl;
    }
  }
}

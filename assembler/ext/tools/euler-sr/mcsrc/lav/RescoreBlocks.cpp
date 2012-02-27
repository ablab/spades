/***************************************************************************
 * Title:          RescoreBlocks.cpp 
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
#include "SeqUtils.h"

const ssize_t REFERENCE = 0;
const ssize_t QUERY     = 1;


int main(int argc, char* argv[]) {

  if (argc < 5) {
    std::cout << "usage: " << argv[0] << " inLavFile refSeq qrySeq outLavFile " << std::endl;
    exit(0);
  }
  std::string inLavName  = argv[1];
  std::string refSeqName = argv[2];
  std::string qrySeqName = argv[3];
  std::string outLavName = argv[4];

  LAVFile lavFile;
  DNASequence refSeq, qrySeq;
  LAVReader::ReadLAVFile(inLavName, lavFile);
  SeqReader::GetSeq(refSeqName, refSeq);
  SeqReader::GetSeq(qrySeqName, qrySeq);
  std::map<char, ssize_t> keywordOptions;
  Score score; // use default scoring matrix
  lavFile.ParseBlastzOpts(lavFile.blastzOpts, score.scoreMat, keywordOptions);
  

  //UNUSED// ssize_t refBlockStart, refBlockEnd;
  //UNUSED// ssize_t qryBlockStart, qryBlockEnd;
  std::set<ssize_t> refCoords, qryCoords;
  ssize_t  a, b;
  LAVAlignedContig *alignedContig;
  LAVBlock *block;
  ssize_t refLength, qryLength;
  
  //UNUSED// double blockScore;
  for (a = 0; a < lavFile.size(); a++) {
    alignedContig = lavFile.alignments[a];
    refLength = alignedContig->refContig.end - 
      alignedContig->refContig.start + 1;
    qryLength = alignedContig->qryContig.end - 
      alignedContig->qryContig.start + 1;


    for (b = 0; b < alignedContig->size(); b++) {
      block = alignedContig->alignments[b];
      block->score = (ssize_t) floor(ScoreBlock(*block, alignedContig->qryContig.strand,
																						refSeq, qrySeq,
																						score.scoreMat, score.gapOpen, score.gapExtend));
    }
  }

  std::ofstream lavOut;
  openck(outLavName, lavOut, std::ios::out);
  LAVPrinter::PrintLAVFile(lavFile, lavOut);
  lavOut.close();
}

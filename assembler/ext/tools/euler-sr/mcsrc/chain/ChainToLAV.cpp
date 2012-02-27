/***************************************************************************
 * Title:          ChainToLAV.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "chain/ChainUtils.h"
#include "chain/ChainReader.h"
#include "chain/Chain.h"
#include "lav/LAVFile.h"
#include "lav/LAVPrinter.h"
#include "DNASequence.h"
#include "SeqReader.h"

#include <string>
#include <vector>

int main(int argc, char* argv[]) {

  std::string chainFileName, refSeqName, qrySeqName, lavFileName;
 
  chainFileName = argv[1];
  refSeqName    = argv[2];
  qrySeqName    = argv[3];
  lavFileName   = argv[4];
  
  DNASequence refSeq, qrySeq;
  SeqReader::GetSeq(refSeqName, refSeq);
  SeqReader::GetSeq(qrySeqName, qrySeq);
  
  std::vector<Chain*> chains;
  ChainReader::ReadChainFile(chainFileName, chains);

  LAVFile lavFile;
  ChainsToAlignFile(chains,lavFile);

  ssize_t a;
  for (a = 0; a < lavFile.size(); a++) {
    lavFile.alignments[a]->refContig.contig = 1;
    lavFile.alignments[a]->qryContig.contig = 1;
    lavFile.alignments[a]->refContig.start = 1;
    lavFile.alignments[a]->refContig.end = refSeq.length;
    lavFile.alignments[a]->qryContig.start = 1;
    lavFile.alignments[a]->qryContig.end = qrySeq.length;
    lavFile.alignments[a]->refContigName = refSeq.namestr;
    lavFile.alignments[a]->qryContigName = qrySeq.namestr;
  }
  
  LAVPrinter::PrintLAVFile(lavFile, lavFileName);

  return 0;
}

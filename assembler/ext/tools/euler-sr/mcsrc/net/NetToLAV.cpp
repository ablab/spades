/***************************************************************************
 * Title:          NetToLAV.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "net/NetReader.h"

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

  std::string chainFileName, netFileName,
    refSeqName, qrySeqName, lavFileName,
    refLenName, qryLenName;

  if (argc != 8) {
    std::cout << "usage: netToLav in.chain in.net refName refLen qryName qryLen out.lav" << std::endl;
    exit(0);
  }
  chainFileName = argv[1];
  netFileName   = argv[2];
  refSeqName    = argv[3];
  refLenName    = argv[4];
  qrySeqName    = argv[5];
  qryLenName    = argv[6];
  lavFileName   = argv[7];
  
  ssize_t refSeqLength, qrySeqLength;
  std::ifstream refSeqLenIn, qrySeqLenIn;
  openck(refLenName, refSeqLenIn, std::ios::in);
  refSeqLenIn >> refSeqLength;
  refSeqLenIn.close();

  openck(qryLenName, qrySeqLenIn, std::ios::in);
  qrySeqLenIn >> qrySeqLength;
  qrySeqLenIn.close();
  
  std::vector<Chain*> chains;
  ChainReader::ReadChainFile(chainFileName, chains);

  NetFile netFile;

  NetReader::ReadNetFile(netFileName, netFile);
  std::vector<Chain*> firstLevelChains;

  ssize_t c, n;
  for (n = 0; n < netFile.size(); n++ ){ 
    if (netFile.nets[n]->level == 1) {
      for (c = 0; c < chains.size(); c++) {
	if (chains[c]->id == netFile.nets[n]->id) {
	  firstLevelChains.push_back(chains[c]);
	  break;
	}
      }
    }
  }

  LAVFile lavFile;
  ChainsToAlignFile(firstLevelChains, lavFile);

  ssize_t a;
  for (a = 0; a < lavFile.size(); a++) {
    lavFile.alignments[a]->refContig.contig = 1;
    lavFile.alignments[a]->qryContig.contig = 1;
    lavFile.alignments[a]->refContig.start = 1;
    lavFile.alignments[a]->refContig.end = refSeqLength;
    lavFile.alignments[a]->qryContig.start = 1;
    lavFile.alignments[a]->qryContig.sequenceName = qrySeqName + ".fasta";
    lavFile.alignments[a]->refContig.sequenceName = refSeqName + ".fasta";
    lavFile.alignments[a]->qryContig.end = qrySeqLength;
    lavFile.alignments[a]->refContigName = ">" + refSeqName;
    if (lavFile.alignments[a]->qryContig.strand == 0)
      lavFile.alignments[a]->qryContigName = ">" + qrySeqName;
    else 
      lavFile.alignments[a]->qryContigName = ">" + qrySeqName + " (reverse complement)";
  }
  
  LAVPrinter::PrintLAVFile(lavFile, lavFileName);

  return 0;
}

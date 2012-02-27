/***************************************************************************
 * Title:          FragmentUtils.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "FragmentUtils.h"
#include "lav/LAVUtils.h"
#include "utils.h"

ssize_t ReadLoci(std::string &locusFileName, Loci &loci) {
  std::ifstream in;
  openck(locusFileName, in, std::ios::in);
  Locus *locus;
  while (in) {
    locus = new Locus;
    if (!(in >> locus->seqid >> locus->start >> locus->end)) {
      delete locus;
      return 0;
    }
    loci.push_back(locus);
  }
  return 0;
}

ssize_t ReadLociSequences(Loci &loci, std::string &seqDir) {
  ssize_t l;
  DNASequence sequence;
  ssize_t locusLength;
  for (l = 0; l < loci.size(); l++ ) {
    std::stringstream nameStrm;
    nameStrm << seqDir << "/" << loci[l]->seqid << ".fasta";
    //    std::cout << "reading " << nameStrm.str() << std::endl;
    sequence.Reset();
    SeqReader::GetSeq(nameStrm.str(), sequence , SeqReader::noConvert  );
    loci[l]->locusSeq.Copy(sequence, loci[l]->start, loci[l]->end+1);
  }
}

ssize_t FindSurroundingRegion(ssize_t start, ssize_t end,
			  std::vector<LAVBlock*> &blocks,
			  ssize_t &regionStart, ssize_t &regionEnd, 
			  ssize_t seqLength) {
  regionStart = -1;
  regionEnd   = -1;

  // this should be a warning, but I don't want to 
  // deal with it.
  if (blocks.size() == 0)
    return 0;
  ssize_t contained;
  ssize_t blockIndex;
  if ( FindBlock(start, blocks, 
		 blockIndex, LAVBlock::ref, contained) ) {
    if (contained) {
      regionStart = FindPosition(start, *blocks[blockIndex]);
    }
    else 
      regionStart = blocks[blockIndex]->qryEnd + 1;
  }
  else  { 
    // some problem occurred, parse problems
    if (blockIndex == -1)
      regionStart = 0;
    if (blockIndex == blocks.size())
      regionStart = blocks[blocks.size()-1]->qryEnd + 1;
  }
    
  if (FindBlock(end, blocks, blockIndex, LAVBlock::ref, contained) ) {
    if (contained) {
      regionEnd = FindPosition(end, *blocks[blockIndex]);
    }
    else {
      regionEnd = blocks[blockIndex+1]->qryBegin - 1;
    }
  }
  else {
    if (blockIndex == -1) {
      regionEnd = blocks[0]->qryBegin-1;
    }
    else {
      regionEnd = seqLength;
    }
  }
  return (regionStart != -1 and regionEnd != -1);
}


void PrintTempSeq(DNASequence &sequence,
		  ssize_t sbjctStart, ssize_t sbjctEnd,
		  std::string &sbjctName) {
  DNASequence sbjctSeq;
  sbjctSeq.seq = &sequence.seq[sbjctStart];
  sbjctSeq.length = sbjctEnd - sbjctStart + 1;
  sbjctSeq.CopyDetails(sequence);
  std::ofstream sbjctOut;
  openck(sbjctName, sbjctOut, std::ios::out);
  sbjctSeq.PrintSeq(sbjctOut);
  sbjctOut << std::endl;
  sbjctOut.close();
}

    
void PrintLoci(Loci &loci, std::string &outName) {
  std::ofstream out;
  openck(outName, out, std::ios::out);
  ssize_t i;
  for ( i = 0; i < loci.size(); i++) {
    out << loci[i]->seqid << "\t" << loci[i]->start << "\t" << loci[i]->end << std::endl;
  }
  out.close();
}

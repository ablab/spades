/***************************************************************************
 * Title:          ExtractGaps.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <iostream>
#include <ostream>
#include <stdio.h>
#include <sstream>

#include "lav/LAVReader.h"
#include "lav/LAVFile.h"
#include "lav/LAVAlignedContig.h"
#include "DNASequence.h"
#include "SeqReader.h"
#include "align/alignutils.h"
#include "SeqUtils.h"
#include "utils.h"

#include "GapFunctions.h"

void InitEnv(int argc, char* argv[], 
	     std::string &refSeqName, 
	     std::string &qrySeqName,
	     std::string &lavFileName,
	     std::string &refGapSeqName,
	     std::string &qryGapSeqName,
	     ssize_t &minGapLen, 
	     ssize_t &strand,
	     double &muchGreater);

void PrintUsage();


int main(int argc, char* argv[]) {
  std::string qrySeqName, refSeqName, lavFileName;
  std::string qryGapFileName, refGapFileName;
  DNASequence refSeq, qrySeq, qrySeqRC;
  qrySeqName = "";
  refSeqName = "";
  lavFileName = "";
  ssize_t strand;
  LAVFile lavFile;
  ssize_t minGapLen;
  double muchGreater;
  std::ofstream refGapFile, qryGapFile;
  std::string refName, qryName;
  std::string regionName;
  strand      = -1;
  minGapLen   = 50;
  muchGreater = 5.0;
  InitEnv(argc, argv, refSeqName, qrySeqName, lavFileName, 
	  refGapFileName, qryGapFileName, minGapLen, strand, muchGreater);
  // Get input.

  //  openck(qryGapFileName, qryGapFile);
  openck(refGapFileName, refGapFile);

  ssize_t foundRefName, foundQryName, foundRegionName;
  foundRefName = GetSpeciesNameFromFile(refSeqName, refName);
  foundQryName = GetSpeciesNameFromFile(qrySeqName, qryName);
  foundRegionName = GetRegionNameFromFile(refSeqName, regionName);

  SeqReader::MaskRepeats();
  SeqReader::GetSeq(refSeqName, refSeq, SeqReader::noConvert);
  SeqReader::GetSeq(qrySeqName, qrySeq, SeqReader::noConvert);

  MakeRC(qrySeq, qrySeqRC);
  LAVReader::ReadLAVFile(lavFileName, lavFile);
  ssize_t a, b, cb;
  LAVAlignedContig *alignedContig;
  
  DNASequence refSubseq, qrySubseq;
  refSubseq._ascii = 1;
  qrySubseq._ascii = 1;
  ssize_t refGapLen, qryGapLen;
  ssize_t refGapNum, qryGapNum;
  refGapNum = 0;
  qryGapNum = 0;
  for (a = 0; a < lavFile.alignments.size(); a++) {
    alignedContig = lavFile.alignments[a];
    LAVBlock *block;
    if (strand == -1 || alignedContig->qryContig.strand == strand) {
      for (b = 0; b < alignedContig->alignments.size(); b++) {
	block = alignedContig->alignments[b];
	refSubseq.seq = &refSeq.seq[block->refBegin];
	refSubseq.length = block->refEnd - block->refBegin  + 1;
	
	if (alignedContig->qryContig.strand == 0)
	  qrySubseq.seq = &qrySeq.seq[block->qryBegin];
	else {
	  qrySubseq.seq = &qrySeqRC.seq[block->qryBegin];
	}
	for (cb = 1; cb < block->size(); cb++) {
	  refGapLen = block->refALBegin[cb] - block->refALEnd[cb-1];
	  qryGapLen = block->qryALBegin[cb] - block->qryALEnd[cb-1];
	  if (refGapLen > minGapLen and MuchGreater(refGapLen, qryGapLen, muchGreater)) {
	    if (foundRefName ) 
	      refGapFile << refName << " ";
	    if (foundQryName)
	      refGapFile << qryName << " ";
	    if (foundRegionName)
	      refGapFile << regionName << " ";
	    refGapFile << block->refALEnd[cb-1] << " " << block->refALBegin[cb] << " " 
		       << block->qryALBegin[cb-1] << " " << block->qryALEnd[cb] << " " 
		       << refGapLen+1 << " ";
	    refSubseq.seq = &refSeq.seq[block->refALEnd[cb-1]];
	    refSubseq.length = refGapLen+1;
	    refSubseq.PrintSeq(refGapFile,0);
	    refGapFile << std::endl;
	  }
	  /*
	    if (qryGapLen > minGapLen and MuchGreater(qryGapLen, refGapLen, muchGreater)) {
	    if (foundQryName)
	      qryGapFile << qryName << " ";
	    if (foundRefName ) 
	      qryGapFile << refName << " ";
	    if (foundRegionName)
	      qryGapFile << regionName << " ";
	    qryGapFile << block->qryALBegin[cb-1] << " " << block->qryALEnd[cb] << " " 
		       << block->refALEnd[cb-1] << " " << block->refALBegin[cb] << " " 
		       << qryGapLen + 1 << " ";
	    qrySubseq.seq = &qrySeq.seq[block->qryALEnd[cb-1]];
	    qrySubseq.length = qryGapLen+1;
	    qrySubseq.PrintSeq(qryGapFile,0);
	    qryGapFile << std::endl;
	  }
	  */
	}
      
	// Done looking inside the block
	if (b > 0) {
	  refGapLen = block->refBegin - alignedContig->alignments[b-1]->refEnd;
	  qryGapLen = block->qryBegin - alignedContig->alignments[b-1]->qryEnd;
	  if (refGapLen > minGapLen 
	      and MuchGreater(refGapLen, qryGapLen, muchGreater) 
	      and refGapLen < 1000 ) {
	    if (foundRefName ) 
	      refGapFile << refName << " ";
	    if (foundQryName)
	      refGapFile << qryName << " ";
	    if (foundRegionName)
	      refGapFile << regionName << " ";
	    refGapFile << block->refBegin << " " << alignedContig->alignments[b-1]->refEnd << " "
		       << block->qryBegin << " " << alignedContig->alignments[b-1]->qryEnd << " "
		       << refGapLen << " ";
	    refSubseq.seq = &refSeq.seq[block->refALEnd[cb-1]];
	    refSubseq.length = refGapLen+1;
	    refSubseq.PrintSeq(refGapFile,0);
	    refGapFile << std::endl;	
	  }
	  /*
	    if (qryGapLen > minGapLen 
	    and MuchGreater(qryGapLen, refGapLen, muchGreater)
	    and qryGapLen < 1000) {
	    if (foundQryName)
	    qryGapFile << qryName << " ";
	    if (foundRefName ) 
	    qryGapFile << refName << " ";
	    if (foundRegionName)
	    qryGapFile << regionName << " ";
	    qryGapFile << block->qryBegin << " " << alignedContig->alignments[b-1]->qryEnd << " "
	    << block->refBegin << " " << alignedContig->alignments[b-1]->refEnd << " " 
	    << qryGapLen << " ";
	    qrySubseq.seq = &qrySeq.seq[block->qryALEnd[cb-1]];
	    qrySubseq.length = qryGapLen+1;
	    qrySubseq.PrintSeq(qryGapFile, 0);
	    qryGapFile << std::endl;
	    }
	  */
	}
      }
    }
  }
  return 0;
}



void InitEnv(int argc, char* argv[], 
	     std::string &refSeqName, 
	     std::string &qrySeqName,
	     std::string &lavFileName,
	     std::string &refGapFileName,
	     std::string &qryGapFileName,
	     ssize_t &minGapLen,
	     ssize_t &strand, 
	     double &muchGreater) {

  ssize_t copt;
  ssize_t i;
  while ( (copt=getopt(argc, argv, "s:m:r:")) != EOF){
    switch(copt) {
    case 's':
      strand = atoi(optarg);
      continue;
    case 'm':
      minGapLen = atoi(optarg);
      continue;
    case 'r':
      muchGreater = atof(optarg);
      continue;
    default:
      PrintUsage();
      exit(1);
    }
  }
  i = optind;
  if (i >= argc) {
    std::cout << "You must specify a reference seq and query seq. " << std::endl;
    PrintUsage();
    exit(1);
  }
  refSeqName = argv[i];
  i++;
  if (i >= argc) {
    std::cout << "You must specify a query seq. " << std::endl;
    PrintUsage();
    exit(1);
  }
  qrySeqName = argv[i];
  i++;
  if (i >= argc) {
    std::cout << "You must specify a lav file. " << std::endl;
    PrintUsage();
    exit(1);
  }
  lavFileName = argv[i];
  i++;
  if (i >= argc) {
    std::cout << "You must specify a ref gap seq file. " << std::endl;
    PrintUsage();
    exit(1);
  }
  refGapFileName = argv[i];
  i++;
  /*
    if (i >= argc) {
    std::cout << "You must specify qry gap seq  file. " << std::endl;
    PrintUsage();
    exit(1);
    }
  qryGapFileName = argv[i];
  */
}


void PrintUsage() {
  std::cout << "extractGaps.  extract the gaps between two sequenes given a lav file. " << std::endl;
  std::cout << "usage: extractGaps [-s strand] [-m minGapLen] [-r ratio] refseq qryseq lavFile refgaps qrygaps " << std::endl;
  std::cout << "    refseq is the first indexed sequence in the lav alignment file " << std::endl;
  std::cout << "    qryseq is the second " << std::endl;
  std::cout << "    lavfile is the alignment file " << std::endl;
  std::cout << "    refgaps and qry gaps are the extracted gapped sequences " << std::endl;
  std::cout << "     -s strand  (0 = plus, 1 = minus  -1=both(default) " << std::endl;
  std::cout << "     -m minGapLen, only extract a gap if it is more than this len " << std::endl;
  std::cout << "     -r ratio, only extract a gap if it the ratio of lengths of the gaps is about this " << std::endl;
}




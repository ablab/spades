/***************************************************************************
 * Title:          lavExtractSeq.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
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


void InitEnv(int argc, char* argv[], 
	     std::string &refSeqName, 
	     std::string &qrySeqName,
	     std::string &lavFileName,
	     ssize_t &strand);

void PrintUsage();

int main(int argc, char* argv[]) {

  std::string qrySeqName, refSeqName, lavFileName;
  DNASequence refSeq, qrySeq, qrySeqRC;
  qrySeqName = "";
  refSeqName = "";
  lavFileName = "";
  ssize_t strand = -1;
  LAVFile lavFile;
  InitEnv(argc, argv, refSeqName, qrySeqName, lavFileName, strand);
  // Get input.

  SeqReader::MaskRepeats();
  SeqReader::GetSeq(refSeqName, refSeq);
  SeqReader::GetSeq(qrySeqName, qrySeq);
  
  MakeRC(qrySeq, qrySeqRC);
  LAVReader::ReadLAVFile(lavFileName, lavFile);
  ssize_t a, b;
  LAVAlignedContig *alignedContig;
  
  DNASequence refSubseq, qrySubseq;
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
    
	qrySubseq.length = block->qryEnd - block->qryBegin + 1;
	std::string refCtgName, qryCtgName;
	std::stringstream rcnstr, qcnstr;

	rcnstr.str(refCtgName);
	rcnstr << refSeqName << block->refBegin << "-" << block->refEnd;
	refSubseq.StoreName((char*)rcnstr.str().c_str());

	qcnstr.str(qryCtgName);
	qcnstr << qrySeqName << block->qryBegin << "-" << block->qryEnd;
	qrySubseq.StoreName((char*) qcnstr.str().c_str());
	refSubseq.PrintSeq(std::cout);
	std::cout << std::endl;
	qrySubseq.PrintSeq(std::cout);
	std::cout << std::endl;
      }
    }
  }
  return 0;
}



void InitEnv(int argc, char* argv[], 
	     std::string &refSeqName, 
	     std::string &qrySeqName,
	     std::string &lavFileName,
	     ssize_t &strand) {

  ssize_t copt;
  ssize_t i;
  while ( (copt=getopt(argc, argv, "s:")) != EOF){
    switch(copt) {
    case 's':
      strand = atoi(optarg);
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
}


void PrintUsage() {
  std::cout << "lavess.  extract sequences defined by lav blocks. " << std::endl;
  std::cout << "usage: laves refseq qryseq lavFile " << std::endl;
  std::cout << "     -s seq  (0 = plus, 1 = minus  -1=both(default) " << std::endl;
}




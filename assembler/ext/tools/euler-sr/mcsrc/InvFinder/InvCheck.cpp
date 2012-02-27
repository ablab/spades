/***************************************************************************
 * Title:          InvCheck.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

#include <vector>
#include <string>

#include "lav/LAVFile.h"
#include "lav/LAVReader.h"
#include "lav/LAVAlignedContig.h"
#include "lav/LAVBlock.h"

#include "SeqUtils.h"
#include "DNASequence.h"
#include "SeqReader.h"
#include "align/alignutils.h"
#include "lav/LAVUtils.h"

ssize_t InitEnv(int argc, char* argv[], 
	    std::string &selfAlignFileName,
	    std::vector <std::string> &alignmentNames,
	    std::string &seqDir,
	    ssize_t &diffThreshold);


void PrintUsage() {
  std::cout << "usage: invcheck [-t threshold] selfAlign  align1 align2 ... \n";
}


int main(int argc, char* argv[]) {

  std::string selfAlignFileName;
  std::vector< std::string> alignFileNames;

  // If a reverse complement alignment is this close to an alignment, 
  // consider it to be a palindrome.
  ssize_t proximityThreshold;
  std::string seqDir;
  ssize_t readRes;
  seqDir = "./";
  proximityThreshold = 200; 

  InitEnv(argc, argv, selfAlignFileName, alignFileNames, seqDir, proximityThreshold);

  LAVFile selfAlign;

  LAVReader::ReadLAVFile(selfAlignFileName, selfAlign);

  DNASequence selfSeq;

  readRes = SeqReader::GetSeq(seqDir + selfAlign.alignments[0]->refContig.sequenceName, 
			      selfSeq,
			      SeqReader::noConvert);

  if (readRes == 0) {
    std::cout << "could not find sequence " << selfAlign.alignments[0]->refContig.sequenceName
	      << " in dir: " << seqDir << std::endl;
    std::cout << "specify dir with \"-d\" option " << std::endl;
    exit(0);
  }
  std::vector<LAVFile*> alignments;

  ssize_t i;
  ssize_t a, ac, rcb;
  LAVAlignedContig *alignedContig;
  LAVBlock *block;
  std::vector<ssize_t> rcBegin, rcEnd;
  // Find out where all of the 
  for (ac = 0; ac < selfAlign.alignments.size(); ac++) {
    alignedContig = selfAlign.alignments[ac];
    if (alignedContig->qryContig.strand == 1) {
      for (a = 0; a < alignedContig->alignments.size(); a++ ) {
	block = alignedContig->alignments[a];
	rcBegin.push_back(block->refBegin);
	rcEnd.push_back(block->refEnd);
      }
    }
  }

  ssize_t alnFile; 
  LAVFile *lavFile, *qrySelfFile;
  ssize_t include;
  ssize_t qryLength;
  DNASequence qrySeq;
  DNASequence qryRC;
  std::vector<ssize_t> qrySelfBegin, qrySelfEnd;
  std::string qryName;
  alnFile = 0;
  while (alnFile < alignFileNames.size()-1) {
    lavFile = new LAVFile;
    qrySelfFile = new LAVFile;
    std::cout << alignFileNames[alnFile] << std::endl;

    LAVReader::ReadLAVFile(alignFileNames[alnFile], *lavFile);
    LAVReader::ReadLAVFile(alignFileNames[alnFile+1], *qrySelfFile);

    qrySelfBegin.clear();
    qrySelfEnd.clear();
    ssize_t orientation, reverse;
    ssize_t  beginPos, endPos;
    orientation = DetermineOrientation(lavFile);
    if (orientation == 0) 
      reverse = 1;
    else
      reverse = 0;
    ssize_t qryLength;
    // Build a list of reverse complement alignments of the query sequence
    for (ac = 0; ac < qrySelfFile->alignments.size(); ac++) {
      alignedContig = qrySelfFile->alignments[ac];
      qryLength = alignedContig->qryContig.end - alignedContig->qryContig.start + 1;
      // qry-self always is aligned with forward strand as 0 and reverse (where 
      // the palindromic would be) as 1.
      if (alignedContig->qryContig.strand == 1) {
	for (a = 0; a < alignedContig->alignments.size(); a++ ) {
	  block = alignedContig->alignments[a];
	  if (orientation == 1) {
	    qrySelfBegin.push_back(block->qryBegin);
	    qrySelfEnd.push_back(block->qryEnd);
	  }
	  else {
	    qrySelfBegin.push_back(qryLength - block->qryBegin + 1);
	    qrySelfEnd.push_back(qryLength - block->qryEnd + 1);
	  }	    
	}
      }
    }

     for (ac = 0; ac < lavFile->alignments.size(); ac++) {
      alignedContig = lavFile->alignments[ac];
      if (alignedContig->qryContig.strand == reverse) {
	qryLength = alignedContig->qryContig.end - alignedContig->qryContig.start + 1;
	qryName = alignedContig->qryContig.sequenceName;
	ssize_t endPos;
	if ((endPos = qryName.find(".fa")) >= 0) {
	  qryName = qryName.substr(0, endPos+3);
	}
	
	SeqReader::GetSeq(seqDir + qryName,
			  qrySeq, SeqReader::noConvert);
	MakeRC(qrySeq, qryRC);
      
	for (a = 0; a < alignedContig->alignments.size(); a++ ) {
	  block = alignedContig->alignments[a];
	  include = 1;
	  for (rcb = 0; rcb < rcBegin.size(); rcb++) {
	    // filter out reverse complements where there
	    // is a match against a reference-self or query-self alignment.
	    if (block->refBegin > (rcBegin[rcb]- proximityThreshold) && 
		block->refBegin < (rcEnd[rcb] + proximityThreshold)) {
	      include = 0;
	    }
	    if (block->refEnd > (rcBegin[rcb]- proximityThreshold) && 
		block->refEnd < (rcEnd[rcb] + proximityThreshold)) {
	      include = 0;
	    }
	  }
	  for (rcb = 0; rcb < qrySelfBegin.size(); rcb++) {
	    if (block->qryBegin > (qrySelfBegin[rcb]- proximityThreshold) && 
		block->qryBegin < (qrySelfEnd[rcb] + proximityThreshold)) {
	      include = 0;
	    }
	    if (block->qryEnd > (qrySelfBegin[rcb]- proximityThreshold) && 
		block->qryEnd < (qrySelfEnd[rcb] + proximityThreshold)) {
	      include = 0;
	    }
	  }
	  if (include) {
	    DNASequence fragment, fragmentRC;
	    fragment._ascii = 1;
	    fragment.seq = &qryRC.seq[block->qryBegin];
	    fragment.length = block->qryEnd - block->qryBegin + 1;

	    /*
	      fragment.PrintSeq(std::cout);
	      std::cout << std::endl;
	    */
	    /*
	    MakeRC(fragment, fragmentRC);
	    ssize_t *locations;
	    locations = new ssize_t[fragment.length];
	    double alignScore;
	    alignScore = AffineAlign(fragment, fragmentRC, 
				     -100, 200, 10000, 400, 30, 
				     locations, NULL);

	    if (alignScore < (-100 *  0.25 * fragment.length)) {
	      include = 0;
	      std::cout << "excluding alignment because it is palindromic " << std::endl;
	    }
	    delete[] locations;
	    */
	  }
	  if (include) {
	    std::cout << block->refBegin << " " << block->refEnd << " " 
		      << qryLength - block->qryBegin << " " 
		      << qryLength - block->qryEnd << std::endl;
	  }
	}
	delete[] qrySeq.seq;
	delete[] qryRC.seq;
	qrySeq.seq = NULL;
	qryRC.seq  = NULL;
      }
    }
    alnFile+=2;
    delete lavFile;
    delete qrySelfFile;
  }
  delete selfSeq.seq;
  return 0;
}


ssize_t InitEnv(int argc, char* argv[], 
	    std::string &selfAlignFileName,
	    std::vector <std::string> &alignmentNames,
	    std::string &seqDir,
	    ssize_t &diffThreshold) {
  ssize_t copt;
  while ( (copt=getopt(argc, argv, "t:d:")) != EOF){
    switch(copt) {
    case 't':
      diffThreshold = atoi(optarg);
      continue;
    case 'd':
      seqDir = optarg;
      continue;
    }
  }

  ssize_t i;
  i = optind;
  if (i >= argc) {
    PrintUsage();
    exit(0);
  }
  selfAlignFileName = argv[i];
  i++;
  if (i >= argc) {
    std::cout << "must specify at least one alignment \n";
    exit(0);
  }
    
  while (i < argc) {
    alignmentNames.push_back(argv[i]);
    i++;
  }
}



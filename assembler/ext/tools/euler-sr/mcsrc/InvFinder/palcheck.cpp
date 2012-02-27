/***************************************************************************
 * Title:          palcheck.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <ostream>
#include <stdio.h>
#include <unistd.h>
#include "DNASequence.h"
#include "SeqReader.h"
#include "align/alignutils.h"
#include "SeqUtils.h"
#include "EnumUtils.h"

void InitEnv(int argc, char* argv[], 
	     std::string &refSeqName, std::string &scoreMatFile,
	     ssize_t &gapOpen, ssize_t &gapExtend, 
	     ssize_t &printAlignment, ssize_t &localAlignment, ssize_t &alignConsecutive);

void PrintUsage();

int main(int argc, char* argv[]) {
  std::string seqName, scoreMatFile;
  DNASequence seq, seqRC, prev;
  FloatMatrix scoreMat;
  seqName ="";
  scoreMatFile = "../align/data/scoremat.txt";
  ssize_t gapOpen, gapExtend;
  ssize_t printAlignment;
  ssize_t alignConsecutive;
  ssize_t localAlign;

  gapOpen = 400;
  gapExtend = 30;
  printAlignment = 0;
  alignConsecutive = 0;
  localAlign = 0;
  InitEnv(argc, argv, seqName, scoreMatFile, gapOpen, gapExtend, 
	  printAlignment, localAlign, alignConsecutive);
  // Get input.
  std::ifstream in;
  in.open(seqName.c_str());
  if (!in.good()) {
    std::cout << "could not open " << seqName << std::endl;
    exit(0);
  }
  ssize_t i;
  if (scoreMatFile != "") {
    ReadScoreMatFile(scoreMatFile, scoreMat);
  }
  SeqReader seqReader(&in);

  prev.seq = NULL;
  ssize_t seqNum = 0;
  while (seqReader.GetSeq(seq) && seq.length > 0) {
    MakeRC(seq, seqRC);
    ssize_t *alignment;
    
    alignment = new ssize_t[seq.length];
    double alignScore;
    if (localAlign) {
      alignScore = AffineLocalAlign(seq, seqRC,
				    -100, // match
				    200,  // general mismatch
				    800,  // general gap
				    gapOpen,  // general gap-open
				    gapExtend,   // general gap-extend
				    alignment,  // resulting alignment
				    scoreMat // custom scoring matrix
				    );
    }
    else {
      alignScore = AffineAlign(seq, seqRC,
			       -100, // match
			       200,  // general mismatch
			       800,  // general gap
			       gapOpen,  // general gap-open
			       gapExtend,   // general gap-extend
			       alignment,  // resulting alignment
			       scoreMat // custom scoring matrix
			       );
    }
      
	std::cout << seqNum++<< "\t" <<  seq.length << "\t" << alignScore << "\t" <<  seq.namestr << std::endl;
    T_Alignment tmpAlignment;
    ssize_t length;
    
    if (printAlignment) {
      length = seq.length;
      ssize_t e = 0;
      ssize_t *refLocations = new ssize_t[length];
      ssize_t *qryLocations = new ssize_t[length];
      ssize_t *enumerations = new ssize_t[length];
      for (i = 0; i < length; i++) {
	if (alignment[i] != -1) {
	  enumerations[e] = e+1;
	  refLocations[e] = i;
	  qryLocations[e] = alignment[i];
	  e++;
	}
      }
      PrintAlignment(&std::cout, seq, seqRC, enumerations, refLocations, qryLocations, e);

      //      delete[] enumerations;
      //      delete[] refLocations;
      //      delete[] qryLocations;
    }
    
    //    delete[] alignment;

    if (alignConsecutive) {
      if (prev.seq == NULL) {
	  prev.seq = seq.seq;
	  prev.length = seq.length;
      }
      else {
	alignment = new ssize_t[prev.length];
	alignScore = Align(prev, seq, -100, 200, 400, alignment, scoreMat);
	std::cout << prev.length << " " << seq.length << "\t" 
		  << alignScore << "\t" <<  seq.namestr << std::endl;
	prev.seq = NULL;
	
	if (printAlignment) {
	  length = seq.length;
	  ssize_t e = 0;
	  ssize_t *refLocations = new ssize_t[length];
	  ssize_t *qryLocations = new ssize_t[length];
	  ssize_t *enumerations = new ssize_t[length];
	  for (i = 0; i < length; i++) {
	    if (alignment[i] != -1) {
	      enumerations[e] = e+1;
	      refLocations[e] = i;
	      qryLocations[e] = alignment[i];
	      e++;
	    }
	  }
	  PrintAlignment(&std::cout, seq, seqRC, enumerations, refLocations, qryLocations, e);
	    
	  //	  delete[] enumerations;
	  //	  delete[] refLocations;
	  //	  delete[] qryLocations;
	}
	//	delete[] alignment;
      }
    }    
  }


  return 0;
}



void InitEnv(int argc, char* argv[], 
	     std::string &seqName,
	     std::string &scoreMatFile,
	     ssize_t &gapOpen, ssize_t &gapExtend,
	     ssize_t &printAlignment, 
	     ssize_t &localAlignment,
	     ssize_t &alignConsecutive) {

  ssize_t copt;
  ssize_t i;
  while ( (copt=getopt(argc, argv, "m:g:e:acl")) != EOF){
    switch(copt) {
    case 'a':
      printAlignment = 1;
      continue;
    case 'c':
      alignConsecutive = 1;
      continue;
    case 'm':
      scoreMatFile = optarg;
      continue;
    case 'g':
      gapOpen = atoi(optarg);
      continue;
    case 'e':
      gapExtend = atoi(optarg);
      continue;
    case 'l':
      localAlignment = 1;
      continue;
    default:
      PrintUsage();
      exit(1);
    }
  }
  i = optind;
  if (i >= argc) {
    std::cout << "You must specify a seq. " << std::endl;
    PrintUsage();
    exit(1);
  }
  seqName = argv[i];
}


void PrintUsage() {
  std::cout << "palcheck. Find  palindromes in sequences. " << std::endl;
  std::cout << "usage:  palcheck [-m scoremat] [-g gapOpen] [-e gapExtend] [-a] seq.fasta " << std::endl;
  std::cout << "      -a will print the alignment " << std::endl;
  std::cout << "      -l performs local alignment " << std::endl;
}

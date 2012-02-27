/***************************************************************************
 * Title:          aalign.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <sstream>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <unistd.h>

#include "mctypes.h"
#include "DNASequence.h"
#include "SeqReader.h"
#include "AlignmentPrinter.h"
#include "align/alignutils.h"


void PrintUsage();
void initenv(int argc, char *argv[], 
	     std::string &inputFileA,
	     std::string &inputFileB,
	     ssize_t &mismatch,
	     ssize_t &gapOpen,
	     ssize_t &gapExtend,
	     ssize_t &nonAffineGap,
	     ssize_t &match);

int main(int argc, char* argv[]) {
  
  DNASequence *seq1;
  DNASequence *seq2;
  
  ssize_t mismatch, gapOpen, gapExtend, nonAffineGap, match;
  match        = -1;
  mismatch     = 1;
  nonAffineGap = 2;
  gapOpen      = 6;
  gapExtend    = 0;
  
  std::string inA, inB;
  initenv(argc, argv, inA, inB, mismatch, gapOpen, gapExtend, nonAffineGap, match);

  std::ifstream ifa, ifb;

  ifa.open(inA.c_str());
  if (! ifa.good() ) { std::cout << "could not open " << inA << std::endl; exit(0);}
  
  ifb.open(inB.c_str());
  if (! ifb.good() ) { std::cout << "could not open " << inB << std::endl; exit(0);}
  
  SeqReader readerA(&ifa);
  SeqReader readerB(&ifb);

  readerA.GetSeq(seq1, SeqReader::noConvert);
  readerB.GetSeq(seq2, SeqReader::noConvert);

  ifa.close(); ifb.close();

  if (seq1 == NULL || seq2 == NULL) {
    std::cout << "Error getting sequences from " << argv[1] << std::endl;
  }
  ssize_t *qryLocations, *refLocations, *enumerations, *locations;
  locations    = new ssize_t[seq1->length];
  qryLocations = new ssize_t[seq1->length];
  refLocations = new ssize_t[seq1->length];
  enumerations  = new ssize_t[seq1->length];
  ssize_t i;
  for (i = 0; i < seq1->length; i++) {
    qryLocations[i] = refLocations[i] = -1;
    enumerations[i] = 0;
  }
  double score;
  score = AffineAlign(*seq1, *seq2, 
		      match, mismatch, nonAffineGap, gapOpen, gapExtend, 
		      locations);
  
	PrintAlignment(*seq1, *seq2, 0, 0, locations, seq1->length + 1, std::cout);

  std::cout << "score: " << score << std::endl;
  
}

void initenv(int argc, char *argv[], 
	     std::string &inputFileA,
	     std::string &inputFileB,
	     ssize_t &mismatch, 
	     ssize_t &gapOpen,
	     ssize_t &gapExtend,
	     ssize_t &nonAffineGap,
	     ssize_t &match) {
  ssize_t copt;
  while ( (copt=getopt(argc, argv, "i:j:m:M:o:g:e:h")) != EOF) {
    switch(copt) {
    case 'i':
      inputFileA = optarg;
      continue;
    case 'j':
      inputFileB = optarg;
      continue;
    case 'm':
      match = atoi(optarg);
      continue;
    case 'M':
      mismatch = atoi(optarg);
      continue;
    case 'o':
      gapOpen = atoi(optarg);
      continue;
    case 'e':
      gapExtend = atoi(optarg);
      continue;
    case 'g':
      nonAffineGap = atoi(optarg);
      continue;
    case 'h':
      PrintUsage();
      exit(0);
    default:
      std::cout << "error with option: " << copt << " arg: " << optarg << std::endl;
      PrintUsage();
      exit(1);
    }
  }
}

void PrintUsage() {
  std::cout << "Affine align two sequences from a file\n" << std::endl;
  std::cout << "usage:  aalign  file -m match -M mismatch "
	    << "-o gapOpen -e gapExtend -g nonAffineGap " << std::endl;
  std::cout << "Alignment corresponds to minimum penalty path (lowest scoring). " << std::endl;
  std::cout << "   -i inputfile input seq" << std::endl;
  std::cout << "   -j query seq " << std::endl;
  std::cout << "   -m match (-1.0) nucleotide match penalty " << std::endl;
  std::cout << "   -M mismatch (1.0) nucleotide mismatch penalty " << std::endl;
  std::cout << "   -o gapOpen (6.0) penalty for opening a gap. " << std::endl;
  std::cout << "   -e gapExtend (0.0) penalty for extending a gap. " << std::endl;
  std::cout << "   -g gap (1.0) penalty for nonaffine gap. " << std::endl;
}

/***************************************************************************
 * Title:          endalign.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  09/07/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <stdlib.h>
#include <string>
#include <unistd.h>

#include "DNASequence.h"
#include "SeqReader.h"
//#include "EnumUtils.h"
#include "AlignmentPrinter.h"
#include "align/alignutils.h"


void PrintUsage();
void InitEnv(int argc, char *argv[], 
	     std::string &inputFileA,
	     std::string &inputFileB,
	     ssize_t &mismatch,
	     ssize_t &gapOpen,
	     ssize_t &gapExtend,
	     ssize_t &nonAffineGap,
	     ssize_t &match,
	     ssize_t &kband,
	     std::string &alignType);

ssize_t EndAlign2Seq(DNASequence &longSeq, DNASequence &shortSeq,
		  ssize_t *&alignmentLocations);
int main(int argc, char* argv[]) {
  
  DNASequence *seq1;
  DNASequence *seq2;
  ssize_t kband;
  std::string alignType;
  ssize_t mismatch, gapOpen, gapExtend, nonAffineGap, match;
  match        = -1;
  mismatch     = 1;
  nonAffineGap = 2;
  gapOpen      = 6;
  gapExtend    = 0;
  kband        = 4;
  ssize_t score;

  std::string inA, inB;
  alignType = "affine";
  inA = "";
  inB = "";

  InitEnv(argc, argv, inA, inB, 
	  mismatch, gapOpen, gapExtend, nonAffineGap, match, kband, alignType);

  if (inA == "" || inB == "") {
    PrintUsage();
    exit(1);
  }
  std::ifstream ifa, ifb;

  ifa.open(inA.c_str());
  if (! ifa.good() ) { std::cout << "could not open " << inA << std::endl; exit(0);}
  
  ifb.open(inB.c_str());
  if (! ifb.good() ) { std::cout << "could not open " << inB << std::endl; exit(0);}
  
  SeqReader readerA(&ifa);
  SeqReader readerB(&ifb);

  readerA.GetSeq(seq1);//readerA.GetSeq(seq1, SeqReader::noConvert);
  readerB.GetSeq(seq2);//readerB.GetSeq(seq2, SeqReader::noConvert);
  

  ifa.close(); ifb.close();

  if (seq1 == NULL || seq2 == NULL) {
    std::cout << "Error getting sequences from " << argv[1] << std::endl;
  }
  ssize_t *qryLocations, *refLocations, *enumerations, *locations;
  ssize_t *optScores;
  locations    = new ssize_t[seq1->length];
  optScores    = new ssize_t[seq1->length];
  qryLocations = new ssize_t[seq1->length];
  refLocations = new ssize_t[seq1->length];
  enumerations  = new ssize_t[seq1->length];
  ssize_t i;
  for (i = 0; i < seq1->length; i++) {
    qryLocations[i] = refLocations[i] = -1;
    enumerations[i] = 0;
  }

  for (i = 0; i < seq1->length; i++) {
    locations[i] = -1;
    optScores[i] = INF;
  }
  score = EndAlign2Seq(*seq1, *seq2, locations);

  /*
    for (i = 0; i < seq1->length; i++) 
    std::cout << std::setw(4) << locations[i]; 
    std::cout << std::endl;
  */
  ssize_t e;
  e = 0;
  for (i = 0; i < seq1->length; i++) {
    if (locations[i] != -1) {
      enumerations[e] = e+1;
      refLocations[e] = i;
      qryLocations[e] = locations[i];
      e++;
    }
  }
	PrintAlignment(*seq1, *seq2, 0, 0, locations, seq1->length + 1, std::cout);
}

void InitEnv(int argc, char *argv[], 
	     std::string &inputFileA,
	     std::string &inputFileB,
	     ssize_t &mismatch, 
	     ssize_t &gapOpen,
	     ssize_t &gapExtend,
	     ssize_t &nonAffineGap,
	     ssize_t &match,
	     ssize_t &kband,
	     std::string &alignType) {
  ssize_t copt;
  while ( (copt=getopt(argc, argv, "i:j:m:M:o:g:e:ht:k:")) != EOF) {
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
    case 'k':
      kband = atoi(optarg);
      continue;
    case 't':
      alignType = optarg;
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
  std::cout << "usage:  aalign  file -m match -M mismatch -o gapOpen -e gapExtend -g nonAffineGap " << std::endl;
  std::cout << "Alignment corresponds to minimum penalty path (lowest scoring). " << std::endl;
  std::cout << "   -i inputfile input seq" << std::endl;
  std::cout << "   -j query seq " << std::endl;
  std::cout << "   -t alignType , where aligntype can be : "<< std::endl;
  std::cout << "                  affine - for affine alignment " << std::endl;
  std::cout << "                  band   - for banded alignment " << std::endl;
  std::cout << "                  fit    - for fitting alignment where full query sequence is fit into " << std::endl
	    << "                           reference sequence " << std::endl;
  std::cout << "   -m match (-1.0) nucleotide match penalty " << std::endl;
  std::cout << "   -M mismatch (1.0) nucleotide mismatch penalty " << std::endl;
  std::cout << "   -o gapOpen (6.0) penalty for opening a gap. " << std::endl;
  std::cout << "   -e gapExtend (0.0) penalty for extending a gap. " << std::endl;
  std::cout << "   -g gap (1.0) penalty for nonaffine gap. " << std::endl;
}


ssize_t EndAlign2Seq(DNASequence &longSeq, DNASequence &shortSeq,
		  ssize_t *&alignmentLocations) {
  // End alignment means finding the optimal alignment of the shorter
  // sequence in the longer sequence such that the ends of the sequences
  // align together, and there is one large gap allowed in the 
  // long sequence.

  ssize_t band = 10;
  ssize_t *forLocations, *revLocations;
  ssize_t *forScores /* and seven years ago */, *revScores;

  forLocations = new ssize_t[shortSeq.length + band];
  revLocations = new ssize_t[shortSeq.length + band];
  forScores  = new ssize_t[shortSeq.length + band];
  revScores  = new ssize_t[shortSeq.length + band];
  
  ssize_t i, j;
  for (i = 0; i < shortSeq.length + band; i++) {
    forLocations[i] = revLocations[i] = -1;
    forScores[i] = revScores[i] = INF;
  }
  // Align the short sequence in the longer one
  BandedAlign(longSeq, shortSeq,
							-4, 5, 5,
							band,
							forLocations, forScores);

  // Check out the result of the alignment
  //UNUSED// ssize_t *tmpLoc = new ssize_t[shortSeq.length+band];
  //UNUSED// ssize_t *tmpEnum = new ssize_t[shortSeq.length+band];
  //UNUSED// ssize_t *tmpRefEnum = new ssize_t[shortSeq.length+band];
  //UNUSED// ssize_t *tmpQryEnum = new ssize_t[shortSeq.length];
  DNASequence tmpLongSeq;
  tmpLongSeq.seq = longSeq.seq;
  tmpLongSeq.length = shortSeq.length + band;
	/*  int enumLength = 0;
  enumLength = AlignmentToEnumeration(tmpLongSeq, 
				      forLocations,
				      1,
				      tmpEnum, tmpRefEnum, tmpQryEnum);
	*/
  // Align in the reverse direction, starting from the ends of 
  // each string.
  DNASequence shortRev(shortSeq.length);
  DNASequence longRev(shortSeq.length+band);

  // Find the reverse of each sequence.  Note this is not 
  // the reverse complement, just the reverse (although the reverse
  // complement would do just same here with a little extra work).
  shortRev.StoreName("shortseq rev");
  shortRev._ascii= 1;
  longRev._ascii = 1;
  for (i = 0; i < shortSeq.length; i++)
    shortRev.seq[i] = shortSeq.seq[shortSeq.length-1-i]; 
  
  for (i = 0; i < shortSeq.length + band; i++) 
    //    longRev.seq[i] = longSeq.seq[longSeq.length - shortSeq.length - band - 1 + i];
    longRev.seq[i] = longSeq.seq[longSeq.length - 1 - i];

  BandedAlign(longRev, shortRev,
	      -4, 5, 5,
	      band,
	      revLocations, revScores);

  // Now find the optimal combination of the two sequences
  ssize_t minScore = INF;
  //UNUSED// ssize_t minIndex = 0;
  //UNUSED+// ssize_t  fordx;
  ssize_t revdx;
  ssize_t minFor, minRev;
  minFor = minRev = 0;
  for (i = 0; i < shortSeq.length+band-2; i++) {
    for (j = 0; j < shortSeq.length+band-2-i-1; j++) {
      revdx = shortSeq.length+band-1-j-1;
      if (forScores[i] + revScores[revdx] < minScore) {
	minScore = forScores[i] + revScores[revdx];
	minFor = i;
	minRev = j;
      }
    }
  }
	/*  enumLength = AlignmentToEnumeration(longRev, 
				      revLocations,
				      -1,
				      tmpEnum, tmpRefEnum, tmpQryEnum);
	*/
  // Now forLocations map from the short (query) sequence, and 
  // revLocations map from the reverse short sequence into the end of the 
  // reverse long sequence.  
  for (i = 0; i < longSeq.length; i++) {
    alignmentLocations[i] = -1;
  }

  for (i = 0; i <= minFor ; i++)
    alignmentLocations[i] = forLocations[i];

  for (i = 0; i < minRev; i++) 
    if (revLocations[i] != -1)
      alignmentLocations[longSeq.length-1-i] = shortSeq.length -1 - revLocations[i];

  return minScore;
}

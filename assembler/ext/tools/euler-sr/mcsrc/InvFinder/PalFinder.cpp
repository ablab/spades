/***************************************************************************
 * Title:          PalFinder.cpp 
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

#include "DNASequence.h"
#include "SeqReader.h"
#include "align/alignutils.h"
#include "StripGen.h"
#include "TupleLib.h"
#include "SeqUtils.h"

#include "SimpleStats.h"

double AlignAroundPosition(DNASequence &refSeq, ssize_t refStartPos, ssize_t refEndPos,
			  DNASequence &qrySeq, ssize_t qryStartPos, ssize_t qryEndPos,
			  ssize_t *&locations, // map of refseq to qryseq
			  ssize_t &length, // length of optimal local alignment
			  ssize_t &refStartAlignPos, ssize_t &qryStartAlignPos,
			  double **scoreMat);

double AlignSubsequences(DNASequence &refSeq, ssize_t refPos, ssize_t refLength,
			DNASequence &qrySeq, ssize_t qryPos, ssize_t qryLength, 
			ssize_t *&locations, // map of refseq to qryseq
			ssize_t &length, // length of optimal local alignment
			ssize_t &refAlignStart, // start of local alignment in ref seq
			ssize_t &qryAlignStart, // start of local alignment in qry seq
			double *scoreMat[5],
			ssize_t  extendAlignment);

void InitEnv(int argc, char* argv[], 
	     std::string &refSeqName, 
	     ssize_t &startPos, ssize_t &endPos,
	     std::string &scoreMatFile,
	     std::string &scoreMatString,
	     double &gapOpen, double &gapExtend,
	     double &initialScore,
	     double &tailEndPct);


void PrintUsage();

int main(int argc, char* argv[]) {

  std::string refSeqName;
  std::string scoreMatFile, scoreMatString;

  DNASequence refSeq, refRC;

  FloatMatrix scoreMat;
  ssize_t i;
  double topPct;
  double gapOpen, gapExtend, initialScore;
  ssize_t startPos, endPos;
  topPct = 0.05;
  scoreMatFile = "";
  scoreMatString = "";

  gapOpen      = 400;
  gapExtend    = 30;
  initialScore = -1500;
  startPos = -1;
  endPos   = -1;
  InitEnv(argc, argv, 
	  refSeqName, 
	  startPos, endPos,
	  scoreMatFile, scoreMatString, gapOpen, gapExtend, initialScore, 
	  topPct);

  /*
  std::ofstream output;
  std::ostream *outputPtr;
  if (outputName == "")
    outputPtr = &std::cout;
  else {
    output.open(outputName.c_str());
    if (!output.good()) {
      std::cout << "Could not open output file " << outputName << std::endl;
      exit(0);
    }
    outputPtr = &output;
  }
  */
  if (scoreMatFile != "") {
    ReadScoreMatFile(scoreMatFile, scoreMat);
  }

  // Get input.
  SeqReader::MaskRepeats();
  SeqReader::GetSeq(refSeqName, refSeq, SeqReader::noConvert);

  if (startPos == -1)
    startPos = 0;
  if (endPos  == -1)
    endPos = refSeq.length;

  // Reverse, but do not complement backwards sequence
  MakeRC(refSeq, refRC);

  
  double palScore;
  ssize_t *palLocations, palLength;
  ssize_t refAlignStart, qryAlignStart;
  DNASequence tmpseq;
  tmpseq._ascii= 1;

  for (i = startPos; i < endPos; i++) {
    palScore= ScoreBandedAffineAlign(refSeq, refRC, // sequences to compare 
				     initialScore,  // score of seed
				     scoreMat, gapOpen, gapExtend, // scoring parameters
				     -2*initialScore,  // negative score is good.
				     palLocations, palLength,
				     refAlignStart, qryAlignStart,
				     i, refSeq.length - i + 1);
    if (palLength > 25) {
      assert(i-palLength > 0);
      tmpseq.seq = &refSeq.seq[i - palLength];
      tmpseq.length = 2*palLength;
      if (CountRepeatMasked(tmpseq) / double(tmpseq.length) < 0.25) {
	std::cout << i << " palindrome length " << palLength << std::endl;
	tmpseq.PrintSeq(std::cout);
	std::cout << std::endl;
      }
      i += palLength;
    }
    delete[] palLocations;
    if (i % 1000 == 0) 
      std::cout << "pos: " << i << std::endl;
  }
  return 0;
}


void InitEnv(int argc, char* argv[], 
	     std::string &refSeqName, 
	     ssize_t &startPos, ssize_t &endPos,
	     std::string &scoreMatFile,
	     std::string &scoreMatString,
	     double &gapOpen, double &gapExtend,
	     double &initialScore,
	     double &tailEndPct) {

  ssize_t copt;
  ssize_t i;
  std::string inpfile;
  while ( (copt=getopt(argc, argv, "g:e:i:m:S:o:a:t:s:f:")) != EOF){
    switch(copt) {
    case 'i':
      initialScore = atof(optarg);
      continue;
    case 'g':
      gapOpen = atof(optarg);
      continue;
    case 's':
      startPos = atoi(optarg);
      continue;
    case 'f':
      endPos  = atoi(optarg);
      continue;
    case 'e':
      gapExtend = atof(optarg);
      continue;
      /*    case 'o':
      outputFile = optarg;
      continue;
      */
    case 't':
      tailEndPct = atof(optarg);
      continue;
    case 'm':
      scoreMatFile = optarg;
      continue;
    case 'S':
      scoreMatString = optarg;
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
}


void PrintUsage() {
  std::cout << "palfinder.  A program to find palindromes in a sequence. " << std::endl;
  std::cout << "usage: invfinder [-otmswe] refseq  " << std::endl;
  std::cout << "-m scoreMatFile : use scoring matrix scoremat " << std::endl;
  std::cout << "-s scoreMatStr  : use scoring matrix string that specifies the 25 values of scoremat" << std::endl;
  std::cout << "-t tailEndPct   : output top t percent of alignment lengths " << std::endl;
}


double AlignAroundPosition(DNASequence &refSeq, ssize_t refStartPos, ssize_t refEndPos,
			  DNASequence &qrySeq, ssize_t qryStartPos, ssize_t qryEndPos,
			  ssize_t *&locations, // map of refseq to qryseq
			  ssize_t &length, // length of optimal local alignment
			  ssize_t &refStartAlignPos, ssize_t &qryStartAlignPos,
			  FloatMatrix &scoreMat) {
  std::cout << "aligning around position: " << refStartPos << std::endl; 
  DNASequence qrySubseqRC;

  MakeRC(qrySeq, qrySubseqRC);

  
  ssize_t pos, i;
  ssize_t *forAlignment, *backAlignment;
  ssize_t forLength, backLength;
  DNASequence refExt, qryExtRC;
  ssize_t startRefAlign, endRefAlign, startQryAlign, endQryAlign;
  pos = 0;
  double alignScore;
  ssize_t forRefAlignStart, forQryAlignStart, backRefAlignStart, backQryAlignStart;
  alignScore = (refEndPos - refStartPos+1)*scoreMat[0][0];

  backAlignment = new ssize_t[refStartPos];
  for (i = 0; i < refStartPos; i++) backAlignment[i] = -1;

  
  ScoreBandedAffineAlign(refSeq, qrySubseqRC, alignScore, scoreMat, 
			 400,  30, alignScore -1500, backAlignment, backLength,
			 backRefAlignStart, backQryAlignStart,
			 refStartPos-1, qryStartPos-1, -1);
  
  ScoreBandedAffineAlign(refSeq, qrySubseqRC, alignScore, scoreMat, 
			 400,  30, alignScore-1500, forAlignment, forLength,
			 forRefAlignStart, forQryAlignStart,
			 refEndPos+1, qryEndPos+1, 1);

  std::cout << "aligning ref seq " << std::endl;
  refSeq.PrintSeq(std::cout);
  std::cout << "qry rc: " << std::endl;
  qrySubseqRC.PrintSeq(std::cout);
  std::cout << " extended alignemnt of score : " << alignScore << " " << forLength << std::endl;
    
  PrintAlignment(&std::cout, refSeq, qrySubseqRC, forAlignment , 1, forLength);

  // Now join the two alignments into one larger one.

  length = backLength + (refEndPos - refStartPos+1) + forLength;

  locations = new ssize_t[length];
  pos = 0;
  // copy the back alignment
  for(i =0; i < backLength; i++, pos++) {
    locations[pos] = backAlignment[i];
  }
  // copy the seed
  for(i = 0; i < refEndPos - refStartPos+1; i++, pos++) {
    locations[pos] = qryStartPos + i;
  }
  // copy the forward alignment
  for (i = 0; i < forLength; i++, pos++) {
    locations[pos] = forAlignment[i];
  }

}

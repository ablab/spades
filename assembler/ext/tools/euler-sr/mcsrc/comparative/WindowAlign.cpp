/***************************************************************************
 * Title:          WindowAlign.cpp 
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
#include "alignutils.h"
#include "SeqUtils.h"
#include "utils.h"
#include "EnumUtils.h"
#include "mctypes.h"

void InitEnv(int argc, char* argv[], 
	     std::string &refSeqName, 
	     std::string &qrySeqName,
	     std::string &coordsFile,
	     std::string &scoreMatFileName,
	     ssize_t &extLen);

void PrintUsage();
ssize_t ReadCoordsFile(std::string coordsFileName, 
		   std::vector<ssize_t> &refStartCoords,
		   std::vector<ssize_t> &refEndCoords,
		   std::vector<char> &refStrands,
		   std::vector<ssize_t> &qryStartCoords,		   
		   std::vector<ssize_t> &qryEndCoords,
		   std::vector<char> &qryStrands);

int main(int argc, char* argv[]) {

  std::string qrySeqName, refSeqName, coordsFileName, scoreMatFileName;
  DNASequence refSeq, qrySeq, refRC, qryRC;
  ssize_t extLen;
  FloatMatrix scoreMat;

  ssize_t i, j;

  qrySeqName = "";
  refSeqName = "";
  scoreMatFileName = "";
  InitEnv(argc, argv, refSeqName, qrySeqName, coordsFileName, scoreMatFileName, extLen);
  // Get input.
  SeqReader::GetSeq(refSeqName, refSeq, SeqReader::noConvert);
  SeqReader::GetSeq(qrySeqName, qrySeq, SeqReader::noConvert);
  
  MakeRC(refSeq, refRC);
  MakeRC(qrySeq, qryRC);

  if (scoreMatFileName != "") 
    ReadScoreMatFile(scoreMatFileName, scoreMat);
  //  else 
  //    InitScoreMat(scoreMatPtr, -100, 300);

  std::vector<ssize_t> refStartCoords;
  std::vector<ssize_t> refEndCoords;
  std::vector<char> refStrands;
  std::vector<ssize_t> qryStartCoords;
  std::vector<ssize_t> qryEndCoords;
  std::vector<char> qryStrands; 
  ReadCoordsFile(coordsFileName,   
		 refStartCoords,
		 refEndCoords,
		 refStrands,
		 qryStartCoords,
		 qryEndCoords,
		 qryStrands);

  ssize_t numCoordPairs = refStartCoords.size();
  DNASequence *refSeqPtr, *qrySeqPtr;
  DNASequence refWinSeq, qryWinSeq;

  double windowScore, upstrScore, downstrScore;
  for (i = 0; i < numCoordPairs; i++) {
    if (refStrands[i] == 1) 
      refSeqPtr = &refSeq;
    else
      refSeqPtr = &refRC;

    if (qryStrands[i] == 1)
      qrySeqPtr = &qrySeq;
    else
      qrySeqPtr = &qryRC;

    // Align inside the window
    refWinSeq.seq = &(refSeqPtr->seq[refStartCoords[i]]);
    refWinSeq.length = refEndCoords[i] - refStartCoords[i] + 1;
    
    qryWinSeq.seq = &(qrySeqPtr->seq[qryStartCoords[i]]);
    qryWinSeq.length = qryEndCoords[i] - qryStartCoords[i] + 1;

    ssize_t *locations, *upstrLocations, *downstrLocations;
    ssize_t upstrLength, downstrLength;
    ssize_t refAlignStart, qryAlignStart; // for now throw away refalign & qryalign start
    locations = new ssize_t[refWinSeq.length];
    for (j = 0; j < refWinSeq.length; j++)
      locations[j] = 0;

    windowScore = AffineAlign(refWinSeq, qryWinSeq, 0,0, 0, 400, 30, 
			      locations, scoreMat);

    upstrScore  = ScoreBandedAffineAlign(*refSeqPtr, *qrySeqPtr, 
					 windowScore,
					 scoreMat, 400, 30, 
					 -windowScore,
					 upstrLocations, upstrLength,
					 refAlignStart, qryAlignStart,
					 refStartCoords[i], qryStartCoords[i], 
					 -1);
    
    downstrScore = ScoreBandedAffineAlign(*refSeqPtr, *qrySeqPtr, 
					  windowScore,
					  scoreMat, 400, 30, 
					  -windowScore,
					  downstrLocations, downstrLength,
					  refAlignStart, qryAlignStart,
					  refEndCoords[i] , qryEndCoords[i]);

    delete[] locations;
    delete[] upstrLocations;
    delete[] downstrLocations;
      
    std::cout << "alignScore: " << windowScore << " upstrscore: " << upstrScore << " uplen: " 
	      << upstrLength 
	      << " downstrscore: " << downstrScore << " downlength: " << downstrLength << std::endl;
  }
  
  return 0;
}

ssize_t ReadCoordsFile(std::string coordsFileName, 
		   std::vector<ssize_t> &refStartCoords,
		   std::vector<ssize_t> &refEndCoords,
		   std::vector<char> &refStrands,
		   std::vector<ssize_t> &qryStartCoords,		   
		   std::vector<ssize_t> &qryEndCoords,
		   std::vector<char> &qryStrands) {
  std::ifstream cf;
  openck(coordsFileName, cf);

  ssize_t refStartPos, refEndPos;
  ssize_t qryStartPos, qryEndPos, refStrand, qryStrand;
  while (cf) {
    if (! (cf >> refStartPos >> refEndPos >> refStrand 
	   >> qryStartPos >> qryEndPos >> qryStrand ))
      break;
    refStartCoords.push_back(refStartPos);
    refEndCoords.push_back(refEndPos);
    refStrands.push_back(refStrand);
    
    qryStartCoords.push_back(qryStartPos);
    qryEndCoords.push_back(qryEndPos);
    qryStrands.push_back(qryStrand);
  }
}




void InitEnv(int argc, char* argv[], 
	     std::string &refSeqName, 
	     std::string &qrySeqName,
	     std::string &coordsFile,
	     std::string &scoreMatFileName,
	     ssize_t &extLen) {
  ssize_t copt;
  ssize_t i;
  while ( (copt=getopt(argc, argv, "m:")) != EOF){
    switch(copt) {
    case 'm':
      scoreMatFileName = optind;
      continue;
    default:
      PrintUsage();
      exit(1);
    }
  }
  i = optind;
  if (i >= argc) {
    PrintUsage();
    exit(1);
  }
  refSeqName = argv[i];
  i++;
  if (i >= argc) {
    PrintUsage();
    exit(1);
  }
  qrySeqName = argv[i];
  i++;
  if (i >= argc) {
    PrintUsage();
    exit(1);
  }
  coordsFile = argv[i];
}


void PrintUsage() {
  std::cout << "winaln align two sequences and their surrounding regions. " << std::endl;
  std::cout << "usage: winaln [-m scorematfile] refseq qryseq coordsfile extendwidth " << std::endl;
}

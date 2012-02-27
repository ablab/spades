/***************************************************************************
 * Title:          MakeGrid.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <string>

// My tools for sequence alignment

#include "align/alignutils.h"
#include "Alignment.h"
#include "AlignmentPrinter.h"
#include "ShuffleAlign.h"

// Sequence modification utils. 
#include "SeqReader.h"
#include "DNASequence.h"
#include "SeqUtils.h"
#include "utils.h"

#include "InversionBins.h"

// Emboss alignment parsing
#include "emboss/EmbossAlign.h"
#include "emboss/EmbossAlignment.h"


#include "SimpleStats.h"

void ExtractSeqNames(std::vector<DNASequence*> &sequences,
		     StringVector &seqNames);

void FindBinToSeqMap(StringVector &binSeqNames,
		     StringVector &allSeqNames,
		     IntVector &binToAllMap,
		     IntVector &allToBinMap);

void PrintCharMat(StringVector &binSeqNames, StringVector &allSeqNames,
		  IntVector &binToAllMap, IntVector &allToBinMap,
		  FloatMatrix &forScores, FloatMatrix &revScores,
		  FloatMatrix &forIdentity, FloatMatrix &revIdentity,
		  FloatMatrix &meanRandScores,
		  double alignThreshold,
		  std::ofstream &charMatOut,
		  std::ofstream &alignScoreOut);


int main(int argc, char* argv[]) {
  std::string fastaBinFileName;
  std::string speciesFileName;
  std::string charsOutName;
  std::string scoreOutName;
  double alignK = 2;
  if (argc < 5) {
    std::cout << "usage: makegrid fastaBinFile species_list_file "
	      << " charsOutFile scoreOutFile [rand_align_factor] [-verbose vbout]" << std::endl;
    std::cout << "	fastabinfile - A list of sequences to find the relative order "<< std::endl;
    std::cout << "	specieslist  - A list of all species so that a m x n grid may be generated " 
	      << std::endl
	      << "                     for the scores" << std::endl;
    std::cout << "	rrand        - Align factor is the multiple above the mean " << std::endl;
    std::cout << " 	               align score that is considered to be significant " << std::endl;
    exit(0);
  }
  fastaBinFileName = argv[1];
  speciesFileName  = argv[2];
  charsOutName     = argv[3];
  scoreOutName     = argv[4];
  int argi = 5;
  if (argc > 5) {
    alignK = atof(argv[argi]);
  }
  std::string verboseFileName = "";
  std::ofstream vbOut;
  
  while (argi < argc) {
    if (strcmp(argv[argi], "-verbose") == 0) {
      verboseFileName = argv[++argi];
      openck(verboseFileName, vbOut, std::ios::out);
    }
    ++argi;
  }
  
  std::ifstream speciesFile;
  openck(speciesFileName, speciesFile);
  StringVector species;
  ReadSpeciesLine(speciesFile, species);
  std::ifstream fastaBinFile;
  openck(fastaBinFileName, fastaBinFile);

  std::ofstream charsOutFile, scoreOutFile;
  openck(charsOutName, charsOutFile, std::ios::out);
  openck(scoreOutName, scoreOutFile, std::ios::out);
  std::vector<DNASequence *> forwardSequences;
  std::vector<DNASequence *> reverseSequences;

  DNASequence *newSeq, *revSeq;
  newSeq = revSeq = NULL;
  SeqReader  seqReader(&fastaBinFile);
  while (seqReader.GetSeq(newSeq,SeqReader::noConvert)) {
    newSeq->HardMask();
    forwardSequences.push_back(newSeq);
    revSeq = new DNASequence;
    revSeq->SetAscii();
    revSeq->HardMask();
    MakeRC(*newSeq, *revSeq);
    reverseSequences.push_back(revSeq);
  }
  StringVector binSpecNames;
  IntVector binToAllMap, allToBinMap;
  ExtractSeqNames(forwardSequences, binSpecNames);
  FindBinToSeqMap(binSpecNames, species, binToAllMap, allToBinMap);

  char *mcsrc = getenv("MCSRC");
  std::string scoreMatName;
  std::string srcDir = mcsrc;
  scoreMatName = std::string(srcDir) + "/comparative/NCBIMAT.short";
  Score scoreMat(scoreMatName, 5, 2);
  

  ssize_t s1, s2;
  FloatMatrix forScores, revScores, forIdentity, revIdentity, meanRandScores;
  CreateMatrix(forScores, forwardSequences.size(), forwardSequences.size());
  CreateMatrix(revScores, forwardSequences.size(), forwardSequences.size());
  CreateMatrix(forIdentity, forwardSequences.size(), forwardSequences.size());
  CreateMatrix(revIdentity, forwardSequences.size(), forwardSequences.size());
  CreateMatrix(meanRandScores, forwardSequences.size(), forwardSequences.size());

  double pctIdentity;

  std::string alignString;
  for (s1 = 0; s1 < forwardSequences.size() - 1; s1++) {
    for (s2 = s1 + 1; s2 < forwardSequences.size(); s2++ ) {
      EmbossAlignment forwardAlignment, reverseAlignment;
      char* home;
      home = getenv("HOME");
     
      if (home == NULL ){ 
        std::cout << "ERROR: environment variable HOME is not defined " << std::endl;
	exit(1);
      }
      /*
	ssize_t *locations = NULL;
	locations = new ssize_t[forwardSequences[s1]->length];
	std::cout << "Affine rrunning local align" << std::endl;
	double myScore;
	myScore = AffineLocalAlign(*forwardSequences[s1], 
	*forwardSequences[s2],
				       scoreMat, locations); 
      delete []locations;
      std::cout << "done " << std::endl;
      */
      std::string homeStr = home;
      alignString = 
	homeStr + "/software/emboss/bin/matcher -datafile=" + 
	homeStr + 
	"/projects/mcsrc/comparative/NCBIMAT  -aformat3=srspair -gapopen=5 -gapextend=2 ";
      EmbossAlign(alignString, 
		  *forwardSequences[s1], *forwardSequences[s2],
		  forwardAlignment);
      
      ssize_t i;
      forIdentity[s1][s2] = forwardAlignment.similarity;
      forIdentity[s2][s1] = forwardAlignment.similarity;
      forScores[s1][s2] = forwardAlignment.alignScore;
      forScores[s2][s1] = forwardAlignment.alignScore;
      if (verboseFileName != "") {
	//vbOut
	std::cout << "forward s1: " << binSpecNames[s1] << " s2: " << binSpecNames[s2]
	      << " sim : " << forwardAlignment.similarity 
	      << " ident: " << forwardAlignment.pctIdentity 
	      << " score: " << forwardAlignment.alignScore 
	      << " length: " << forwardAlignment.length << std::endl;
      }

      alignString = 
	homeStr + "/software/emboss/bin/matcher -datafile=" + 
	homeStr +
	"/projects/mcsrc/comparative/NCBIMAT  -aformat3=srspair -gapopen=5 -gapextend=2 "; 
      EmbossAlign(alignString, 
		  *forwardSequences[s1], *reverseSequences[s2],
		  reverseAlignment);
     
      std::vector<double> alignScores;
      ssize_t nShuffles = 5;
      ShuffleAlign(*forwardSequences[s1], *reverseSequences[s2],
		   forwardSequences[s1]->length/2, nShuffles, alignScores);
      revIdentity[s1][s2] = reverseAlignment.similarity;
      revIdentity[s2][s1] = reverseAlignment.similarity;
      revScores[s1][s2] = reverseAlignment.alignScore;                                            
      revScores[s2][s1] = reverseAlignment.alignScore;					   
      double meanRandAlign, varRandAlign = 0;
      GetMeanVar(alignScores, meanRandAlign, varRandAlign);
      meanRandScores[s1][s2] = meanRandAlign;
      meanRandScores[s2][s1] = meanRandAlign;
      if (verboseFileName != "") {
	//	vbOut
	std::cout << "reverse s1: " << binSpecNames[s1] << " s2: " << binSpecNames[s2] 
	      << " sim : " << reverseAlignment.similarity 
	      << " ident: " << reverseAlignment.pctIdentity 
	      << " score: " << reverseAlignment.alignScore
	      << " length: " << reverseAlignment.length 
	      << " random: " << meanRandAlign << std::endl;
	vbOut << std::endl;
      }

      //      locations = NULL;
    }
  }
  // Build a distribution on the scores.
  
  // Complie a list of the top scores
  FloatVector maxScores;
  ssize_t i, j;
  double maxScore = -1;
  FloatVector maxRowScores;
  for (i = 0; i < forScores.size(); i++) {
    maxRowScores.push_back(-1);
    for ( j = 0; j < forScores[i].size(); j++ ) {
      if (i != j) {
	if (forScores[i][j] > revScores[i][j])
	  maxScores.push_back(forScores[i][j]);
	else
	  maxScores.push_back(revScores[i][j]);
	if (forScores[i][j] > maxScore) 
	  maxScore = forScores[i][j];
	if (revScores[i][j] > maxScore) 
	  maxScore = revScores[i][j];
	if (forScores[i][j] > revScores[i][j])
	  if ( maxRowScores[i] < forScores[i][j])
	    maxRowScores[i] = forScores[i][j];
	  else if (maxRowScores[i] < revScores[i][j])
	    maxRowScores[i] = revScores[i][j];
      }
    }
  }

  double meanMaxScore, varMaxScore, stdevMaxScore;
  double meanMaxRowScore, varMaxRowScores, stddevMaxRowScores;
  GetMeanVar(maxScores, meanMaxScore, varMaxScore);
  GetMeanVar(maxRowScores, meanMaxRowScore, varMaxRowScores);
  stdevMaxScore = sqrt(varMaxScore);
  stdevMaxScore = sqrt(varMaxScore);
  double alignThreshold;
  PrintCharMat(binSpecNames, species,
	       binToAllMap, allToBinMap,
	       forScores, revScores,
	       forIdentity, revIdentity,
	       meanRandScores,
	       alignK,// arbitrary minimum score.  This should be changed.
	       charsOutFile,  scoreOutFile);
}

void ExtractSeqNames(std::vector<DNASequence*> &sequences,
		     StringVector &seqNames) {

  ssize_t i;
  ssize_t dot;
  std::string seqName;
  for (i = 0; i < sequences.size(); i++ ) {
    dot = sequences[i]->namestr.find(".", 0);
    if (dot != sequences[i]->namestr.npos) {
      seqName = sequences[i]->namestr.substr(0, dot);
    }
    else {
      seqName = sequences[i]->namestr;
    }
    seqNames.push_back(seqName);
  }
}

void FindBinToSeqMap(StringVector &binSeqNames,
		     StringVector &allSeqNames,
		     IntVector &binToAllMap,
		     IntVector &allToBinMap) {
  ssize_t b, s;
  allToBinMap.resize(allSeqNames.size());
  for (s = 0; s < allSeqNames.size(); s++ ) {
    allToBinMap[s] = -1;
  }

  for (b = 0; b < binSeqNames.size(); b++ ) {
    for (s = 0; s < allSeqNames.size(); s++ ) {
      if (allSeqNames[s] == binSeqNames[b]) {
	binToAllMap.push_back(s);
	allToBinMap[s] = b;
	break;
      }
    }
  }
}

void PrintCharMat(StringVector &binSeqNames, StringVector &allSeqNames,
		  IntVector &binToAllMap, IntVector &allToBinMap,
		  FloatMatrix &forScores, FloatMatrix &revScores,
		  FloatMatrix &forIdenties, FloatMatrix &revIdenties,
		  FloatMatrix &meanRandScores,
		  double alignThreshold,
		  std::ofstream &charMatOut,
		  std::ofstream &alignScoreOut) {
  // Print the character matrix for each species as the reference 
  // species.  
  ssize_t binIndex, seqIndex;
  double forScore, revScore;
  double forIdentity, revIdentity;
  ssize_t seqBinIndex;
  ssize_t binSeqIndex;
  ssize_t binIndex2, seqIndex2;

  // Print the m x n character matrix, where m is the number of species
  // where the inversion was found and n is the total number of species.
  for (seqIndex = 0; seqIndex < allSeqNames.size(); seqIndex++ ) {
    binIndex = allToBinMap[seqIndex];
    if (binIndex == -1)
      continue;
    charMatOut << binSeqNames[binIndex] << "\t";
    for (seqIndex2 = 0; seqIndex2 < allSeqNames.size(); seqIndex2++ ) {
      binIndex2 = allToBinMap[seqIndex2];
      if (binIndex2 == binIndex) {
	charMatOut << " 1";
      }
      else if (binIndex2 == -1) {
	charMatOut << " 2";
      }
      else {
	forScore = forScores[binIndex][binIndex2];
	revScore = revScores[binIndex][binIndex2];
	// case 1. Both alignments are high-scoring, this may happen if there are 
	// multiple rearrangements
	if (forScore > alignThreshold*meanRandScores[binIndex][binIndex2] and
	    revScore > alignThreshold*meanRandScores[binIndex][binIndex2] ) {
	  charMatOut << " 2";
	}
	else if (forScore > alignThreshold*meanRandScores[binIndex][binIndex2]) {
	  charMatOut << " 1";
	}
	else if (revScore > alignThreshold*meanRandScores[binIndex][binIndex2]) {
	  charMatOut << " 0";
	}
	else {
	  charMatOut << " 2";
	}
      }
    }
    charMatOut << std::endl;
  }
  // print the m x m matrix of scores. We don't need to worry about printing
  // the alignment scores of the sequences not found for the inversion.
  for (seqIndex = 0; seqIndex < allSeqNames.size(); seqIndex++ ) {
    binIndex = allToBinMap[seqIndex];
    if (binIndex == -1) 
      continue;
    for (seqIndex2 = 0; seqIndex2 < allSeqNames.size(); seqIndex2++ ) {
      binIndex2 = allToBinMap[seqIndex2];
      if (binIndex2 == -1) continue;
      if (binIndex == binIndex2) {
	alignScoreOut << " 1";
      }
      else {
	forScore = forScores[binIndex][binIndex2];
	revScore = revScores[binIndex][binIndex2];
	forIdentity = forIdenties[binIndex][binIndex2];
	revIdentity = revIdenties[binIndex][binIndex2];
	// case 1. Both alignments are high-scoring, this may happen if there are 
	// multiple rearrangements
	/*
	  if (forScore > alignThreshold*meanRandScores[binIndex][binIndex2] and
	    revScore > alignThreshold*meanRandScores[binIndex][binIndex2] ) {
	  alignScoreOut << " " << -99;
	}
	// Case 2.  Forward score is higher.
	else if (forScore > alignThreshold*meanRandScores[binIndex][binIndex2]) {
	  alignScoreOut << " " << forScore;
	}
	// Case 3.  Reverse score is higher.
	else if (revScore > alignThreshold*meanRandScores[binIndex][binIndex2] ) {
	  alignScoreOut << " " << -revScore;
	}
		else {
	  // Both scores are too low to call.
	  alignScoreOut << " " << -99;
	  }*/
	if (forScore > revScore) {
	  alignScoreOut << " " << forScore;
	}
	// Case 3.  Reverse score is higher.
	else {
	  alignScoreOut << " " << -revScore;
	}
      }
    }
    alignScoreOut << std::endl;
  }
}
		  


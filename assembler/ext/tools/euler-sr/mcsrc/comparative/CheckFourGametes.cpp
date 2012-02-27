/***************************************************************************
 * Title:          CheckFourGametes.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <unistd.h>

#include <set>

#include "utils.h"
#include "InversionChars.h"
#include "InversionBins.h"
#include "FourGameteUtils.h"

ssize_t verbose;

ssize_t FindFourGameteViolatingSpecies(IntMatrix & charMat, std::vector<std::string> &species,  PosVect &inv1, PosVect &inv2, 
				   PosVect &s1, PosVect &s2, PosVect &s3, PosVect &s4);


ssize_t NextViolation(std::vector<ssize_t> &vind, ssize_t &curInd);

void InitEnv(int argc, char* argv[], 
	     std::string &binFileName,
	     ssize_t &findViolatingSpecies,
	     ssize_t &verbose,
	     ssize_t &checkConsensus,
	     std::string &fixedBinFile);

void FixBins(BinMap &binMap, 
	     std::vector<ssize_t> &startPos, std::vector<ssize_t> &endPos, 
	     std::vector<ssize_t> &vind1, std::vector<ssize_t> &vind2);


void PrintUsage();

int main(int argc, char* argv[]) {
  std::string binFileName, fixedBinFile;
  ssize_t findViolatingSpecies = 0;
  ssize_t checkConsensus = 1;
  ssize_t verbose = 0;

  binFileName = "";
  fixedBinFile = "";
  InitEnv(argc, argv, binFileName, findViolatingSpecies, verbose, checkConsensus, fixedBinFile);

  // Read in the table of species.
  StringVector species;
  std::vector<ssize_t> startPos, endPos;
  BinMap binnedInversions;
  ReadBinFile(binFileName, species, binnedInversions);

  std::vector<ssize_t> inv1, inv2, spec1, spec2, spec3, spec4;
  //  FindFourGameteViolatingCharacters(charMat, species, inv1, inv2);
  BinMap consensusMap;

  if (checkConsensus) {
    InversionList consensus;
    std::vector<ssize_t> mult;
    
    GetBinConsensus(binnedInversions, species.size(), startPos, endPos, consensus, mult);
    
    ssize_t invIdx;
    for (invIdx = 0; invIdx < consensus.size(); invIdx++)
      consensusMap[0].push_back(consensus[invIdx]);
    
    IntMatrix conMat;
    BuildCharMat(consensusMap, conMat);
    std::string consensusFileName = "consensus.txt";
    std::ofstream charMatFile;
    openck(consensusFileName, charMatFile, std::ios::out);
    PrintBins(species, consensusMap, charMatFile);
    charMatFile.close();

    std::cout << "checking four gamete with consensus " << std::endl;
    if (findViolatingSpecies)
      FindFourGameteViolatingSpecies(conMat, species, inv1, inv2, spec1, spec2, spec3, spec4);
    /*
      else
      FindFourGameteViolatingCharacters(conMat, inv1, inv2);
    */
    
    if (fixedBinFile != "") {
      startPos.resize(conMat.size());
      endPos.resize(conMat.size());
    }
  }
  else {
    IntMatrix charMat;
    BuildCharMat(binnedInversions, charMat);
    if (findViolatingSpecies)
      FindFourGameteViolatingSpecies(charMat, species, inv1, inv2, spec1, spec2, spec3, spec4);
    /*    else
	  FindFourGameteViolatingCharacters(charMat, inv1, inv2);
    */
  }
  /*
    ssize_t curInd = 0;
    while (curInd < inv1.size() ) {
    NextViolation(inv1, curInd);
    std::cout << inv1[curInd] << " " << inv2[curInd] << std::endl;
    ++curInd;
    }
    

    ssize_t invIdx;
    for (invIdx = 0; invIdx < inv1.size(); invIdx++) {
    std::cout << inv1[invIdx] << " " << inv2[invIdx] << std::endl;
    }
  */


  if (fixedBinFile != "" ) {
    std::ofstream binout;
    std::cout << "writing binned inversions " << fixedBinFile << std::endl;
    openck(fixedBinFile, binout, std::ios::out);

    if (checkConsensus) {
      FixBins(consensusMap, startPos, endPos, inv1, inv2);
      PrintBins(species, consensusMap, binout, 0);
    }
    else {
      FixBins(binnedInversions, startPos, endPos, inv1, inv2);
      PrintBins(species, binnedInversions, binout, 0);
    }
    binout.close();
  }
  return 0;
}



ssize_t FindFourGameteViolatingSpecies(IntMatrix & charMat, std::vector<std::string> &species, PosVect &inv1, PosVect &inv2, 
				   PosVect &spec1, PosVect &spec2, PosVect &spec3, PosVect &spec4) {

  ssize_t numInversions = charMat.size();
  if (numInversions <= 0) 
    return 0;

  ssize_t numSpecies = charMat[0].size();
  
  // Find species/traits that violate the 4-gamete test.

  ssize_t i1, i2;  // Inversion features

  ssize_t s1, s2, s3, s4;
  
  ssize_t aa, ab, ba, bb;
  // Look through every pair of inversions.
  for (i1 = 0; i1 < numInversions-1; i1++) {
    for (i2 = i1 + 1; i2 < numInversions; i2++) {
      // For every pair of inversions, try to find four species that violate the 4-gamete test.

      for (s1 = 0; s1 < numSpecies - 3; s1++ ) {
	for (s2 = s1+1; s2 < numSpecies-2; s2++) {
	  for (s3 = s2 +1; s3 < numSpecies-1; s3++) {
	    for (s4 = s3 + 1; s4 < numSpecies; s4++) {
	      aa = 0; ab = 0; ba = 0; bb = 0;
	      if ((charMat[i1][s1] == 0 and charMat[i2][s1] == 0) or 
		  (charMat[i1][s2] == 0 and charMat[i2][s2] == 0) or 
		  (charMat[i1][s3] == 0 and charMat[i2][s3] == 0) or 
		  (charMat[i1][s4] == 0 and charMat[i2][s4] == 0))
		aa = 1;	          		           
	      if ((charMat[i1][s1] == 0 and charMat[i2][s1] == 1) or 
		  (charMat[i1][s2] == 0 and charMat[i2][s2] == 1) or 
		  (charMat[i1][s3] == 0 and charMat[i2][s3] == 1) or 
		  (charMat[i1][s4] == 0 and charMat[i2][s4] == 1))
		ab = 1;	          		           
	      if ((charMat[i1][s1] == 1 and charMat[i2][s1] == 0) or 
		  (charMat[i1][s2] == 1 and charMat[i2][s2] == 0) or 
		  (charMat[i1][s3] == 1 and charMat[i2][s3] == 0) or 
		  (charMat[i1][s4] == 1 and charMat[i2][s4] == 0))
		ba = 1;	          		           
	      if ((charMat[i1][s1] == 1 and charMat[i2][s1] == 1) or 
		  (charMat[i1][s2] == 1 and charMat[i2][s2] == 1) or 
		  (charMat[i1][s3] == 1 and charMat[i2][s3] == 1) or 
		  (charMat[i1][s4] == 1 and charMat[i2][s4] == 1))
		bb = 1;
	      //	      std::cout << aa << " " << ab << " " << ba << " " << bb << std::endl;
	      if (aa and ab and ba and bb ) {
		// Found a violation of 4-gamete condition. Store it.
		inv1.push_back(i1); inv2.push_back(i2);
		spec1.push_back(s1); spec2.push_back(s2); spec3.push_back(s3); spec4.push_back(s4);
	      }
	    }
	  }
	}
      }
    }
  }
  return 0;
}

void InitEnv(int argc, char* argv[], 
	     std::string &binFileName,
	     ssize_t &findViolatingSpecies,
	     ssize_t &verbose,
	     ssize_t &checkConsensus,
	     std::string &fixedBinFile) {
  ssize_t copt;
  verbose = 0;
  while ( (copt=getopt(argc, argv, "c:C:vf:")) != EOF){
    switch(copt) {
    case 's':
      if (strcmp(optarg, "s") == 0)
	findViolatingSpecies = 1;
      else if (strcmp(optarg, "c") == 0)
	findViolatingSpecies = 0;
      continue;
    case 'v':
      ++verbose;
      continue;
    case 'C':
      checkConsensus = atoi(optarg);
      continue;
    case 'f':
      fixedBinFile = optarg;
      std::cout << "fxed bin file: " << fixedBinFile << std::endl;
      continue;
    default:
      PrintUsage();
      exit(1);
    }
  }
  ssize_t i = optind;
  if (i >= argc) {
    std::cout <<" not enough arguments " << optind << std::endl;
    PrintUsage();
    exit(0);
  }
  binFileName = argv[i];
}


void PrintUsage() {
  std::cout << "usage: fgc binFileName [-v] [-C {0|1}] [-c {s|c}] [-f fixedBinFile]" << std::endl;
  std::cout << "-v adds verbossity per -v " << std::endl
	    << "-C {0,1} (1) find consensus of bins before checking validity.  " << std::endl
	    << "-c {s,c} (c) find the individual species (s) or characters (c) responsible for violation. " << std::endl;
  std::cout << "-f fixedBinFile.  Once violating chars are found, remove them, and output the bins" << std::endl
	    << "   to fixedBinFile " << std::endl;
}

ssize_t NextViolation(std::vector<ssize_t> &vind, ssize_t &curInd) {
  while (curInd < vind.size()-1 and 
	 vind[curInd] == vind[curInd+1])
    curInd++;
  return curInd;
}

void FixBins(BinMap &binMap, 
	     std::vector<ssize_t> &startPos, std::vector<ssize_t> &endPos, 
	     std::vector<ssize_t> &vind1, std::vector<ssize_t> &vind2) {
  ssize_t binIndex, vIndex;
  
  binIndex = 0; vIndex = 0;
  // vind1 and vind2 contain the indices of the violating vertices.  
  // don't bother looking through vind2 for now, since it is totally 
  // contained in vind1.

  // Look through all of the feature rows in binmap, they are enumerated 1 .. max bin int.
  // vind1 indexes.

  BinMap::iterator binIt;
  InversionList::iterator invIt;
  ssize_t totalErased = 0;
  // Look through the features for each bin.
  for (binIt = binMap.begin(); binIt != binMap.end(); ) {
    for (invIt = (*binIt).second.begin(); 
	 invIt != (*binIt).second.end();) {
      if (binIndex == vind1[vIndex]) {
	// This bin is listed in the violating indices list
	// Get rid of it.
	invIt = (*binIt).second.erase(invIt);
	vIndex = NextViolation(vind1, vIndex) + 1;
	//	startPos.erase(startPos.begin() + binIndex);
	//	endPos.erase(endPos.begin() + binIndex);
	++totalErased;
      }
      else {
	++invIt;
      }
      ++binIndex;
    }
    if ((*binIt).second.size() == 0) {
      BinMap::iterator erasedBin;
      erasedBin = binIt;
      ++binIt;
      binMap.erase(erasedBin);
      ssize_t erasedIndex = (*binIt).first;
      //      startPos.erase(startPos.begin() + erasedIndex);
      //		endPos.erase(endPos.begin() + eras);
    }
    else
      ++binIt;
  }
}

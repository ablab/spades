/***************************************************************************
 * Title:          mergeblocks.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/28/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <stack>
#include <list>
#include <string>
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>

#include "StripGen.h"

void printusage();
void initenv(int argc, char *argv[], ssize_t &startDist, ssize_t &endDist, ssize_t &stepSize, 
	     std::string &inpfile, std::string &outfile, 
	     ssize_t &inType, ssize_t &withSeq, ssize_t &withDir);

int main(int argc, char* argv[]) {

  std::string inFileName, outFileName;
  ssize_t startDist, endDist;
  ssize_t distThresh;
  ssize_t stepSize = 1;

  std::ifstream inFile;
  std::ofstream outFile;

  startDist = endDist = 1;
  char seq[100]; // sequence.
  T_Strips strips;
  ssize_t inType, withSeq, withDir;
  initenv(argc, argv, startDist, endDist, stepSize,
	  inFileName, outFileName, inType, withSeq, withDir);

  inFile.open(inFileName.c_str());
  if (!inFile.good()) {
    std::cout << "Could not open " << inFileName << std::endl;
    exit(1);
  }
  ssize_t numStrips, numDiscardedStrips;
  numStrips = numDiscardedStrips = 0;
  ssize_t lineNumber = 0;
  ssize_t stripSize = 0;
  ssize_t numSkipped = 0;
  ssize_t numRead = 0;
  numStrips = GetStrips(inFile, strips, numRead);

  inFile.close();
  T_Strips::iterator cur, prev, end;
  std::cout << "read: " << numRead << " counted : " << numStrips << " strips and discarded "
	    << numDiscardedStrips << std::endl;
  for (distThresh = startDist; distThresh <= endDist; distThresh += stepSize) {
    numDiscardedStrips = DiscardSmallStrips(strips, distThresh);
    T_StripQQry stripsSortedByQry;
    BuildStripPqueue(strips, stripsSortedByQry);
    AssignEnumerations(stripsSortedByQry);
    std::cout << "iteration " <<  distThresh - startDist + 1 << " size: " << strips.size() << " discarded: " << numDiscardedStrips << std::endl;
  }
  outFile.open(outFileName.c_str());
  if (!outFile.good() ){
    std::cout << "could not open " << outFileName << std::endl;
    exit(1);
  }

  PrintEnumerations(strips, outFile);
  outFile.close();

  return 0;
}

void initenv(int argc, char *argv[], 
	     ssize_t &startDist, ssize_t &endDist, ssize_t &stepSize,
	     std::string &inpfile, std::string &outfile, ssize_t &intype, ssize_t &withSeq, ssize_t &withDir){
  ssize_t copt;
  intype = 0;
  withSeq = 1;
  withDir = 1;
  while ( (copt=getopt(argc, argv, "i:o:d:D:tsuS:")) != EOF) {
    switch(copt) {
      case 'i':
	inpfile = optarg;
	continue;
      case 'o':
	outfile = optarg;
	continue;
    case 't':
      intype = 1;
      continue;
    case 'S':
      sscanf(optarg, "%d", &stepSize);
      continue;
    case 's':
      withSeq = 0;
      continue;
    case 'u':
      withDir = 0;
      continue;
    case 'd':
	sscanf(optarg,"%d", &startDist);
	continue;
    case 'D':
	sscanf(optarg,"%d", &endDist);
	continue;
    default:
      printusage();
      exit(1);
    }
  }
  if (inpfile == "") {
    printusage();
    exit(1);
  }
  /*  if (startDist == -1 ) {
    std::cout << " Start distance not specified " << std::endl;
    printusage();
    exit(0);
    }*/
  if (endDist == 1)
    endDist = startDist;
}
void printusage() {
  std::cout << "usage:  mergeblocks  -i seq1  [-d distThresh] -o output [-t] " << std::endl;
  std::cout << "   -i input sequence1 ... : sequences to find unique probes in." << std::endl;
  std::cout << "   -d startDistSize       : smallest size strip to delete (default 1). " << std::endl;
  std::cout << "   -D endDist             : largest distance threshold to check for " << std::endl;
  std::cout << "   -S stepSize            : stepSize to use between distances (default 1) " << std::endl;
  std::cout << "   -o output              : output file " << std::endl;
  std::cout << "   -t                     : iTerate, use the output of a previous run as input. " << std::endl;
  std::cout << "   -s                     : do not read a sequence at the end of each line." << std::endl;
  std::cout << "   -u                     : unsigned strips (default signed)." << std::endl;
}

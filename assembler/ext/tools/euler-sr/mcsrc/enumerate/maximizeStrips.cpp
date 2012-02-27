/***************************************************************************
 * Title:          maximizeStrips.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "MaximizeStrip.h"
#include "StripGen.h"
#include "TupleLib.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <queue>
#include <stack>
#include <list>
#include <string>
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>


void printusage();
void initenv(int argc, char *argv[], 
	     std::string &inpfile, 
	     std::string &outfile,
	     ssize_t &window, ssize_t &threshold, ssize_t &merge);


int main(int argc, char* argv[]) {

  
  std::string inputFileName, outputFileName;
  ssize_t window, threshold, mergeStrips;
  mergeStrips = 0;
  window = -1;
  initenv(argc, argv, inputFileName, outputFileName, window, threshold, mergeStrips);
  T_Strips *strips;
  
  strips = new T_Strips;
  std::ifstream inFile;
  
  inFile.open(inputFileName.c_str());
  if ( ! inFile.good() ) {
    std::cout << "Could not open infile " << std::endl;
    exit(0);
  }
  ssize_t numRead;
  // Read the strips in from the file.  Condense adjacent strips.
  GetStrips(inFile, *strips, numRead, mergeStrips);
  if (window == -1)
    window = strips->size();

  T_StripQQry stripsSortedByQry;
  // Re-enumerate strips
  
  PrintEnumerations(*strips, std::cout);
  MaximizeStrips(strips, window, threshold);
  BuildStripPqueue(*strips, stripsSortedByQry);
  AssignEnumerations(stripsSortedByQry);
  if (outputFileName != "") {
    std::ofstream out;
    out.open(outputFileName.c_str());
    if (!out.good()) {
      std::cout << "Could not open " << outputFileName << std::endl;
      exit(0);
    }
    PrintEnumerations(*strips, out);
    out.close();
  }
  else 
    PrintEnumerations(*strips, std::cout);
  return 0;
}


void initenv(int argc, char *argv[], std::string &inpfile, std::string &outfile, ssize_t &window, ssize_t &threshold, ssize_t &mergeStrips) {
  ssize_t copt;
  while ( (copt=getopt(argc, argv, "i:o:w:t:m")) != EOF) {
    switch(copt) {
    case 'i':
      inpfile = optarg;
      continue;
    case 'o':
      outfile = optarg;
      continue;
    case 'w':
      sscanf(optarg, "%d", &window);
      continue;
    case 't':
      sscanf(optarg, "%d", &threshold);
      continue;
    case 'm':
      mergeStrips = 1;
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
}
void printusage() {
  std::cout << "usage:  mergeblocks  -i seq1  [-d distThresh] -o output [-t] " << std::endl;
  std::cout << "   -i input sequence1 ... : sequences to find unique probes in." << std::endl;
  std::cout << "   -o output              : output file " << std::endl;
  std::cout << "   -w windowsize " << std::endl;
  std::cout << "   -t distthreshold " << std::endl;
  std::cout << "   -m    (merge adjacent strips in input file) " << std::endl;
}




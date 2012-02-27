/***************************************************************************
 * Title:          stripsToGlue.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
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
#include <cmath>
#include <string>
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <math.h>

#include "StripGen.h"

void printusage();
void initenv(int argc, char *argv[], 
	     std::string &inpfile, 
	     std::string &outfile,
	     std::string &indexA, 
	     std::string &indexB);

int main(int argc, char* argv[]) {
  std::string inputFileName, outputFileName;
  std::string indexA("@");
  std::string indexB("!");
  initenv(argc, argv, inputFileName, outputFileName, indexA, indexB);

  T_Strips strips;
  
  std::ifstream inFile;
  
  inFile.open(inputFileName.c_str());
  if ( ! inFile.good() ) {
    std::cout << "Could not open infile " << std::endl;
    exit(0);
  }
  ssize_t numRead;
  GetStrips(inFile, strips, numRead, 0);

  std::ofstream outFile;
  outFile.open(outputFileName.c_str());
  if (!outFile.good()) {
    std::cout << "could not open " << outputFileName << std::endl;
    exit(0);
  }
  T_Strips::iterator it,end;
  end = strips.end();
  ssize_t step, qryPos, refPos;

  double stripRatio;
  double bothRatio, refRatio, qryRatio;
  double bothDiff, refDiff, qryDiff;
  std::string sign;
  for (it = strips.begin(); it != end; ++it) {
    if ((*it).endQryPos < (*it).startQryPos)
      step = -1;
    else
      step = 1;

    stripRatio = double((*it).endRefPos - (*it).startRefPos) /
             (abs((*it).endQryPos - (*it).startQryPos));
		// TODO: if we expand size, change abs to szabs

    // Work out for positive step direction
    refPos = (*it).startRefPos;
    qryPos = (*it).startQryPos;
    outFile << " " << refPos << "\t" << indexA << ";\t" << qryPos << "\t" << indexB << ";" << std::endl;
    refPos++;
    qryPos += step;
    if (Sign((*it).startQryEnum) == -1) 
      sign = "-1";
    else
      sign = "";
    while (refPos < (*it).endRefPos && 
	   ((step == 1 && qryPos < (*it).endQryPos) || (step == -1 && qryPos > (*it).endQryPos))) {
      // Decide what steps to make
			// TODO: casts & arithmetic could be simplfied
      bothRatio = (double(refPos+1) - (*it).startRefPos) / std::fabs((double(qryPos + step) - (*it).startQryPos));
      refRatio  = (double(refPos+1) - (*it).startRefPos) / std::fabs(double(qryPos - (*it).startQryPos));
      qryRatio  = (double(refPos) - (*it).startRefPos) / std::fabs((double(qryPos + step) - (*it).startQryPos));
      bothDiff = std::fabs(stripRatio - bothRatio);
      refDiff  = std::fabs(stripRatio - refRatio);
      qryDiff  = std::fabs(stripRatio - qryRatio);
      
      // Move to the next step
      if (bothDiff <= refDiff && bothDiff <= qryDiff) {
	refPos++;
	qryPos += step;
	outFile << " " << refPos << " " << indexA << "; " << qryPos <<  " " << sign << indexB << ";" << std::endl;
      }
      else if (refDiff <= bothDiff && refDiff <= qryDiff) {
	refPos++;
      }
      else if (qryDiff <= refDiff && qryDiff <= bothDiff) {
	qryPos += step;
      }
    }
  }
  outFile.close();
  return 0;
}

void initenv(int argc, char *argv[], std::string &inpfile, std::string &outfile, std::string &indexA, std::string &indexB) {
  ssize_t copt;
  while ( (copt=getopt(argc, argv, "i:o:a:b:")) != EOF) {
    switch(copt) {
      case 'i':
	inpfile = optarg;
	continue;
      case 'o':
	outfile = optarg;
	continue;
    case 'a':
      indexA = optarg;
      continue;
    case 'b':
      indexB = optarg;
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
  std::cout << "usage:  mergeblocks  -i seq1 -o output [-a index] [-b index]" << std::endl;
  std::cout << "   -i strips file." << std::endl;
  std::cout << "   -o glue file " << std::endl;
  std::cout << "   -a index, index of first column " << std::endl;
  std::cout << "   -b index, index of second column " << std::endl;
}








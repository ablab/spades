/***************************************************************************
 * Title:          AlignEnds.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "emboss/EmbossAlign.h" 
#include "DNASequence.h"
#include "SeqReader.h"
#include <iostream>
#include <string>
#include <stdlib.h>
#include <limits.h>

int main(int argc, char* argv[]) {

  std::string file1, file2;

  ssize_t step, intv, il1, il2;  // interval lengths
  ssize_t nits;

  if (argc <= 5) {
    std::cout << "usage: alignends sequence " << std::endl;
    std::cout << "  look for self-similarity on the ends of a sequence " 
	      << std::endl; 
    exit(1);
  }

  file1 = argv[1];
   
  ssize_t verbose = 0;
  if (argc == 7) {
    if (strcmp(argv[6],"-v") == 0) {
      verbose = 1;
    }
  }
      
  DNASequence seq;
  SeqReader::GetSeq(file1, seq, SeqReader::noConvert);

  
  ScoreBandedAffineAlign(refSeq, qrySeq, // sequences to compare 
			    double initialScore,  // score of seed
			    double *scoreMat[5], double gapOpen, double gapExtend, // scoring parameters
			    double maxScore,  // negative score is good.
			    ssize_t *&locations, ssize_t &length, 
			    ssize_t &refAlignStart, ssize_t &qryAlignStart, // used for storing the resulting map
			    ssize_t refStartPos=0, ssize_t qryStartPos=0, 
			    ssize_t extDir=1
			    );

  

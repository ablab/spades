/***************************************************************************
 * Title:          LAVFile.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  10/03/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "LAVFile.h"
#include <sstream>
#include "align/alignutils.h"

const char* LAVFile::standardBlastzOpts = "blastz.v7 mouse.fasta hedgehog.fasta M=3\n     A    C    G    T\n    91 -114  -31 -123\n  -114  100 -125  -31\n   -31 -125  100 -114\n  -123  -31 -114   91\n  O = 400, E = 30, K = 3000, L = 3000, M = 3";


void LAVFile::ParseBlastzOpts(std::string options,
			      IntMatrix &scoreMat,
			      std::map<char, ssize_t> &keywordOptions) {

  /* parse a string of the form
"blastz.v7 /Users/mchaisso/projects/inversions/ENCODE/ENm001/mouse.fasta /Users/mchaisso/projects/inversions/ENCODE/ENm001/human.fasta
     A    C    G    T
    91 -114  -31 -123
  -114  100 -125  -31
   -31 -125  100 -114
  -123  -31 -114   91
  O = 400, E = 30, K = 3000, L = 3000, M = 0"

  */
  std::istringstream optstream(options);

  // get rid of the first line
  while(optstream) {
    if (optstream.get() == '\n') break;
  }
  if (!optstream) {
    std::cout << "Failed to parse blastz options: " << std::endl;
    std::cout << options << std::endl;
    exit(1);
  }
  ssize_t i, j;
  // pitch the acgt 
  unsigned char nucs[4];
  for (i = 0; i < 4; i++) {
    optstream >> nucs[i];
  }
  ssize_t score;
  for (i = 0; i < 4; i++) {
    for (j = 0; j < 4; j++) {
      optstream >> score;
      scoreMat[nuc_index[nucs[i]]][nuc_index[nucs[j]]] = score;
    }
  }

   
  char optchar;
  ssize_t optvalue;
  std::string junk;
  while (optstream) {
    optstream >> optchar >> junk >> optvalue;
    keywordOptions[optchar] = optvalue;
  }
}

/***************************************************************************
 * Title:          AlignRand.cpp 
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
#include "compatibility.h"
#include <iostream>
#include <string>
#include <algorithm>
#include <stdlib.h>
#include <limits.h>

int main(int argc, char* argv[]) {

  std::string file1, file2;

  ssize_t step, intv, il1, il2;  // interval lengths
  ssize_t nits;

  if (argc <= 5) {
    std::cout << "usage: alignrand sequence  intval1 intval2 step nits " << std::endl;
    std::cout << "  sample nonoverlapping fragments of len intval1 from sequence, " 
	      << " intval2 from sequence, and align them. " << std::endl; 
    std::cout << "  do this nits times " << std::endl;
    exit(1);
  }

  file1 = argv[1];
  il1   = atoi(argv[2]);
  il2   = atoi(argv[3]);
  step  = atoi(argv[4]);
  nits  = atoi(argv[5]);
   
  ssize_t verbose = 0;
  if (argc == 7) {
    if (strcmp(argv[6],"-v") == 0) {
      verbose = 1;
    }
  }
      
  DNASequence seq;
  SeqReader::GetSeq(file1, seq, SeqReader::noConvert);

  ssize_t i;
  ssize_t pos1, pos2;
  DNASequence frag1, frag2;
  frag1._ascii = 1;
  frag2._ascii = 1;
  EmbossAlignment alignment;
  std::string alignCommand;
  ssize_t s;
  double totScore;
  for (intv = il1; intv < il2; intv+= step) {
    totScore = 0;
    for (i = 0; i < nits; i++) {
      // Create the two positions

      pos1 = (ssize_t) floor((double(random()) / RANDOM_MAX) * (seq.length - intv));
      pos2 = pos1;
      while ( (pos2 <= pos1 and (pos2 + il2) >= pos1) or 
	      (pos1 <= pos2 and (pos1 + intv) >= pos2) ) {
				pos2 = (ssize_t) floor((double(random()) / RANDOM_MAX) * (seq.length - il2));
      }
    
      frag1.seq = &seq.seq[pos1];
      frag1.length = intv;

      frag2.seq = &seq.seq[pos2];
      frag2.length = il2;

      alignCommand = "~/projects/sw/bin/water -aformat srspair -gapOpen 10.0 -gapExtend 0.5 ";
      if (!EmbossAlign(alignCommand, frag1, frag2, alignment) ){ 
	std::cout <<" error with align command: " << alignCommand << std::endl;
	exit(1);
      }
      totScore += alignment.alignScore;
      if (verbose == 1) {
	std::cout << il1 << "\t" << pos1 << "\t" << il2 << "\t" << pos2 << "\t" 
	<< alignment.alignScore << std::endl;
      }
    }
    if (verbose == 0) {
      std::cout << intv << "\t" << totScore / nits << std::endl;
    }
  }
}

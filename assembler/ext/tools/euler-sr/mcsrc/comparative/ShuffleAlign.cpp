/***************************************************************************
 * Title:          ShuffleAlign.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "ShuffleAlign.h"
#include "compatibility.h"

#include <stdlib.h>
#include <cmath>
#include "emboss/EmbossAlign.h"
#include "emboss/EmbossAlignment.h"
#include "SeqUtils.h"
void MakeShuffledOrder(ssize_t nblocks, std::vector<ssize_t> &order);

void ShuffleAlign(DNASequence &seqA, DNASequence &seqB,
		  ssize_t nblocks,
		  ssize_t niter,
		  std::vector<double> &scores) {
  scores.clear();
  scores.resize(niter);
  ssize_t i;
  DNASequence shuffled;
  DNASequence rc;
  shuffled = seqA;
  std::vector<ssize_t> order;
  // Prepare the order

  // Shuffle the order
  ssize_t sbegin, send, slen, rbegin, rend;
  ssize_t step = seqA.length / nblocks;
  ssize_t b, j;
  char *home = getenv("HOME");
  std::string homeStr = home;
  std::string scoreMatName;
  scoreMatName = std::string(home) + "/projects/mcsrc/align/data/scoremat.txt";
  shuffled.namestr = "shuffled";
  rc._ascii = 1;
  shuffled._ascii = 1;
  MakeRC(seqA, rc);
  ssize_t useReverse;
  for (i = 0; i < niter; i++ ) {
    // scramble seq A
    rbegin = 0;
    rend   = std::min(seqA.length, rbegin + step);
    order.clear();
    MakeShuffledOrder(nblocks, order);
    useReverse = random() & 1;
    /*
    std::cout << "order: ";
    for (b = 0; b < nblocks; b++) {
      std::cout << " " << order[b];
    }
    std::cout << std::endl;
    */
    for (b = 0; b < nblocks; b++) {
      sbegin = order[b] * step;
      send   = std::min(sbegin + step, seqA.length);
      slen   = send - sbegin;
      for (j = 0; j < slen; j++) {
	if (useReverse) {
	  shuffled.seq[sbegin+j] = rc.seq[rbegin+j];
	}
	else {
	  shuffled.seq[sbegin+j] = seqA.seq[rbegin+j];
	}
      }
      rbegin += step;
      rend   = std::min(seqA.length, rbegin + step);
    }
    std::string alignString;
    EmbossAlignment alignment;
    char *home = getenv("HOME");
    alignString = 
      std::string(home) + "/software/emboss/bin/matcher -datafile=" + 
      homeStr +
      "/projects/mcsrc/comparative/NCBIMAT  -aformat3=srspair -gapopen=5 -gapextend=2 ";
     //std::cout << "aligning with " << alignString << std::endl;
    EmbossAlign(alignString, shuffled, seqB, alignment);
    //    std::cout << "got alignscore: " << alignment.alignScore << std::endl;
    scores[i] = alignment.alignScore;
  }
}

void MakeShuffledOrder(ssize_t nblocks, std::vector<ssize_t> &order) {
  
  ssize_t i,j;
  ssize_t found;
  std::vector<ssize_t> identity;
  identity.resize(nblocks);
  for (i = 0; i < nblocks; i++) identity[i] = i;
  ssize_t index;
  for (i = 0; i < nblocks; i++) {
    index = (ssize_t) std::floor((double(random()) / RANDOM_MAX) * identity.size()); 
    order.push_back(identity[index]);
    identity.erase(identity.begin() + index);
  }
}
  

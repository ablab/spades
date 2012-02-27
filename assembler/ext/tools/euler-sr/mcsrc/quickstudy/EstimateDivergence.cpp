/***************************************************************************
 * Title:          EstimateDivergence.cpp 
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
#include "mctypes.h"
#include "alignutils.h"
#include "SeqUtils.h"
#include "compatibility.h"
#include <iostream>
#include <string>
#include <set>
#include <stdlib.h>
#include <limits.h>
const char other_nucs[12] = {'a', 'c', 't', 
			     'g', 'c', 't', 
			     'g', 'a', 't', 
			     'g', 'a', 'c'};

char OtherNuc(char old);
ssize_t Random(ssize_t max);
void CopySeq(DNASequence &src, DNASequence &dest, ssize_t pos, ssize_t length);
void MutateSequence(DNASequence &original,
		    DNASequence &mutated,
		    double pctDivergence);

int main(int argc, char* argv[]) {

  std::string file1, file2;

  if (argc <= 3) {
    std::cout << "usage: alignrand sequence length pct_divergence [nits] " << std::endl;
    std::cout << "  Try to mutate a sequence with pct_divergence, and " 
	      << " and compare the expected and new alignment scores.";
    std::cout << "  do this nits times " << std::endl;
    std::cout << "outputs one line per sample with: " << std::endl;
    std::cout << "conserved mutated expected" << std::endl;
    exit(1);
  }
  ssize_t length;
  double pctDivergence;
  ssize_t nIts;
  file1 = argv[1];
  length= atoi(argv[2]);
  pctDivergence  = atof(argv[3]);
  nIts = 1000;
  if (argc > 4) {
    nIts = atoi(argv[4]);
  }
  DNASequence seq;
  SeqReader::GetSeq(file1, seq, SeqReader::noConvert);

  DNASequence orig, mutated;

  ssize_t it, p;
  ssize_t sourcePos;
  orig.Reset(length);
  mutated.Reset(length);
  ssize_t conservedScore, mutatedScore, expectedScore;
  // Initialize a default score matrix
  char* home;
  home = getenv("HOME");
  std::string scorematName = 
    std::string(home) + "/projects/mcsrc/align/data/blastzscoremat.txt"; 
  Score score(scorematName, 400,50);
  // calculate the total penalty
  ssize_t i, j;
  double totalPenalty = 0.0;
  for (i = 0; i < 4; i++) 
    for (j = 0; j < 4; j++) 
      if (i != j) totalPenalty+= score.scoreMat[i][j] - score.scoreMat[i][i];
  
  double averagePenalty = totalPenalty / 12.0;
  for (it = 0; it < nIts; it++) {
    sourcePos = Random(seq.length - length - 1);
    CopySeq(seq, orig, sourcePos, length);
    CopySeq(seq, mutated, sourcePos, length);
    conservedScore = 0;
    for (p = 0; p < orig.length; p++ ) {
      conservedScore += (ssize_t) score.scoreMat[nuc_index[orig.seq[p]]][nuc_index[orig.seq[p]]];
    }
    MutateSequence(orig, mutated, pctDivergence);
    ssize_t *locations = NULL;
    mutatedScore = (ssize_t) AffineAlign(orig, mutated,
				     0.0,0.0,0.0,
				     score.gapOpen, score.gapExtend, 
				     locations, score.scoreMat);

    expectedScore = (ssize_t) (conservedScore + 
			   averagePenalty * (pctDivergence*orig.length));
    std::cout << sourcePos << "\t" 
	      << conservedScore << "\t" 
	      << mutatedScore << "\t" 
	      << expectedScore << std::endl;
  }
}

void MutateSequence(DNASequence &original,
		    DNASequence &mutated,
		    double pctDivergence) {
  ssize_t nMutations;
  nMutations = (ssize_t) floor(pctDivergence * original.length);
  std::set<ssize_t> mutatedPositions;
  ssize_t pos;
  while (mutatedPositions.size() < nMutations) {
    pos = Random(original.length);
    mutatedPositions.insert(pos);
  }
  std::set<ssize_t>::iterator posIt, posEnd;
  posEnd = mutatedPositions.end();
  posIt = mutatedPositions.begin();
  for (; posIt != posEnd; ++posIt) {
    mutated.seq[*posIt] = OtherNuc(original.seq[*posIt]);
  }
}

ssize_t Random(ssize_t max) {
  return (ssize_t) floor((double(random())/RANDOM_MAX)*max);
}

char OtherNuc(char old) {
  ssize_t r = Random(3);
  ssize_t ni = nuc_index[old];
  ssize_t offset = ni*3+r;
  char other = other_nucs[offset + r];
  return other;
}

void CopySeq(DNASequence &src, DNASequence &dest, ssize_t pos, ssize_t length) {

  ssize_t p;
  assert (pos + length < src.length);
  assert (length <= dest.length);
  
  ssize_t end = pos + length;
  ssize_t dp = 0;
  for (p = pos; p < end; p++, dp++) {
    dest.seq[dp] = src.seq[p];
  }
}

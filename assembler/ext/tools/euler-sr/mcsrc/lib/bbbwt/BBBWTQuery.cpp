/***************************************************************************
 * Title:          BBBWTQuery.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "BBBWTQuery.h"
#include "utils.h"
#include <string>

ssize_t BW::Store(SimpleSequence &seq, BBBWT &bbb, ssize_t start, ssize_t length) {
  // Store 'N' to start a new string
  bbb.append(1).append(1).append(1).append(1).append(0);
  
  ssize_t i;
  ssize_t end = start + length;
   for (i = start; i < end; i++ ) {
    switch(seq.seq[i]) {
    case 't': case 'T': case 3 : case 255 : bbb.append(1); 
    case 'g': case 'G': case 0 : case 252 : bbb.append(1); 
    case 'c': case 'C': case 2 : case 254 : bbb.append(1); 
    case 'a': case 'A': case 1 : case 253 : bbb.append(0); 
      break;
    default: // treat like an X sign
      bbb.append(1).append(1).append(1).append(1).append(0);
      break;
    }
  }
  // Store 'N' to end a new string
  bbb.append(1).append(1).append(1).append(1).append(0);

  return 1;
}

ssize_t BW::Store(SimpleSequence &seq, BBBWT &bbb) {
  return Store(seq, bbb, 0, seq.length);
}

ssize_t BW::Query(SimpleSequence &seq, BBBWT &bbb, ssize_t &low, ssize_t &high) {
  return Query(seq, 0, seq.length, bbb, low, high);
}
    
ssize_t BW::Query(SimpleSequence &seq, ssize_t start, ssize_t length, BBBWT &bbb, ssize_t &low, ssize_t &high) {
  word_t lo=bbb.forward(bbb.begin(),0);
  word_t hi=bbb.forward(bbb.end(),0);
    // now let's do some queries
  ssize_t i = 0;
  word_t F=lo,L=hi;
  for(unsigned char *p=&seq.seq[start]; i < length && F<L; p++, ++i) {
    switch(*p) {
    case 't': case 'T': case 3 : case 255 : F = bbb.forward(F,1); L = bbb.forward(L,1);
      if (L<=F) break;
    case 'g': case 'G': case 0 : case 252 : F = bbb.forward(F,1); L = bbb.forward(L,1);
      if (L<=F) break;
    case 'c': case 'C': case 2 : case 254 : F = bbb.forward(F,1); L = bbb.forward(L,1);
      if (L<=F) break;
    case 'a': case 'A': case 1 : case 253 : F = bbb.forward(F,0); L = bbb.forward(L,0);
      break;
    default:
      L = F;
      break;
    }
  }
  low = F;
  high = L;
  return high - low;
}

void BW::Write(BBBWT &bbb, std::string &outFileName) {
  FILE *out;
  out = fopen(outFileName.c_str(), "w");
  if (! bbb.fprint(out) ) {
    std::cout << "writing of bbb failed " << std::endl;
    exit(1);
  }
}

void BW::Read(std::string &inFileName, BBBWT &bbb) {
  FILE *in;
  in = fopen(inFileName.c_str(), "r");
  if (!in) {
    std::cout << "could not open " << inFileName << std::endl;
    exit(1);
  }
  if (! bbb.fscan(in)) {
    std::cout << "could not read " << inFileName << std::endl;
    exit(1);
  }
}
  

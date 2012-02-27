/***************************************************************************
 * Title:          transseq.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/25/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "compatibility.h"

int main(int argc, char* argv[]) {

  if (argc != 2) {
    printf("usage: transseq sequence\n");
    exit(-1);
  }
  char nucs[400];
  nucs['g'] = 0;
  nucs['a'] = 1;
  nucs['c'] = 2;
  nucs['t'] = 3;
  nucs['G'] = 0;
  nucs['A'] = 1;
  nucs['C'] = 2;
  nucs['T'] = 3;

  ssize_t len = strlen(argv[1]);
  ssize_t i;
  LongTuple result;
  for (i = 0; i < len-1; i++) {
    result += nucs[argv[1][i]];
    result <<=2;
  }
  result+= nucs[argv[1][i]];

  printf(PRId_longTuple "\n", result);
  
  return 0;
}
  
  

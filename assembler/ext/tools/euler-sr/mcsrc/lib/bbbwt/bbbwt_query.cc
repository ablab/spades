/*
** Copyright 2004 Ross Lippert
** 
** This file is part of bbbwt.
** 
** bbbwt is free software; you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation; either version 2 of the License, or
** (at your option) any later version.
** 
** bbbwt is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with bbbwt; if not, write to the Free Software
** Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
**
*/
/**
   bbbwt_query - a very generic querying program
**/

#include <sys/stat.h>
#include <stdlib.h>
#include <string.h>
#include <cstdio>
#include <iostream>

using namespace std;

#include "BBBWT.h"
#include "bitdictionary.h"
#include "BWdictionary_static.h"

typedef bitdictionary_T<4>                bitdict;
typedef BW::dictionary_static_T<bitdict>  dictionary;
typedef BW::bbbwt_T<dictionary>           bbbwt;

const int bufsize=10000;   // I think 10000-mers is as high as we will go
char buf[bufsize+1];

int main(int argc,char **argv) {
  if (argc<2) {
    cerr << "bbbwt_query {bbbwt} [query]" << endl;
    exit(2);
  }
  char *s=argv[1];

  FILE *f = fopen(s,"r");

  if (!f) {
    cerr << "Unable to open file " << s << endl;
    exit(3);
  }

  bbbwt bbb;
  if (!bbb.fscan(f)) {
    cerr << "Unable to parse input binary file" << endl;
    exit(13);
  }
  fclose(f);

  const word_t lo=bbb.forward(bbb.begin(),0),hi=bbb.forward(bbb.end(),0);
  for(;;) {
    if (argc>2)
      strcpy(buf,argv[2]);
    else {
      if (scanf("%s",buf)<0) break;
    }

    // now let's do some queries
    word_t F=lo,L=hi;
    for(char *p=buf; *p && F<L; ++p) {
      switch(*p) {
      case 't': case 'T': F = bbb.forward(F,1); L = bbb.forward(L,1);
	if (L<=F) break;
      case 'g': case 'G': F = bbb.forward(F,1); L = bbb.forward(L,1);
	if (L<=F) break;
      case 'c': case 'C': F = bbb.forward(F,1); L = bbb.forward(L,1);
	if (L<=F) break;
      case 'a': case 'A': F = bbb.forward(F,0); L = bbb.forward(L,0);
	break;
      default:
	L = F;
	break;
      }
    }
    printf("%s #" word_FMT "\n",buf,L-F); // return the count of the word

    if (argc>2) break;
  }

}



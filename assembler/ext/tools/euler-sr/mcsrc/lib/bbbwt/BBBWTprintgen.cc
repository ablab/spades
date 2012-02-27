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
#include <stdlib.h>
#include <string.h>
#include <cstdio>
#include <iostream>
using namespace std;

#include "nathans_time_hacks.h"

#include "bitdictionary.h"
#include "BWdictionary_static.h"
#include "BBBWT.h"

typedef bitdictionary_T<2>                bitdict;
typedef BW::dictionary_static_T<bitdict>  dictionary;
typedef BW::bbbwt_T<dictionary>           bbbwt;

static unsigned VERBOSE=0;
static char *oname = NULL;
static char *fname = NULL;
static char *whoami;

// print modes
static bool Prefixes     =false;
static bool InverseString=false;
static bool BWTAscii     =false;
static bool BBBWTAscii   =false;

void usage(FILE *f) {
  fprintf(f,
	  "%s [-h] [-v] [-S|-P|-A|-N] [-i inputfile] [-o outputfile]\n",
	  whoami);
}

static int get_args(int argc, char **argv) {
  int errflg = 0;
  whoami = strdup(*argv);
  argv++;
#define OPT_ARG ( (argv[0][1])?(++(*argv)):(*(++argv)) )
  while(!errflg && *argv) {
    if (**argv == '-') {
      (*argv)++;
      while (!errflg && **argv) {
        switch (**argv) {
        case 'v':
          VERBOSE++;
          goto loopin;
	case 'o':
	  oname = strdup(OPT_ARG);
	  goto loopout;
	case 'i':
	  fname = strdup(OPT_ARG);
	  goto loopout;
	case 'P':
	  Prefixes = true;
	  goto loopin;
	case 'S':
	  InverseString = true;
	  goto loopin;
	case 'N':
	  BWTAscii = true;
	  goto loopin;
	case 'A':
	  BBBWTAscii = true;
	  goto loopin;
        case 'h':
          usage(stdout);
          exit(0);
          goto loopin;
        default:
          cerr << whoami << ": unknown flag '-" << *argv <<"'"<< endl;
          errflg++;
        }
      loopin:
        (*argv)++;
      }
    loopout:
      argv++;
    }
    else {
      errflg++;
    }
  }
  return errflg;
}




int main(int argc,char **argv) {
  if (get_args(argc,argv)) {
    usage(stderr);
    exit(1);
  }

  TIC

  FILE *f = stdin;
  FILE *o = stdout;
  FILE *verb = stderr;

  if (fname) {
    f = fopen(fname,"r");
    if (!f) {
      cerr << "Unable to open file " << fname << endl;
      exit(3);
    }
  }
  if (oname) {
    o = fopen(oname,"w");
    if (!o) {
      cerr << "Unable to open file " << oname << endl;
      exit(4);
    }
    verb = stdout;
  }
  
  if (VERBOSE) {
    TOC(verb) fprintf(verb," loading BBBWT\n");
    fflush(verb);
  }

  bbbwt bbb;
  if (!bbb.fscan(f)) {
    cerr << "Unable to parse input binary file" << endl;
    exit(13);
  }
  fclose(f);

  if (VERBOSE) {
    TOC(verb) fprintf(verb," bbbwt loaded\n");
    fflush(verb);
  }

  const word_t lo=bbb.forward(bbb.begin(),0);
  const word_t hi=bbb.forward(bbb.end(),0);
  if (BBBWTAscii) {  
    for(word_t i=bbb.begin(); i < bbb.end(); ++i) {
      word_t j = bbb.forward(i);

      if (j < lo)
	fputc('$',o);
      else if (j < hi)
	fputc('0',o);
      else
	fputc('1',o);
    }
  } 
  else if (BWTAscii) {
    for(word_t i=lo; i < hi; ++i) {
      unsigned cnt=0;
      word_t j=i;

      for(j = bbb.forward(j); hi <= j; j = bbb.forward(j)) ++cnt;
      if (j < lo) {
	fputc('$',o);
	j = bbb.forward(j); // this skips the padding 0 bit
      }
      else {
	switch(cnt) {
	case 0: fputc('A',o); break;
	case 1: fputc('C',o); break;
	case 2: fputc('G',o); break;
	case 3: fputc('T',o); break;
	case 4: fputc('N',o); break;
	default:
	  fprintf(verb,"internal problem: run of 0's seen as long as %u\n",cnt);
	}
      }
    }
  }
  else if (Prefixes) {
    for(word_t i=lo; i < hi; ++i) {
      word_t j=i;

      do {
	unsigned cnt=0;

	for(j=bbb.backward(j); hi <= j; j = bbb.backward(j)) ++cnt;
	if (j<lo) {
	  fputc('$',o);
	  j = bbb.backward(j); // this skips the padding 0 bit
	}
	else {
	  switch(cnt) {
	  case 0: fputc('A',o); break;
	  case 1: fputc('C',o); break;
	  case 2: fputc('G',o); break;
	  case 3: fputc('T',o); break;
	  case 4: fputc('N',o); break;
	  default:
	    fprintf(verb,"internal problem: run of 0's seen as long as %u\n",cnt);
	  }
	}
      }while (j!=i);
      fputc('\n',o);
    }
  }  
  else if (InverseString) {
    unsigned cnt=0;
    // skip initial 0 bit with two forwards
    for(word_t i=bbb.forward(bbb.forward(0)); i; i = bbb.forward(i)) {
      if (lo <= i && i < hi) {
	switch(cnt) {
	case 0: fputc('A',o); break;
	case 1: fputc('C',o); break;
	case 2: fputc('G',o); break;
	case 3: fputc('T',o); break;
	case 4: fputc('N',o); break;
	default:
	  fprintf(verb,"internal problem: run of 0's seen as long as %u\n",cnt);
	}
	cnt = 0;
      }
      else cnt++;
    }
  }
  fclose(o);

  if (VERBOSE) {
    TOC(verb) fprintf(verb," done\n");
    fflush(verb);
  }
  
}



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
#include <cstdio>
#include <iostream>
using namespace std;

#include "nathans_time_hacks.h"
#include <stdlib.h>
#include <string.h>
#include "BBBWT.h"
#include "BWdictionary_static.h"
#include "bitdictionary.h"
#include "Progress.h"

typedef bitdictionary_T<2>                bitdict;
typedef BW::dictionary_static_T<bitdict>  dictionary;
typedef BW::bbbwt_T<dictionary>           bbbwt;

static unsigned VERBOSE=0;
static bool Prog=false;
static char *oname = NULL;
static char *iname = NULL;
static char *bname = NULL;
static char *whoami;

void usage(FILE *f) {
  fprintf(f,
	  "%s "
	  "[-v] [-h] [-C] [-o outputfile] -b bwt [-i indices]\n"
	  "-h\tthis message\n"
	  "-v\tincrease verbosity (0)\n"
	  "-b\tbbbwt file\n"
	  "-i\tinput file of indices\n"
	  "-o\toutput file\n"
	  "-C\tshow a progress count on stderr (false)\n"
	  ,whoami
	  );
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
	case 'C':
	  Prog = true;
	  goto loopin;
	case 'o':
	  oname = strdup(OPT_ARG);
	  goto loopout;
	case 'b':
	  bname = strdup(OPT_ARG);
	  goto loopout;
	case 'i':
	  iname = strdup(OPT_ARG);
	  goto loopout;
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
  if (!bname) ++errflg;
  return errflg;
}


int main(int argc,char **argv) {
  if (get_args(argc,argv)) {
    usage(stderr);
    exit(1);
  }
  
  TIC
    ;

  FILE *f = stdin;
  FILE *o = stdout;
  FILE *b = fopen(bname,"r");
  FILE *verb = stderr;

  if (iname) {
    if (!(f=fopen(iname,"r"))) {
      cerr << "Unable to open file " << iname << endl;
      exit(3);
    }
  }
  if (oname) {
    if (!(o=fopen(oname,"w"))) {
      cerr << "Unable to open file " << oname << endl;
      exit(4);
    }
    verb = stdout;
  }
  if (!b) {
    cerr << "Unable to open file " << bname << endl;
    exit(5);
  }

  bbbwt bbb;
  if (!bbb.fscan(b)) {
    cerr << "Unable to parse input bbbwt file" << endl;
    exit(13);
  }
  fclose(b);

  if (VERBOSE) {
    TOC(verb) fprintf(verb," bbbwt loaded\n");
    fflush(verb);
  }
  word_t Oend = bbb.forward(bbb.end(),0);

  std::vector<word_t> BitMap;

  BitMap.resize( (Oend>>6) +1);
  for(word_t p=0; p < BitMap.size(); ++p) BitMap[p] = 0;

  const word_t MASK6=0x3f;
  const word_t ONE  =1;
  word_t I;
  while(fread(&I,sizeof(word_t),1,f)) {
    if (Oend <= I) {
      cerr << "Input error: I(" << I << ") >= Oend(" << Oend << ")" << endl;
      exit(49);
    }
    BitMap[I >> 6] |= ONE <<(I & MASK6);
  }
  if (iname) fclose(f);
  
  if (VERBOSE) {
    TOC(verb) fprintf(verb," bitmap initialized\n");
    fflush(verb);
  }

  word_t nchars=0; 
  word_t X=bbb.begin();
  Progress<word_t> prog;
  if (Prog) prog.file = stderr;
  prog.maximum = Oend;
  while((X=bbb.forward(X))!=bbb.begin()) {
    if (X < Oend) {
      if ((BitMap[X>>6]>>(X & MASK6)) & ONE) {
	prog.set(nchars);
	if (!fwrite(&X,sizeof(word_t),1,o) ||
	    !fwrite(&nchars,sizeof(word_t),1,o) ) {
	  cerr << "positions: fwrite error" << endl;
	  exit(55);
	}
      }
      ++nchars;
    }
  }

  if (VERBOSE) {
    TOC(verb) fprintf(verb," done\n");
    fflush(verb);
  }

  if (oname) fclose(o);
}

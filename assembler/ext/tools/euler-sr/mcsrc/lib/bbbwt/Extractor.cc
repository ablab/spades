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
#include <stdlib.h>
#include <string.h>

#include "BBBWT.h"

static unsigned VERBOSE=0;
static char *oname = NULL;
static char *iname = NULL;
static word_t StartOff = 0;
static word_t PrintCount = 1;
static word_t SkipCount = 0;
static char *whoami;

void usage(FILE *f) {
  fprintf(f,
	  "%s %s %s\n",whoami,
	  "[-v] [-h] [-i input] [-o output]",
	  "[-b start_offset] [-p print_count] [-s skip_count]"
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
	case 'o':
	  oname = strdup(OPT_ARG);
	  goto loopout;
	case 'i':
	  iname = strdup(OPT_ARG);
	  goto loopout;
	case 'b':
	  StartOff = atoi(OPT_ARG);
	  goto loopout;
	case 'p':
	  PrintCount = atoi(OPT_ARG);
	  goto loopout;
	case 's':
	  SkipCount = atoi(OPT_ARG);
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
  return errflg;
}


int main(int argc,char **argv) {
  if (get_args(argc,argv)) {
    usage(stderr);
    exit(1);
  }

  FILE *f = stdin;
  FILE *o = stdout;
  FILE *verbout = stderr;

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
    verbout = stdout;
  }

  word_t buf;
  for(word_t i=0; i < StartOff; ++i)
    if (!fread(&buf,sizeof(word_t),1,f)) {
      cerr << "Premature end of input while finding start offset" << endl;
      exit(10);
    }
  bool done=false;
  for(;!done;) {
    if (!done)
      for(word_t i=0; i < PrintCount; ++i) {
	if (!fread(&buf,sizeof(word_t),1,f)) done = true;
	else if (!fwrite(&buf,sizeof(word_t),1,o)) {
	  cerr << "Unable to write to output" << endl;
	  exit(12);
	}
      }
    if (!done)
      for(word_t i=0; i < SkipCount; ++i) {
	if (!fread(&buf,sizeof(word_t),1,f)) done = true;
      }
  }
  if (iname) fclose(f);
  if (oname) fclose(o);
}

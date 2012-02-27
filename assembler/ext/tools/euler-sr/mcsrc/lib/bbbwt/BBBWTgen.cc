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
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
using namespace std;

#include "nathans_time_hacks.h"

#include "bitdictionary.h"
#include "BBBWT.h"
#include "BWdictionary.h"
#include "Progress.h"

typedef BW::bbbwt_T<BW::dictionary_T<bitdictionary_T<48> > >  bbbwt;

static unsigned VERBOSE=0;
static char *oname = NULL;
static char *fname = NULL;
static char *whoami;
static bool Fasta=false;
static bool Prog=false;

void usage(FILE *f) {
  fprintf(f,
	  "%s [-h] [-v] [-C] [-F] [-i inputfile] [-o outputfile]\n"
	  "-h\tthis message\n"
	  "-v\tincreased verbosity (0)\n"
	  "-C\tshow a progress count on stderr (false) [incompatible with input from stdin]\n"
	  "-F\tinput file is fasta format (false)\n"
	  "-i\tinput file (stdin)\n"
	  "-o\toutput file (stdout)\n"
	  ,
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
	case 'F':
	  Fasta = true;
	  goto loopin;
	case 'C':
	  Prog = true;
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
  if (Prog && !fname) errflg++;
  return errflg;
}




int main(int argc,char **argv) {
  if (get_args(argc,argv)) {
    usage(stderr);
    exit(1);
  }
  
  TIC

  FILE *f=stdin;
  FILE *o=stdout;
  FILE *verb=stderr;

  if (oname) {
    o = fopen(oname,"w");
    if (!o) {
      cerr << "Unable to open file " << oname << endl;
      exit(4);
    }
    verb = stdout;
  }
  Progress<off_t> prog;
  if (fname) {
    struct stat st;
    if (stat(fname,&st) < 0) {
      cerr << "Unable to stat file " << fname << endl;
      exit(2);
    }
    prog.file = stderr;
    prog.maximum = st.st_size;
    f = fopen(fname,"r");
    if (!f) {
      cerr << "Unable to open file " << fname << endl;
      exit(3);
    }
  }
  
  if (VERBOSE) {
    TOC(verb) fprintf(verb," loading sequence\n");
    fflush(verb);
  }

  bbbwt bbb;

  // various kinds of test conditions
#if (0)
  if (!bbb.fscan(f)) {// scan in self for copy test
    cerr << "Unable to parse input binary file" << endl;
    exit(13);
  }
#else
  // read sequence from file
  bbb.append(0);

  prog.reset();
  if (Fasta) {
    // remove fasta-style headers
    for(int c; (c = fgetc(f))!= EOF; prog.add(1) ) {
      if (c=='>') {
	for( ; (c=fgetc(f))!='\n' ; prog.add(1) );
	bbb.append(1).append(1).append(1).append(1).append(0);
      }
      if (c=='\n') continue;
      switch(c) {
      case 't': case 'T': bbb.append(1);
      case 'g': case 'G': bbb.append(1);
      case 'c': case 'C': bbb.append(1);
      case 'a': case 'A': bbb.append(0);
	break;
      default: // treat like an X sign
	bbb.append(1).append(1).append(1).append(1).append(0);
	break;
      }
    }
    // terminate with an X
    bbb.append(1).append(1).append(1).append(1).append(0);
  }
  else {
    // read sequence raw  
    for(int c; (c = fgetc(f))!= EOF; prog.add(1) ) {
      switch(c) {
      case 't': case 'T': bbb.append(1);
      case 'g': case 'G': bbb.append(1);
      case 'c': case 'C': bbb.append(1);
      case 'a': case 'A': bbb.append(0);
	break;
      case '>': // let's give a warning that it looks like fasta
	cerr << whoami << ": input contains '>' characters.";
	cerr << "  Use the -F flag to process fasta files" << endl;
      default: // treat like an X sign
	bbb.append(1).append(1).append(1).append(1).append(0);
	break;
      }
    }
  }
  // reads done
#endif

  fclose(f);
  if (VERBOSE) {
    TOC(verb) fprintf(verb," bbbwt formed\n");
    fflush(verb);
  }

  if (!bbb.fprint(o)) {
    fprintf(stderr,"BBBWT::fprint failed\n");
    exit(77);
  }
  fclose(o);

#if (0)
  // reread to do some testing...
  //putchar('\n'), bbb.fdump(stdout), putchar('\n');
  bbb.clear();
  o = fopen(oname,"r");
  if (!o) {
    cerr << "trouble re-openning " << oname << endl;
    exit(16);
  }
  if (!bbb.fscan(o)) {
    fprintf(stderr,"BBBWT::fscan failed\n");
    exit(78);
  }
  if (VERBOSE) {
    TOC(verb) fprintf(verb," bbbwt reread\n");
    fflush(verb);
  }
  fclose(o);
  //putchar('\n'), bbb.fdump(stdout), putchar('\n');
#endif

  if (VERBOSE) {
    TOC(verb) fprintf(verb," done\n");
    fflush(verb);
  }
  
}



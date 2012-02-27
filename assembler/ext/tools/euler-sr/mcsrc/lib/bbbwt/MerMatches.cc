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
#include <string.h>
#include "nathans_time_hacks.h"

#include "bitdictionary.h"
#include "BWdictionary_static.h"
#include "BBBWT.h"
#include "Mers.h"
#include "Mers2.h"

typedef bitdictionary_T<4>                bitdict;
typedef BW::dictionary_static_T<bitdict>  dictionary;
typedef BW::bbbwt_T<dictionary>           bbbwt;
typedef BW::Mers_T<bbbwt>       Mers;
typedef BW::Mers2_T<bbbwt>      Mers2;

static unsigned VERBOSE=0;
static unsigned MerSize=0;
static word_t  AlphabetSize=4; // default to N = 1111
static word_t MinCount =1;
static word_t MaxCount =word_LIT(0xffffffffffffffff);
static word_t MinCount2=1;
static word_t MaxCount2=word_LIT(0xffffffffffffffff);
static bool Transitive=false;
static bool DoingPairs=false;
static bool Prog=false;
static char *oname = NULL;
static char *iname = NULL;
static char *iname2 = NULL;
static char *whoami;

void usage(FILE *f) {
  fprintf(f,
	  "%s "
	  "[-v] [-h] [-T] [-a int] [-C] "
          "-m mersize -o outputfile "
	  "[-n mincount] [-x maxcount] -i inputbwt "
	  "[[-N mincount2] [-X maxcount2] -I inputbwt2]\n"
	  "-h\tthis message\n"
	  "-v\tincrease verbosity (0)\n"
	  "-i\tinput bbbwt file\n"
	  "-I\tsecond input bbbwt file (default to matching within one file)\n"
	  "-o\toutput matches file\n"
	  "-m\tlength of pattern to match\n"
	  "-C\tshow a progress count on stderr (false)\n"
	  "-T\toutput all matches instead of a transitive subset (false) [you probably want this]\n"
	  "-n\tminimum count of a matched pattern in the first file (1)\n"
	  "-N\tminimum count of a matched pattern in the second file (1)\n"
	  "-x\tmaximum count of a matched pattern in the first file (INF)\n"
	  "-X\tmaximum count of a matched pattern in the second file (INF)\n"
	  "-a\talphabetsize (4) [don't touch this!]\n"
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
	case 'i':
	  iname = strdup(OPT_ARG);
	  goto loopout;
	case 'a':
	  AlphabetSize = atoi(OPT_ARG);
	  goto loopout;
	case 'm':
	  MerSize = atoi(OPT_ARG);
	  goto loopout;
	case 'n':
	  MinCount = atoi(OPT_ARG);
	  goto loopout;
	case 'x':
	  MaxCount = atoi(OPT_ARG);
	  goto loopout;
	case 'T':
	  Transitive = true;
	  goto loopin;
	case 'I':
	  DoingPairs = true;
	  iname2 = strdup(OPT_ARG);
	  goto loopout;
	case 'N':
	  DoingPairs = true;
	  MinCount2 = atoi(OPT_ARG);
	  goto loopout;
	case 'X':
	  DoingPairs = true;
	  MaxCount2 = atoi(OPT_ARG);
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
  if (!MerSize) ++errflg;
  if (!oname) ++errflg;
  if (!iname) ++errflg;
  if (DoingPairs && !iname2) ++errflg;
  return errflg;
}


static void critique_params(word_t Min,word_t Max) {
  if (Min > Max) {
    fprintf(stderr,"Warning: min count " word_FMT " is greater than max count " word_FMT "\n",
	    Min,Max);
  }
  if (Min == 1 && !DoingPairs) {
    fprintf(stderr,"Warning: min count " word_FMT " when not doing pairs will generate a lot of output\n",Min);
  }
  if (Min < 1) {
    fprintf(stderr,"Warning: min count " word_FMT " is nonsensical, but maybe you know what you're doing.\n",Min);
  }
}


int main(int argc,char **argv) {
  if (get_args(argc,argv)) {
    usage(stderr);
    exit(1);
  }
  
  TIC
    ;


  critique_params(MinCount,MaxCount);
  if (DoingPairs) 
    critique_params(MinCount2,MaxCount2);

  FILE *f = fopen(iname,"r");
  FILE *f2= NULL;
  FILE *o = fopen(oname,"w");

  if (!f) {
    cerr << "Unable to open file " << iname << endl;
    exit(3);
  }
  if (!o) {
    cerr << "Unable to open file " << oname << endl;
    exit(4);
  }
  if (DoingPairs) {
    f2= fopen(iname2,"r");

    if (!f2) {
      cerr << "Unable to open file " << iname2 << endl;
      exit(3);
    }
  }

  if (VERBOSE) {
    TOC(stdout) fprintf(stdout," loading bwt\n");
    fflush(stdout);
  }

  bbbwt bbb;
  bbbwt bbb2; // this one will sit around empty if not doing pairs

  if (!bbb.fscan(f)) {
    cerr << "Unable to parse input binary file" << endl;
    exit(13);
  }

  if (DoingPairs) {
    if (!bbb2.fscan(f2)) {
      cerr << "Unable to parse second input binary file" << endl;
      exit(13);
    }
  }


  if (VERBOSE) {
    if (!DoingPairs) {
      TOC(stdout) fprintf(stdout," bbbwt loaded\n");
      fflush(stdout);
    }
    else {
      TOC(stdout) fprintf(stdout," bbbwt's loaded\n");
      fflush(stdout);
    }
  }

  fclose(f);
  if (DoingPairs) fclose(f2);

  if (!DoingPairs) {
    Mers mers(bbb,AlphabetSize,MerSize);

    if (VERBOSE) {
      TOC(stdout) fprintf(stdout," computing mers\n");
      fflush(stdout);
    }
    mers.set_transitive(Transitive);
    if (Prog) {
      mers.set_progress( stderr );
    }
    mers.compute(o,MinCount,MaxCount);
    fclose(o);
    
    if (VERBOSE) {
      TOC(stdout) fprintf(stdout," computing positions\n");
      fflush(stdout);
    }
  }
  else {
    Mers2 mers(bbb,bbb2,AlphabetSize,MerSize);

    if (VERBOSE) {
      TOC(stdout) fprintf(stdout," computing mers\n");
      fflush(stdout);
    }
    mers.set_transitive(Transitive);
    if (Prog) {
      mers.set_progress( stderr );
    }
    mers.compute(o,MinCount,MaxCount,MinCount2,MaxCount2);
    fclose(o);
  }

  if (VERBOSE) {
    TOC(stdout) fprintf(stdout," done\n");
    fflush(stdout);
  }

}

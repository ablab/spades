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
#include <sys/stat.h>
#include <sys/types.h>
#include <cstdio>
#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;
#include <string.h>
#include <stdlib.h>

#include "nathans_time_hacks.h"
#include "BBBWT.h"

static unsigned VERBOSE=0;
static unsigned MerSize=0;
static bool Ascii=0;
static char *oname = NULL;
static char *iname = NULL;
static char *pname = NULL;
static char *Pname = NULL;
static char *whoami;

void usage(FILE *f) {
  fprintf(f,
	  "%s [-h] [-v] [-m mersize] [-a] -p pos_input [-P pos2_input] -i matches [-o outputfile]\n",
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
	case 'p':
	  pname = strdup(OPT_ARG);
	  goto loopout;
	case 'P':
	  Pname = strdup(OPT_ARG);
	  goto loopout;
	case 'i':
	  iname = strdup(OPT_ARG);
	  goto loopout;
	case 'm':
	  MerSize = atoi(OPT_ARG);
	  goto loopout;
	case 'a':
	  Ascii = true;
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
  //if (!MerSize) ++errflg;   // MerSize==0 indicates we are dealing with MEMs
  //Pname indicates that we are doing pair matches
  if (!pname) ++errflg;
  if (!iname) ++errflg;
  return errflg;
}

struct idxpos {
  word_t idx;
  word_t      pos;
  inline bool operator<(const idxpos &a) const {
    return idx < a.idx || (idx == a.idx && pos < a.pos);
  }
  inline bool operator<(const word_t &a) const {
    return idx < a;
  }
};
struct match {
  word_t      A,B;
};

int main(int argc,char **argv) {
  if (get_args(argc,argv)) {
    usage(stderr);
    exit(1);
  }
  
  TIC

  struct stat st;
  if (stat(pname,&st) < 0) {
    cerr << "unable to stat file " << pname << endl;
    exit(1);
  }
  if (st.st_size % sizeof(idxpos)) {
    cerr << "file " << pname << " is not a multiple of " << sizeof(idxpos) << " bytes" << endl;
    exit(1);
  }
  word_t Npos = st.st_size/ sizeof(idxpos);
  word_t Npos2=0;
  if (Pname) {
    if (stat(Pname,&st) < 0) {
      cerr << "unable to stat file " << Pname << endl;
      exit(1);
    }
    if (st.st_size % sizeof(idxpos)) {
      cerr << "file " << Pname << " is not a multiple of " << sizeof(idxpos) << " bytes" << endl;
      exit(1);
    }
    Npos2 = st.st_size/ sizeof(idxpos);
  }

  FILE *f = fopen(iname,"r");
  FILE *o = stdout; 
  FILE *p = fopen(pname,"r");
  FILE *p2=NULL;
  FILE *verb = stderr;

  if (!f) {
    cerr << "Unable to open file " << iname << endl;
    exit(3);
  }
  if (oname) {
    o = fopen(oname,"w");
    if (!o) {
      cerr << "Unable to open file " << oname << endl;
      exit(4);
    }
    verb = stdout;
  }

  if (!p) {
    cerr << "Unable to open file " << pname << endl;
    exit(5);
  }
  if (Pname) {
    p2 = fopen(Pname,"r");
    if (!p2) {
      cerr << "Unable to open file " << Pname << endl;
      exit(12);
    }
  }
  
  if (VERBOSE) {
    TOC(verb) fprintf(verb," loading positions\n");
    fflush(verb);
  }

  vector<idxpos> I(Npos);
  vector<idxpos> I2(Npos2);
  if (fread(&I[0],sizeof(idxpos),Npos,p) < Npos) {
    cerr << "fread of position file went bad" << endl;
    exit(6);
  }
  if (Pname) {
    if (fread(&I2[0],sizeof(idxpos),Npos2,p2) < Npos2) {
      cerr << "fread of position file went bad" << endl;
      exit(6);
    }
  }

  if (VERBOSE) {
    TOC(verb) fprintf(verb," sorting positions\n");
    fflush(verb);
  }

  sort(I.begin(),I.end());
  if (Pname) sort(I2.begin(),I2.end());

  if (VERBOSE) {
    TOC(verb) fprintf(verb," streaming matches\n");
    fflush(verb);
  }

  word_t Len;
  match mat;
  vector<idxpos>::const_iterator i;
  while(fread(&mat,sizeof(match),1,f)) {

    // read match length if MerSize not given
    if (MerSize)
      Len = MerSize;
    else {
      if(!fread(&Len,sizeof(Len),1,f)) {
	cerr << "Unable to read match length from file" << endl;
	exit(23);
      }
    }
    i=lower_bound(I.begin(),I.end(),mat.A);
    if (i == I.end()) {
      cerr << "Went looking for " << mat.A << " and fell off the end" << endl; exit(9);
    }
    if (i->idx != mat.A) {
      cerr << "Went looking for " << mat.A << " but only found " << i->idx << endl; exit(10);
    }
    mat.A = i->pos - Len ;

    if (Pname) {
      i=lower_bound(I2.begin(),I2.end(),mat.B);
      if (i == I2.end()) {
	cerr << "Went looking for 2nd position " << mat.B << " and fell off the end" << endl; exit(11);
      }
      if (i->idx != mat.B) {
	cerr << "Went looking for 2nd position " << mat.B << " but only found " << i->idx << endl;
	exit(12);
      }
      mat.B = i->pos - Len ;
    }
    else {
      i=lower_bound(I.begin(),I.end(),mat.B);
      if (i == I.end()) {
	cerr << "Went looking for " << mat.B << " and fell off the end" << endl; exit(11);
      }
      if (i->idx != mat.B) {
	cerr << "Went looking for " << mat.B << " but only found " << i->idx << endl;
	exit(12);
      }
      mat.B = i->pos - Len ;
    }
    if (Ascii) {
      fprintf(o,word_FMT" "word_FMT,mat.A,mat.B);
      if (MerSize) fprintf(o,"\n");
      else fprintf(o," "word_FMT"\n",Len);
    }	
    else { 
      if (!fwrite(&mat,sizeof(mat),1,o)) {
	cerr << "Postmatches: write failure" << endl;
	exit(23);
      }
      if (!MerSize && !fwrite(&Len,sizeof(Len),1,o)) {
	cerr << "Postmatches: write failure" << endl;
	exit(23);
      }
    }
  }

  fclose(o);

  if (VERBOSE) {
    TOC(verb) fprintf(verb," done\n");
    fflush(verb);
  }
}
  

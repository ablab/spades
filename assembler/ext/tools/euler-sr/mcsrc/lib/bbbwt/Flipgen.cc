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
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <cstdio>
#include <vector>

using namespace std;

static bool Reverse=false;
static bool Complement=false;
static char *fname=NULL;
static char *oname=NULL;
static char *whoami=NULL;

static unsigned VERBOSE=0;

void usage(FILE *f) {
  fprintf(f,
	  "%s [-h] [-v] [-r] [-c] [-i inputfile] [-o outputfile]\n",
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
	case 'r':
	  Reverse = true;
	  goto loopin;
	case 'c':
	  Complement = true;
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

  FILE *f=stdin;
  FILE *o=stdout;
  if (fname) {    f = fopen(fname,"r");
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
  }

  char comp_map[256];
  for(unsigned i=0; i < 256; ++i) comp_map[i] = char(i);
  comp_map[(unsigned char)'a'] = 't';
  comp_map[(unsigned char)'t'] = 'a';
  comp_map[(unsigned char)'c'] = 'g';
  comp_map[(unsigned char)'g'] = 'c';
  comp_map[(unsigned char)'A'] = 'T';
  comp_map[(unsigned char)'T'] = 'A';
  comp_map[(unsigned char)'C'] = 'G';
  comp_map[(unsigned char)'G'] = 'C';
  vector<unsigned> linelen;
  vector<char> buf;

  int c;
  unsigned cur_len=0;
  do {
    c = fgetc(f);
    if (c==EOF || (cur_len==0 && c=='>')) {
      if (!linelen.empty()) {
	if (Reverse) {
	  char tmp;
	  for(unsigned j=0;j < buf.size()/2; ++j) {
	    tmp = buf[j];
	    buf[j] = buf[buf.size()-j-1];
	    buf[buf.size()-j-1] = tmp;
	  }
	}
	if (Complement) {
	  for(unsigned j=0; j < buf.size(); ++j)
	    buf[j] = comp_map[(unsigned char)buf[j]];
	}
	unsigned j=0;
	for(unsigned i=0; i < linelen.size(); ++i) {
	  for(unsigned k=0; k < linelen[i]; ++k) 
	    fputc(buf[j+k],o);
	  j += linelen[i];
	  fputc('\n',o);
	}
      }
      if (c=='>') {
	fputc(char(c),o);
	while( (c=fgetc(f))!='\n' ) fputc(char(c),o);
	fputc(char(c),o);
      }
      linelen.clear();
      buf.clear();
      cur_len = 0;
    }
    else if (c=='\n') {
      linelen.push_back(cur_len);
      cur_len = 0;
    }
    else {
      cur_len++;
      buf.push_back(char(c));
    }
  } while(c!=EOF);

  if (fname) fclose(f);
  if (oname) fclose(o);
}

/***************************************************************************
 * Title:          Mers.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  02/18/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
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
#ifndef BBBWT_MERS_H
#define BBBWT_MERS_H
#include <stdlib.h>
#include <cstdio>
#include <vector>
#include "word.h"
#include "Progress.h"

namespace BW {
  // for doing mers within a bbbwt's

  // templated on the implementation of BBBWT
  template<class BBBWT> class Mers_T {
    const BBBWT &bbb;
    word_t MerSize;

    // during the search phase these numbers indicate the minimum
    //  and maximum counts
    word_t MinCount,MaxCount;
    // N's give mismatches.  This holds the index of the first
    //  N prefix in bbb
    word_t Nbeg,Obeg,Oend; // N avoidance zone, and 0 end zone

    FILE *Output; // file to send output

    // if true, output it put out in space inefficient transitive form.
    // if false, a subset of the matches is produced, from which the
    //  rest are inferrable
    bool transitive;

    // having a progress counter seems like such a good idea
    Progress<word_t> prog;

  public:
    // bbb, N-length and mersize given
    //  computes the Nbeg
    Mers_T(const BBBWT &b,word_t Nlen,word_t m):
      bbb(b),MerSize(m),transitive(false),
      prog() {

      Nbeg = 0;
      for(word_t i=0; i < Nlen; ++i) Nbeg = bbb.forward(Nbeg,1);

      Obeg = bbb.forward(bbb.begin(),0);
      Oend = bbb.forward(bbb.end(),0);

      prog.maximum = Oend - Obeg;
    }
    ~Mers_T() {}

    inline void set_transitive(bool b) {transitive = b;}
    inline void set_progress(FILE *f) {
      prog.file = f;
    }

    // the head of the recursive compute.  Matches are written as
    //  index pairs into FILE*
    inline void compute(FILE *o,word_t min,word_t max) {
      MinCount = min;
      MaxCount = max;
      Output = o;

      prog.reset();
      compute_rec(Obeg,Oend,0,1);
    }
  private:
    // recursive compute: given an interval, the sum of the bits
    //  of the prefixes in the interval up to size
    //  and the current size to which they are guarrenteed to match.
    // if size-sum >= Mersize, reports the matches in this interval,
    //  else extends by both 0 and 1 and recurs
    void compute_rec(word_t beg,word_t end,
		     word_t sum, word_t size) {
      word_t beg0, end0, beg1, end1;
      
      bbb.forward(beg, beg0, beg1);
      bbb.forward(end, end0, end1);
      
      // recur on the 0 extensions, possibly reporting a Mer hit
      if (end0-beg0 >= MinCount) {
	if (size-sum < MerSize) {
	  compute_rec(beg0,end0,sum,size+1);
	}
	else {
	  if (end0-beg0 <= MaxCount)
	    process(beg0,end0,sum,size);
	  prog.add(end0-beg0);
	}
      }
      else { prog.add(end0-beg0); }

      // recur on the 1 extensions
      if (end1-beg1 >= MinCount) {
	if (beg1 < Nbeg) { // avoid N's
	  compute_rec(beg1,end1,sum+1,size+1);
	}
	else { prog.add(end1-beg1); }
      }
      else { prog.add(end1-beg1); }
    }

    // The reporting function
    void process(word_t beg,word_t end,
		 word_t sum, word_t size) {
      if (!transitive) {
	// output matches in star form
	for(word_t p=beg+1; p < end; ++p) {
	  if (!fwrite(&beg,sizeof(word_t),1,Output) ||
	      !fwrite(&p,sizeof(word_t),1,Output) ) {
	    fprintf(stderr,"Mers::process: fwrite failure\n");
	    exit(55);
	  }
	}
      }
      else {
	for(word_t p=beg; p < end; ++p) {
	  for(word_t p2=p+1; p2 < end; ++p2) {
	    if (!fwrite(&p,sizeof(word_t),1,Output) ||
		!fwrite(&p2,sizeof(word_t),1,Output) ) {
	      fprintf(stderr,"Mers::process: fwrite failure\n");
	      exit(55);
	    }
	  }
	}
      }
      /** debug **
       // write out the bits for checking
       word_t j;
       word_t     sz;
       
       for(word_t i=beg; i < end; ++i) {
       j = i;
       sz= 1+size;
       do {
       fputc( '0'+(Oend<j) , stdout);
       j = bbb.backward(j);
       } while (sz--);
       fputc('\n',stdout);
       }
       fputc('\n',stdout);
      **/
    }
  };
};
#endif

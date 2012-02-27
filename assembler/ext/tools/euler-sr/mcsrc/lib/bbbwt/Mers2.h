/***************************************************************************
 * Title:          Mers2.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
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
#ifndef BBBWT_MERS2_H
#define BBBWT_MERS2_H

#include <cstdio>
#include <vector>
#include "word.h"
#include "Progress.h"

namespace BW {
  // for doing mers between pairs of bbbwt's

  // template on bbbwt implementation
  template <class BBBWT> class Mers2_T {
    const BBBWT &bbbA, &bbbB; // the two bbbwt's being matched
    word_t MerSize;

    // min and max counts for either side to prune searches
    word_t MinCountA,MaxCountA;
    word_t MinCountB,MaxCountB;
    word_t NbegA,OendA,ObegA; // N avoidance zone for A and 0 end for A
    word_t NbegB,OendB,ObegB; // N avoidance zone for B and 0 end for B

    FILE *Output; // file to send output

    // if transitive then all matches are output
    // else a subset of matches from which the rest can be infered
    bool transitive;

    // having a progress counter seems like such a good idea
    Progress<word_t> prog;

  public:
    // bind the two bbbwt's and compute the Nbegs
    Mers2_T(const BBBWT &bA,const BBBWT &bB,word_t Nlen,word_t m):
      bbbA(bA),bbbB(bB),MerSize(m), transitive(false),
      prog() {

      NbegA = 0;
      for(word_t i=0; i < Nlen; ++i) NbegA = bbbA.forward(NbegA,1);
      NbegB = 0;
      for(word_t i=0; i < Nlen; ++i) NbegB = bbbB.forward(NbegB,1);

      ObegA = bbbA.forward(bbbA.begin(),0);
      ObegB = bbbB.forward(bbbB.begin(),0);
      OendA = bbbA.forward(bbbA.end(),0);
      OendB = bbbB.forward(bbbB.end(),0);

      prog.maximum = (OendA-ObegA)+(OendB-ObegB);
    }
    ~Mers2_T() {}

    inline void set_transitive(bool b) {transitive = b;}
    inline void set_progress(FILE *f) {
      prog.file = f;
    }

    // front end for the recursive process
    inline void compute(FILE *o,
			word_t mina,word_t maxa,word_t minb,word_t maxb) {
      MinCountA = mina;
      MaxCountA = maxa;
      MinCountB = minb;
      MaxCountB = maxb;
      Output = o;

      prog.reset();

      compute_rec(ObegA,OendA,ObegB,OendB,0,1);
    }
  private:
    // given an interval of bbbA and one of bbbB with matches
    //  guarranteed up to size with total sum of sum, will either
    //  report the matches between the intervals (if (size-sum>=Mersize))
    //  or perform both 0 and 1 extensions of the intervals and recurse
    void compute_rec(word_t begA,word_t endA,
		     word_t begB,word_t endB,
		     word_t sum, word_t size) {
      word_t beg0A, end0A, beg1A, end1A;
      word_t beg0B, end0B, beg1B, end1B;
      
      bbbA.forward(begA, beg0A, beg1A);
      bbbA.forward(endA, end0A, end1A);
      
      bbbB.forward(begB, beg0B, beg1B);
      bbbB.forward(endB, end0B, end1B);
      
      // recur on the 0 extensions, possibly reporting a Mer hit
      if (end0A-beg0A >= MinCountA && end0B-beg0B >= MinCountB) {
	if (size-sum < MerSize) {
	  compute_rec(beg0A,end0A,beg0B,end0B,sum,size+1);
	}
	else {
	  if (end0A-beg0A <= MaxCountA && end0B-beg0B <= MaxCountB)
	    process(beg0A,end0A,beg0B,end0B,sum,size);
	  prog.add(end0A-beg0A+end0B-beg0B);
	}
      }
      else {prog.add(end0A-beg0A+end0B-beg0B);}
      // recur on the 1 extensions
      if (end1A-beg1A >= MinCountA && end1B-beg1B >= MinCountB) {
	if (beg1A < NbegA && beg1B < NbegB) // avoid N's
	  compute_rec(beg1A,end1A,beg1B,end1B,sum+1,size+1);
	else {prog.add(end1A-beg1A+end1B-beg1B);}
      }
      else {prog.add(end1A-beg1A+end1B-beg1B);}
    }
    // The reporting function
    void process(word_t begA,word_t endA,
		 word_t begB,word_t endB,
		 word_t sum, word_t size) {
      if (!transitive) {
	// output matches in star form
	if (begA < endA && begB < endB) {
	  for(word_t p=begB; p < endB; ++p) {
	    if (!fwrite(&begA,sizeof(word_t),1,Output) ||
		!fwrite(&p,sizeof(word_t),1,Output) ) {
	      fprintf(stderr,"Mers2::process: fwrite failure\n");
	      exit(55);
	    }
	  }
	  for(word_t p=begA+1; p < endA; ++p) {
	    if (!fwrite(&p,sizeof(word_t),1,Output) ||
		!fwrite(&begB,sizeof(word_t),1,Output) ) {
	      fprintf(stderr,"Mers2::process: fwrite failure\n");
	      exit(55);
	    }
	  }
	}
      }
      else { // if transitive
	for(word_t p=begA; p < endA; ++p) {
	  for(word_t p2=begB; p2 < endB; ++p2) {
	    if (!fwrite(&p,sizeof(word_t),1,Output) ||
		!fwrite(&p2,sizeof(word_t),1,Output) ) {
	      fprintf(stderr,"Mers2::process: fwrite failure\n");
	      exit(55);
	    }
	  }
	}
      }
    }
  };
};
#endif

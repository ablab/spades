/***************************************************************************
 * Title:          BWdictionary_static.h 
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
#ifndef BWDICTIONARY_SLAB_H
#define BWDICTIONARY_SLAB_H

#include <cstdio>
#include "word.h"
#include "bitio.h"

namespace BW {

  // definition of dictionary, just an vnode which is
  //  responsible for its deletion
  template <class BitDict>
  class dictionary_static_T {
    word_t _sum,_size,_slab_size;
    word_t *counts;
    typename BitDict::bits *words;

    // internal print function for outputting bitstreams
    bool fprint(bitio &) const;

    // internal scan to read from bitstreams
    bool fscan(bitio &);

  public:
    dictionary_static_T():
      _sum(0),_size(0),_slab_size(0),
      counts(NULL),words(NULL) {}
    // delete all allocations
    ~dictionary_static_T() { 
      clear();
    }
    void clear() {
      if (_slab_size) delete[] counts;
      if (_slab_size) delete[] words;
      _size = 0;
      _sum = 0;
      _slab_size = 0;
    }

    bool fscan(FILE *);
    bool fprint(FILE *) const;
    void fdump(FILE *f) const;

    inline word_t size() const {return _size;}
    inline word_t sum() const {return _sum;}

    inline void lookup(word_t pos,word_t &sum,bool &b) const {
      word_t off = pos/BitDict::bit_size;
      sum=counts[off];
      BitDict::lookup(words+off*BitDict::word_size,
			    BitDict::bit_size,
			    pos - off*BitDict::bit_size,
			    sum,b);
    }
    inline void lookup(word_t pos,word_t &sum) const {
      bool b=0;
      lookup(pos,sum,b);
    }
    inline void lookup(word_t pos,bool &b) const {
      word_t off = pos/BitDict::bit_size;
      word_t sum=0;
      BitDict::lookup(words+off*BitDict::word_size,
			    BitDict::bit_size,
			    pos - off*BitDict::bit_size,
			    sum,b);
    }
    inline word_t find1(word_t s) const {
      word_t lo=0,hi=_slab_size;
      for(word_t i=(hi-lo)/2; hi-lo>1; i=(hi+lo)/2) {
	if (counts[i] < s)
	  lo = i;
	else
	  hi = i;
      }
      return lo*BitDict::bit_size +
	BitDict::find1(words+lo*BitDict::word_size,
			     BitDict::bit_size,
			     s - counts[lo]);
    }
    inline word_t find0(word_t s) const {
      const word_t bsz = BitDict::bit_size;
      word_t lo=0,hi=_slab_size;
      for(word_t i=(hi-lo)/2; hi-lo>1; i=(hi+lo)/2) {
	if ( (i*bsz-counts[i]) < s)
	  lo = i;
	else
	  hi = i;
      }
      return lo*BitDict::bit_size +
	BitDict::find0(words+lo*BitDict::word_size,
			     BitDict::bit_size,
			     s - (lo*bsz-counts[lo]));
    }
  }; // class dictionary
};

template <class BitDict>
bool BW::dictionary_static_T<BitDict>::fscan(bitio &bf) {
  word_t w = 0,remainder=_size;

  counts[0] = 0;
  for(word_t i=0;i<_slab_size;++i,w+=BitDict::word_size) {
    word_t  sz = remainder;
    word_t  sm = 0;
    if (!BitDict::fscan(words+w,bf,sz,sm)) {
      return false;
    }
    remainder -= sz;
    if (sz < BitDict::bit_size && remainder) {
      fprintf(stderr,"Weirdness #1\n");
      return false;
    }
    counts[i+1] = counts[i]+sm;
  }
  if (remainder) {
    fprintf(stderr,"Weirdness #2\n");
    return false;
  }
  return true;
}

template <class BitDict>
bool BW::dictionary_static_T<BitDict>::fscan(FILE *f) {
  word_t s ;
  word_t sm;
  if (!fread(&s,   sizeof(word_t),    1,f)) return false;
  if (!fread(&sm,  sizeof(word_t),    1,f)) return false;

  _size = s;
  _sum  = sm;

  _slab_size = ((_size-1)/BitDict::bit_size)+1;
  counts = new word_t[_slab_size+1];
  words  = new typename BitDict::bits[_slab_size*BitDict::word_size];

  bitio bf(f,"r");

  return fscan(bf);
}

template <class BitDict>
bool BW::dictionary_static_T<BitDict>::fprint(bitio &bf) const {
  word_t w=0,remainder=_size;
  for(word_t i=0; i < _slab_size; ++i,w+=BitDict::word_size) {
    if (!BitDict::fprint(words+w,bf,remainder)) {
      return false;
    }
    remainder -= BitDict::bit_size;
  }
  return true;
}

template <class BitDict>
bool BW::dictionary_static_T<BitDict>::fprint(FILE *f) const {
  word_t s = size();
  word_t sm= sum();

  if (!fwrite(&s,    sizeof(word_t),    1,f)) return false;
  if (!fwrite(&sm,   sizeof(word_t),    1,f)) return false;

  bitio bf(f,"w");

  if (!fprint(bf)) return false;

  return bf.flush();
}

template <class BitDict>
void BW::dictionary_static_T<BitDict>::fdump(FILE *f) const {
  fprintf(f,"(" word_FMT ":" word_FMT ":" word_FMT ")\n",
	  _sum,_size,_slab_size);

  word_t w=0,remainder=_size;
  for(word_t i=0; i < _slab_size; ++i,w+=BitDict::word_size) {
    BitDict::fdump(words+w,f,remainder);
    remainder -= BitDict::bit_size;
    fprintf(f,"(count=" word_FMT " ",counts[i]);
    BitDict::fdump(words+w,f,BitDict::bit_size);
    fprintf(f,")\n");
  }
  fprintf(f,"(count=" word_FMT " ",counts[_slab_size]);
};

#endif

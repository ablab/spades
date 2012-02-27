/*
** Copyright 2004 Ross Lippert and Brian Walenz
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
#ifndef BITDICTIONARY_H
#define BITDICTIONARY_H

#include <cstdio>
#include <vector>

#include "word.h"
#include "bitio.h"

//  Bri's bitpacked implementation.
//
template<word_t WordSize>
struct bitdictionary_T {

  //  I assume all over that you have 64-bit words.  It will be Real Work to
  //  change this assumption.
  //
  typedef word_t bits;

  //  _words is the number of words to use in this dictionary.  _size
  //  is the number of bits in this dictionary.
  //
  static
  const word_t  word_size    = WordSize;

  static
  const word_t  bit_size     = 64 * WordSize;

  // allocate (delete) all memory required to create a new fixed size
  // array element.  size of element is set at compile time
  //
  static
  inline bits *new_bits() {
    return new bits [word_size];
  }

  static
  inline void delete_bits(bits *b) {
    delete [] b;
  }



  // return the count of the number of bits set in a word
  //
  static
  inline
  word_t  countNumberOfSetBits64(word_t x) {
    x = ((x >>  1) & word_LIT(0x5555555555555555)) + (x & word_LIT(0x5555555555555555));
    x = ((x >>  2) & word_LIT(0x3333333333333333)) + (x & word_LIT(0x3333333333333333));          
    x = ((x >>  4) & word_LIT(0x0f0f0f0f0f0f0f0f)) + (x & word_LIT(0x0f0f0f0f0f0f0f0f));
    x = ((x >>  8) & word_LIT(0x00ff00ff00ff00ff)) + (x & word_LIT(0x00ff00ff00ff00ff));              
    x = ((x >> 16) & word_LIT(0x0000ffff0000ffff)) + (x & word_LIT(0x0000ffff0000ffff));
    x = ((x >> 32) & word_LIT(0x00000000ffffffff)) + (x & word_LIT(0x00000000ffffffff));
    return(x);
  };



  //  The size parameter is the number of bits stored in this
  //  dictionary.  It is always available from the caller, and
  //  sometimes updated by the functions.

  //  Returns the sum of the bits before 'position'.
  //
  static
  inline
  word_t count(const bits *b,
                 word_t size,
                 word_t position) {
    //  Sum all the bits in whole words before our word.
    //
    word_t x=0, s=0;
    for (; x<(position >> 6); x++)
      s += countNumberOfSetBits64(b[x]);

    //  Count the number of set bits before position
    //
    if (position & 0x3f) {
      x  = b[x] << (64 - (position & 0x3f));
      s += countNumberOfSetBits64(x);
    }

    return(s);
  };



  //  Return bit value stored at some position.
  //
  static
  inline
  void lookup(const bits *b,
	      word_t size,
	      word_t position,
	      bool &bit) {
    bit = ( (b[ position >> 6 ] >> (position & 0x3f)) & 0x1);
  };

  //  Add to sum the sum of bit values before some position
  //
  static
  inline
  void lookup(const bits *b,
	      word_t size,
	      word_t position,
	      word_t &sum) {
    sum += count(b,size,position);
  };

  //  Add in the sum of bit values before some position
  //  and bit value at position
  static
  inline
  void lookup(const bits *b,
	      word_t size,
	      word_t position,
	      word_t &sum,bool &bit) {
    sum += count(b,size,position);
    bit = ( (b[ position >> 6 ] >> (position & 0x3f)) & 0x1);
  };


  //  Find the rank'th 0 in b.  Assumes that there actually is a
  //  rank'th 0 in b.
  //
  static
  inline
  word_t find0(const bits *b,
		 word_t size,
		 word_t rank) {
    word_t word=0, posn=0, s=0;

    //  Find the word that the r'th set bit is in.
    //
    s = 64 - countNumberOfSetBits64(b[word]);
    while (rank > s) {
      rank -= s;
      posn += 64;
      word++;
      s     = 64 - countNumberOfSetBits64(b[word]);
    }

    //  we must now find the position of the rank'th 0 in the current
    //  word.  XXX  There is certainly a better way to do this.
    //
    s = ~b[word];
    rank -= s & 0x1;
    while (rank > 0) {
      posn++;
      s >>= 1;
      rank -= s & 0x1;
    }
      
    return(posn);
  };


  //  Similar to find0() above.
  //
  static
  inline
  word_t find1(const bits *b,
		 word_t size,
		 word_t rank) {
    word_t word=0, posn=0, s=0;

    //  Find the word that the rank'th set bit is in.
    //
    s = countNumberOfSetBits64(b[word]);
    while (rank > s) {
      rank -= s;
      posn += 64;
      s     = countNumberOfSetBits64(b[++word]);
    }

    //  we must now find the position of the rank'th 0 in the current
    //  word.  XXX  There is certainly a better way to do this.
    //
    s = b[word];
    rank -= s & 0x1;
    while (rank > 0) {
      posn++;
      s >>= 1;
      rank -= s & 0x1;
    }

    return(posn);
  };


  // If possible, insert a bit into position, moving current occupant
  // of position and what follows to the right, returning true.  If
  // unable to accomodate the request, change nothing and return
  // false.
  //
  // If returning true, then sum is to have added to it the sum of the
  // bits to the left of position.
  //
  static
  inline
  bool insert(bits      *b,
	      word_t   size,
	      word_t   position,
	      bool       newbit,
	      word_t  &sum) {

    if (size >= bit_size)
      return(false);

    word_t  word;
    word_t  destword;
    word_t  tmpw;

    //  Starting at the last word in the vector (that's 'size'),
    //  make space for the new bit at 'position'.
    //
    //  If the 'size'th bit is in the last word, we just shift, otherwise,
    //  we spill the last bit into the next word.

    word = size >> 6;
    destword = position >> 6;

    if (word > destword) {
      if (word < (word_size - 1)) {
        b[word+1] <<= 1;
        b[word+1]  |= b[word] >> 63;
      }

      while (word > destword) {
        b[word] <<= 1;
        b[word]  |= b[word-1] >> 63;
        word--;
      }
    }

    //  Do the final shift in the word that gets the new bit

    //  Two annoying cases here.  If we are inserting into the last
    //  bit, 'tmpw' gets no piece of b[word].  Likewise, if we are
    //  inserting into the first bit, 'b[word] (below) should be
    //  erased.  We'd like to do these without if tests.

    //  If we are inserting into the last bit, then tmpw gets no
    //  piece of b[word].  (position+1) & 0x3f will be zero in this
    //  case.  If not the last bit, then we shift the extra bits off
    //  to the right, then shift the bits we want to keep back to
    //  their new spot, opening a hole for the new bit.
    //

    //  It would be nice to get rid of this if, but we leave that for
    //  a later exercise.
    //
    tmpw      = b[word];
    tmpw    >>= position & 0x3f;
    tmpw    <<= position & 0x3f;
    tmpw    <<= 1;

    //  Easy.  Insert the new bit.
    //
    tmpw     |= ((word_t)newbit) << (position & 0x3f);

    //  Like the bit above, if we are inserting into the first bit
    //  (bit #0) we first shift all the bits out through the left,
    //  then do no shift.  If not the first bit, well, just read the
    //  last comment.
    //
    b[word] <<= 63 - (position & 0x3f);
    b[word] <<= 1;
    b[word] >>= 63 - (position & 0x3f);
    b[word] >>= 1;

    b[word]  |= tmpw;

    sum += count(b, size, position);

    return(true);
  };


  // Assume that b2 contains no elements (size = 0) and move about
  // half the bits from b1 to b2.  Store new size of b1 in size.  Store
  // sum of bits in b1 in sum.
  //
  //  b2 gets the second half of b1.  b1_old = b1 + b2.
  //
  static
  inline
  void divide(bits *b1,
	      word_t &size,
	      bits *b2,
	      word_t &sum) {
    if (word_size == 1) {
      b2[0] = b1[0];
      b2[0] >>= 32;
      b1[0] <<= 32;
      b1[0] >>= 32;
      size  = 32;
      sum   = count(b1, size, size + 1);
    } else {
      word_t  word = word_size / 2;

      //  We keep the first 'word' words in b1, and copy the rest into b2.
      for (word_t p=word; p < word_size; p++) {
	b2[p-word] = b1[p];
	b1[p]      = 0;
      }

      size = word * 64;
      sum  = count(b1, size, size + 1);
    }
  }


  // If possible, absorb some of the bits of r into b, taking them
  // from r's left and sliding the elements of r down appropriately,
  // returning true.  if true is returned then sz and szr are updated
  // to the new sizes of b and r as are the sums.
  //
  //  Shift some bits from r onto b, preserve the order given by b+r.
  //
  static
  inline
  bool merge(bits *b,
	     word_t &szb,
	     word_t &smb,
	     bits *r,
	     word_t &szr,
	     word_t &smr) {

    if ((szb == bit_size) || (szr == 0))
      return(false);

    word_t  word = szb >> 6;
    word_t  i;

    //  b is probably not word aligned initially, but r is.  If there
    //  is a half-word in b, fill it with the start of r, then shift
    //  r.  Result is that both b and r are word aligned.
    //
    if (szb & 0x3f) {
      //  'inb' the number of bits in b currently,
      //  'oub' the number of bits needed to align b to words.
      //
      word_t  inb =      (szb & 0x3f);
      word_t  oub = 64 - (szb & 0x3f);
      word_t  tw;

#if 0
      fprintf(stderr, "word=%d szb=%d szr=%d inb=%d oub=%d\n",
              (int)word, (int)szb, (int)szr, (int)inb, (int)oub);
      fprintf(stderr, "%d 0x%016llx 0x%016llx\n", (int)word, b[word], r[0]);
#endif

      //  Clear the bits in b that will receive the bits from r
      //
      b[word] <<= oub;
      b[word] >>= oub;

#if 0
      fprintf(stderr, "%d 0x%016llx 0x%016llx\n", (int)word, b[word], r[0]);
#endif

      //  Grab the first 'oub' bits from r, put them onto the end of b.
      //
      b[word]  |= r[0] << inb;

#if 0
      fprintf(stderr, "%d 0x%016llx 0x%016llx\n", (int)word, b[word], r[0]);
#endif

      //  Update sizes
      //
      if (oub > szr)
        oub = szr;

      szb += oub;
      szr -= oub;

      //  Shift all of r to the left
      //
      if (szr > 0)
        for (i=1; i<word_size; i++) {
          r[i-1] >>= oub;
          tw       = r[i];
          tw     <<= inb;
          r[i-1]  |= tw;
        }

      r[word_size-1] >>= oub;

      word++;
    }

    //  Everybody is aligned, so do word copies

    for (i=word; (szr > 0) && (i < word_size); i++) {
#if 0
      fprintf(stderr, "%d 0x%016llx 0x%016llx\n", (int)i, b[i], r[i-word]);
#endif
      b[i] = r[i-word];

      if (szr >= 64) {
        szb += 64;
        szr -= 64;
      } else {
        szb += szr;
        szr  = 0;
      }
    }

    //  Finish by shifting r

    for (i=0; i < word; i++)
      r[i] = r[word_size - word + i];

    smb = count(b, bit_size, szb);
    smr = count(r, bit_size, szr);

    return true;
  }


  //
  // some sort of io for printing and reading bits from bitstreams
  //
  static
  bool fprint(const bits *b,
	      bitio &bf,
	      word_t size)  {
    if (size > bit_size) size = bit_size;
    word_t destword = size >> 6;
    word_t posnword = size & 0x3f;
    for(word_t i=0; i < destword; ++i) {
      for(int j=0; j < 64; ++j)
	if (!bf.printbit( b[i] & (word_LIT(1) << j) )) return false;
      //if (!bf.printword(b[i])) return false;
    }
    for(word_t i=0; i < posnword; ++i) {
      if (!bf.printbit( b[destword] & (word_LIT(1) << i) )) return false;
    }
    return true;
  }

  // assume size reflects how much caller desires to be read
  // sum will be overwritten
  // in other implementations may modify size
  static
  bool fscan( bits *b,
	      bitio &bf,
	      word_t &size,
	      word_t &sum)  {
    if (size > bit_size) size = bit_size;
    word_t destword = size >> 6;
    word_t posnword = size & 0x3f;
    for(word_t i=0; i < destword; ++i) {
      if (!bf.scanword(b[i])) return false;
    }
    if (posnword) {
      b[destword] = 0;
      for(word_t i=0; i < posnword; ++i) {
	bool bol;
	if (!bf.scanbit( bol )) return false;
	b[destword] |= (word_t(bol) << i);
      }
    }
    sum = count(b, size, size);
    return true;
  }


  //
  // something useful for debugging, which prints ASCII
  //
  static
  void fdump(const bits *b,
	     FILE *file,
	     word_t size) {
    word_t x;
    word_t s = 0;
  
    for (word_t i=0; i<word_size && s<size; i++) {
      x = b[i];
      for (word_t j=0; j<64 && s<size; j++) {
	fprintf(file, "%d", (int)(x & 0x1));
	x >>= 1;
	s++;
      }
      fprintf(stdout, " ");
    }
  }

};

#endif

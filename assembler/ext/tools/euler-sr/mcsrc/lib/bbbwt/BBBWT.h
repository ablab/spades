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
#ifndef BBBWT_H
#define BBBWT_H

#include <cstdio>
#include "word.h"

namespace BW {

  // to prevent various tragedies when reading and writing binary
  // files, a magic number if affixed to the beginning of each
  // file printed by fprint and checked for by fscan
  const word_t magic_word=word_LIT(0x0bbb000000000001);

  // The class is templated on the underlying DICITONARY implementation
  // to facilitate experiments with other implementations.
  template <class DICTIONARY> class bbbwt_T {
    // DICTIONARY + dollar == bbbwt
    DICTIONARY bd;  // hold bits
    word_t   _dollar; // hold next insert position and mark start/stop
   
  public:
	bbbwt_T(): bd(), _dollar(0) {}
    ~bbbwt_T() {}
   
    // put a new bit in the bitstring represented by this.
    // inserts b into the dictionary at dollar and recomputes
    // dollar.
    bbbwt_T<DICTIONARY> &append(bool b) { 
      word_t s=0;
      bd.insert(_dollar,b,s);
      if (!b)
				_dollar = 1+_dollar-s;
      else
				_dollar = 1+size0()+s;
      return (*this);
    }
    
    // passthrough methods for freeze and thaw
    //  freeze: reclaim memory and condense things, no appends after a freeze
    inline void freeze() { bd.freeze(); }
    //  thaw: checks/restores internals so that further appends are possible 
    inline void thaw() { bd.thaw(); }
    
    
    // The basic interval of all prefixes in the bbbwt
    //  intervals should start at a prefix which marks the first
    //  occurrence of the subword and end just AFTER the prefix
    //  which is the last occurrence.
    inline word_t begin() const {return 0;}
    inline word_t end()   const {return 1+bd.size();}

    // location of the dollar in the bbbwt
    inline word_t dollar() const {return _dollar;}
    
    // counts of things
    inline word_t  size() const {return end()-begin();} // size of the string rep-ed
  private:
    inline word_t  size0() const {return bd.size()-bd.sum();} // number of 0's
    inline word_t  size1() const {return bd.sum();}           // number of 1's
    
    // bbbwt functions:

  private:
    // counts the number of 0's and 1's before position in the bbbwt
    inline void occ(word_t idx,word_t &size0,word_t &size1) const {
      const bool dol = (_dollar < idx);
      size1 = 0;
      bd.lookup(idx-dol,size1);
      size0 = idx-dol-size1;
    }

  public:
    // takes position idx to the first occurrence of the best possible match
    // to the prefix at idx if it were extended by the bit b
    inline word_t forward(word_t idx,bool b) const {
      word_t s0,s1;
      occ(idx,s0,s1);
      return 1 + (b?(size0()+s1):(s0));
    }
    // computes both "forwards" at once
    inline void forward(word_t idx,word_t &idx0,word_t &idx1) const {
      word_t s0,s1;
      occ(idx,s0,s1);
      idx0 = 1 + s0;
      idx1 = 1 + size0() + s1;
    }

    // follow prefix represented by position in the bbbwt to the
    // prefix which extends it by one character
    word_t forward(word_t idx) const {
      word_t suml=0;
      bool b;
      if (idx < dollar()) {
				bd.lookup(idx,suml,b);
				if (!b)
					return 1+idx-suml;
				else
					return 1+size0()+suml;
      }
      else if (idx > dollar()) {
				bd.lookup(idx-1,suml,b);
				if (!b)
					return 1+idx-1-suml;
				else
					return 1+size0()+suml;
      }
      else {
				return 0;
      }
    }
    // follow prefix represented by position in the bbbwt to the
    // prefix which results from truncated by a character
    word_t backward(word_t idx) const {
      if (idx == 0) {
				return dollar();
      }
      else {
				if (idx + bd.sum() <= bd.size())
					idx = bd.find0(idx);
				else
					idx = bd.find1(idx-size0());
				return idx+(idx>=dollar());
      }
    }
    
    // reset state of bbbwt to that of representing an empty string
    inline void clear() {
      _dollar=0;
      bd.clear();
    }

    // a dumper to support debugging
    void fdump(FILE *f) const {
      fprintf(f,"%d:%d:%d:",bd.size(),bd.sum(),dollar());
      bd.fdump(f);
    }
    // print state of the bbbwt to a stream.  return true on success
    bool fprint(FILE *f) const {
      word_t mnum = magic_word;
      if (!fwrite(&mnum,sizeof(word_t),1,f)) return false;
      if (!fwrite(&_dollar,sizeof(word_t),1,f)) return false;
      return bd.fprint(f);
    }
    // print state of the bbbwt from a stream.  returns true on success
    // only use with new bbbwt's or "clear"-ed ones
    bool fscan(FILE *f) {
      clear();
      word_t mnum;
      if (!fread(&mnum,sizeof(word_t),1,f)) return false;
      if (mnum != magic_word) {
				fprintf(stderr,"bbbwt::fscan : magic number mismatch\n");
				return false;
      }
      if (!fread(&_dollar,sizeof(word_t),1,f)) return false;
      return bd.fscan(f);
    }
    
    // useful tools for checking up on the tree
    inline bool rb_check() const { word_t bh=0; return bd.rb_check(bh);}
    inline word_t depth() const {return bd.depth();}
    inline word_t num_nodes() const {return bd.num_nodes();}
  }; // bbbwt

}; // namespace

#endif

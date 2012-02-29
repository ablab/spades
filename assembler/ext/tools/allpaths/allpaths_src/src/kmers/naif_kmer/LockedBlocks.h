///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef KMERS__NAIF_KMER__LOCKED_BLOCKS_H
#define KMERS__NAIF_KMER__LOCKED_BLOCKS_H

#include "system/SpinLockedData.h"


// ---------- class to manage blocks of data in a thread-safe way --------


class LockedBlocks : public vec<uint64_t>, private SpinLockedData
{

public:
  LockedBlocks(const size_t n = 0)
    : vec<uint64_t>(n, 0), SpinLockedData()
  {}

  size_t lock_some_block(const std::set<size_t> & i_blocks)
  {
    ForceAssertGt(i_blocks.size(), 0ul);
 
    SpinLocker lock(*this);  // destructor of 'lock' unlocks mutex
    
    std::set<size_t>::iterator it = i_blocks.begin();
    while (1) {

      if (!(*this)[*it]) {   // block is not locked
        (*this)[*it] = 1;
        return *it;
      }
      
      it++;
      if (it == i_blocks.end()) it = i_blocks.begin();      
    }
  }

  // we can unlock the block without a Locker if blocks 
  // are associated with a vec<uint64_t>
  // if you have a vec<bool>, alignment could be a problem 
  void unlock_block(const size_t iblk)
  {
    ForceAssert((*this)[iblk]);
    (*this)[iblk] = 0;
  }

  String str() const
  {
    String s = "";
    for (size_t i = 0; i != this->size(); i++)
      s += ((*this)[i] ? "1" : "0");
    return s;
  }

};


#endif

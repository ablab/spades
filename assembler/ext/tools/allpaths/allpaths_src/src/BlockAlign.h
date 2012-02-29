///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef BLOCK_ALIGN_H
#define BLOCK_ALIGN_H

#include "Alignment.h"
#include "Block.h"
#include "lookup/LookAlign.h"



/*
 * Class block_align
 *
 * It makes easier to deal with the aligning portions of a given align
 * (blocks). Each align is broken into its gapless units (block), with the
 * option of saving also mutation rates (on a per block basis).
 *
 * Constructor and SetFromAlign are passed an alignment, and a vector of
 * mutations (as pointers). The latter may be null, but if not it must
 * have the same size as the alignment.
 */
class block_align {
  
public:

  block_align( ) { }
  
  block_align( const align *al, const vec<int> *mut = 0 );
  
  block_align( const look_align_plus *hit );
  
  void SetFromAlign( const align *al, const vec<int> *mut = 0 );
  
  void SetFromLookAlign( const look_align_plus *hit );
  
  size_t size( ) const { return blocks_.size( ); }
  
  const block &operator[] ( int ii ) const { return blocks_[ii]; }
  
  
private:
  
  void Setup( const align *al, const vec<int> *mut = 0 );
  
  
private:
  
  vec<block> blocks_;
  
};
  


#endif

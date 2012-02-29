///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Alignment.h"
#include "BlockAlign.h"
#include "lookup/LookAlign.h"



/*
 * BlockAlign
 * Constructor
 */
block_align::block_align( const align *al, const vec<int> *mut )
{
  this->Setup( al, mut );
}



/*
 * BlockAlign
 * Constructor
 */
block_align::block_align( const look_align_plus *hit )
{
  this->Setup( &( hit->a ), &( hit->mutations_by_block ) );
}



/*
 * BlockAlign
 * SetFromAlign
 */
void block_align::SetFromAlign( const align *al, const vec<int> *mut )
{
  this->Setup( al, mut );
}



/*
 * BlockAlign
 * SetFromLookAlign
 */
void block_align::SetFromLookAlign( const look_align_plus *hit )
{
  this->Setup( &( hit->a ), &( hit->mutations_by_block ) );
}



/*
 * BlockAlign
 * Setup
 */
void block_align::Setup( const align *al, const vec<int> *mut )
{
  if ( mut )
    ForceAssert( (int)mut->size( ) == al->Nblocks( ) );
  blocks_.clear( );
  blocks_.reserve( al->Nblocks( ) );
  int p1 = al->pos1( );
  int p2 = al->pos2( );
  for (int ii=0; ii<al->Nblocks( ); ii++) {
    int mutations = ( mut ) ? (*mut)[ii] : 0;
    int length = al->Lengths( ii );
    int gap = al->Gaps( ii );
    if ( gap < 0 )
      p1 += -gap;
    if ( gap > 0 )
      p2 += gap;
    block new_block( length, mutations, p1, p2 );
    blocks_.push_back( new_block );
    p1 += length;
    p2 += length;
  }
}




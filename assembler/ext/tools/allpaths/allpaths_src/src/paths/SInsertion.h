///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef S_INSERTION_H
#define S_INSERTION_H

#include "SupersHandler.h"
#include "paths/CSublink.h"

/**
 * struct SInsertion
 *
 * Structure encapsulating an insertion of a small scaffold inside a
 * bigger one.
 */
struct SInsertion {
  
public:
  
  SInsertion( );
  
  SInsertion( int small_id, int big_id, bool small_rc, int pos_in_big,
	      pair<int,int> gap_before, pair<int,int> gap_after);
  
  SInsertion( const vec<superb> &scaffolds, const CSublink &link,
	      const int gap_id, const int small_start, const int gap_len );
  
  void PrintInfo( const vec<superb> &scaffolds, ostream &out ) const;

  friend bool operator< ( const SInsertion &left, const SInsertion &right );
  
  
public:

  int small_id_;              // id of small scaffold
  int big_id_;                // id of large scaffold
  bool small_rc_;             // if small scaffold is inserted rc in big
  int pos_in_big_;            // insert small in big after this contig
  pair<int,int> gap_before_;  // new gap before small (after insertion)
  pair<int,int> gap_after_;   // new gap after small (after insertion)

};

/**
 * SInsertion_sorter_Big_StartOnBig
 *
 * Sort by big_id_, start of insertion on big.
 */
struct SInsertion_sorter_Big_StartOnBig :
  public binary_function< const SInsertion &, const SInsertion &, bool >
{
  bool operator( ) ( const SInsertion &left, const SInsertion &right ) {
    if ( left.big_id_ < right.big_id_ ) return true;
    if ( left.big_id_ > right.big_id_ ) return false;
    if ( left.pos_in_big_ < right.pos_in_big_ ) return true;
    if ( left.pos_in_big_ > right.pos_in_big_ ) return false;
    return ( left.gap_before_.first < right.gap_before_.first );
  }
};

#endif

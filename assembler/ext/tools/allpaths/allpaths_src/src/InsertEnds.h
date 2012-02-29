// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#ifndef INSERT_ENDS_H
#define INSERT_ENDS_H

#include <iostream>
#include <set>
#include "ReadLocation.h"



/*
 * class insert_ends
 *
 * Simple class for the storage of the two ends of an insert for which
 * both ends have been assembled.
 */
class insert_ends {
  
public:
  
  insert_ends( ) :
    loc1_ ( -1 ), loc2_ ( -1 ), pair_id_ ( -1 ) { }
  
  insert_ends( int loc1, int loc2, int pair_id ) :
    loc1_ ( loc1 ), loc2_ ( loc2 ), pair_id_ ( pair_id ) { }
  
  void Set( int loc1, int loc2, int pair_id ) {
    loc1_ = loc1; loc2_ = loc2; pair_id_ = pair_id; }
  
  int Loc1( ) const { return loc1_; }
  
  int Loc2( ) const { return loc2_; }

  int PairId( ) const { return pair_id_; }

  friend ostream& operator<< ( ostream &out, const insert_ends &ins ) {
    out << ins.loc1_ << "\t"
	<< ins.loc2_ << "\t"
	<< ins.pair_id_;
    return out;
  }
  
  friend istream& operator>> ( istream &in, insert_ends &ins ) {
    in >> ins.loc1_
       >> ins.loc2_
       >> ins.pair_id_;
    return in;
  }
  
  
private:
  
  int loc1_;      // position of first end in vec of read_locations.
  int loc2_;      // position of second end in vec of read_locationss.
  int pair_id_;   // position in the vector of read_pairings.

};



/*
 * order_InsertEnds_Contigs
 * ordering functor (sort by contig_id).
 */
struct order_InsertEnds_Contigs
  : public binary_function <const insert_ends&, const insert_ends&, bool>
{
private:
  const vec<read_location> &locs_;

public:
  order_InsertEnds_Contigs( const vec<read_location> &locs ) :
    locs_( locs ) { }
  
  bool operator() ( const insert_ends &left, const insert_ends &right ) {
    int cg_left_1 = locs_[left.Loc1( )].Contig( );
    int cg_left_2 = locs_[left.Loc2( )].Contig( );
    int cg_right_1 = locs_[right.Loc1( )].Contig( );
    int cg_right_2 = locs_[right.Loc2( )].Contig( );
 
    int cg_ll = ( cg_left_1 < cg_left_2 ) ? cg_left_1 : cg_left_2;
    int cg_lr = ( cg_left_1 < cg_left_2 ) ? cg_left_2 : cg_left_1;
    int cg_rl = ( cg_right_1 < cg_right_2 ) ? cg_right_1 : cg_right_2;
    int cg_rr = ( cg_right_1 < cg_right_2 ) ? cg_right_2 : cg_right_1;

    if ( cg_ll == cg_rl && cg_lr == cg_rr ) {
      if ( left.Loc1( ) == right.Loc1( ) )
	return ( left.Loc2( ) < right.Loc2( ) );
      return ( left.Loc1( ) < right.Loc1( ) );
    }
    if ( cg_ll == cg_rl )
      return ( cg_lr < cg_rr );
    return ( cg_ll < cg_rl );
  }
};



/*
 * order_InsertEnds_Locs
 * ordering functor (sort by loc1_, loc2_).
 */
struct order_InsertEnds_Locs
  : public binary_function <const insert_ends&, const insert_ends&, bool>
{
  bool operator() ( const insert_ends &left, const insert_ends &right ) {
    if ( left.Loc1( ) == right.Loc1( ) )
      return ( left.Loc2( ) < right.Loc2( ) );
    return ( left.Loc1( ) < right.Loc1( ) );
  }
};



/*
 * order_InsertEnds_PairId
 * ordering functor (sort by pair_id).
 */
struct order_InsertEnds_PairId
  : public binary_function <const insert_ends&, const insert_ends&, bool>
{
  bool operator() ( const insert_ends &left, const insert_ends &right ) {
    if ( ( left.PairId( ) == right.PairId( ) ) &&
	 ( left.Loc1( ) == right.Loc1( ) ) )
      return ( left.Loc2( ) < right.Loc2( ) );
    if ( left.PairId( ) == right.PairId( ) )
      return ( left.Loc1( ) < right.Loc1( ) );
    return ( left.PairId( ) < right.PairId( ) );
  }
};



#endif

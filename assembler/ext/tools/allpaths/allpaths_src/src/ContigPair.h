// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#ifndef CONTIG_PAIR_H
#define CONTIG_PAIR_H

#include "InsertEnds.h"
#include "Vec.h"



/*
 * class contig_pair
 *
 * Two contigs and their links. Contig ids are always stored so that
 * cg1 < cg2.
 */
class contig_pair {
  
public:
  contig_pair (int super, int cg1, int cg2 ) : super_ ( super ) {
    cg1_ = ( cg1 < cg2 ) ? cg1 : cg2;
    cg2_ = ( cg1 < cg2 ) ? cg2 : cg1;
  }
  
  void AddInsert( const insert_ends *ins ) { vecp_ins_.push_back( ins ); }
  
  int Super( ) const { return super_; }

  int Cg1( ) const { return cg1_; }
  
  int Cg2( ) const { return cg2_; }

  const insert_ends *operator[] ( int ii ) const { return vecp_ins_[ii]; }

  unsigned int size ( ) const { return vecp_ins_.size( ); }
  
  friend bool operator== ( const contig_pair &l, const contig_pair &r ) {
    return ( l.cg1_ == r.cg1_ && l.cg2_ == r.cg2_ );
  }

  friend bool operator< ( const contig_pair &left, const contig_pair &right ) {
    if ( left.super_ == right.super_ ) {
      if ( left.cg1_ == right.cg1_ )
	return ( left.cg2_ < right.cg2_ );
      return ( left.cg1_ < right.cg1_ );
    }
    return ( left.super_ < right.super_ );
  }
  
  
private:
  int super_;                             // supercontig
  int cg1_;                               // first contig (cg1_<cg2_)
  int cg2_;                               // second contig
  vec< const insert_ends* > vecp_ins_;    // the actual inserts
  
};



/*
 * order_contig_pair_LinksCount
 * ordering functor
 *
 * Sort contig_pair's by their number of links (left<right means left has
 * less links than right).
 */
struct order_contig_pair_LinksCount
  : public binary_function<const contig_pair&, const contig_pair&, bool>
{
public:
  bool operator() ( const contig_pair &left, const contig_pair &right ) {
    return ( left.size( ) < right.size( ) );
  }
};



#endif

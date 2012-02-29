///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef C_BUNDLE_H
#define C_BUNDLE_H

#include "Vec.h"
#include "math/Functions.h"
#include "paths/reporting/COffset.h"

/**
 * class CBundle
 *
 * Container for a bundle of links between contigs or supers. It
 * stores the offset between id1_ and id2_, rather than the gap
 * between them.
 */
class CBundle {

public:

  CBundle( ) :
    id1_ ( -1 ),
    id2_ ( -1 ),
    rc2_ ( false ),
    weight_ ( 0 ),
    score_ ( 0.0 ),
    offset_ ( make_pair( 0, 0 ) ),
    handle1_ ( make_pair( 0, 0 ) ),
    handle2_ ( make_pair( 0, 0 ) ) { }
    
  CBundle( int id1,
	   int id2,
	   bool rc2,
	   int weight,
	   float score,
	   pair<int,int> offset,
	   pair<int,int> handle1,
	   pair<int,int> handle2 ) :
    id1_ ( id1 ),
    id2_ ( id2 ),
    rc2_ ( rc2 ),
    weight_ ( weight ),
    score_ ( score ),
    offset_ ( offset ),
    handle1_ ( handle1 ),
    handle2_ ( handle2 ) { }
  
  void SetId1( int id1 )                 { id1_ = id1; }
  void SetId2( int id2 )                 { id2_ = id2; }
  void SetOffset( pair<int,int> offset ) { offset_ = offset; }

  int           Id1( )     const { return id1_; }
  int           Id2( )     const { return id2_; }
  Bool          Fw2( )     const { return ! rc2_; }
  Bool          Rc2( )     const { return rc2_; }
  int           Weight( )  const { return weight_; }
  float         Score( )   const { return score_; }
  pair<int,int> Offset( )  const { return offset_; }
  pair<int,int> Handle1( ) const { return handle1_; }
  pair<int,int> Handle2( ) const { return handle2_; }
  
  void Print( int len1, int len2, ostream &out ) const {
    out << id1_ << " -> " << id2_ << ( rc2_ ? "[-]" : "[+]" ) << "   "
	<< "w=" << weight_ << "   s=" << ToString( score_, 2 ) << "   "
	<< "<" << offset_.first << "," << offset_.second << ">   "
	<< "h1=[" << handle1_.first
	<< "," << handle1_.second << ")_" << len1 << "   "
	<< "h2=[" << handle2_.first
	<< "," << handle2_.second << ")_" << len2 << "\n";
  }

  friend bool operator< ( const CBundle &left, const CBundle &right ) {
    if ( left.Id1( ) < right.Id1( ) ) return true;
    if ( left.Id1( ) > right.Id1( ) ) return false;

    if ( left.Id2( ) < right.Id2( ) ) return true;
    if ( left.Id2( ) > right.Id2( ) ) return false;

    if ( left.Rc2( ) < right.Rc2( ) ) return true;
    if ( left.Rc2( ) > right.Rc2( ) ) return false;

    if ( left.Weight( ) < right.Weight( ) ) return true;
    if ( left.Weight( ) > right.Weight( ) ) return false;
    
    return ( left.Offset( ).first > right.Offset( ).first );   // notice >
    
    // Unused: score_, handle1_, handle2_;
  }
  
  
private:
  
  int    id1_;         // id of first contig or super
  int    id2_;         // id of second contig or super
  bool   rc2_;         // orientation of id2_ wrt id1_
  
  int            weight_;   // number of links
  float          score_;    // score of bundle
  pair<int,int>  offset_;   // offset (mean/stdev)
  pair<int,int>  handle1_;  // where bundles land on id1_ (interval)
  pair<int,int>  handle2_;  // where bundles land on id1_ (interval)
  
};

#endif

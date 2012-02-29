///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef C_LINK_BUNDLE_H
#define C_LINK_BUNDLE_H

#include "CoreTools.h"

/**
 * class CLinkBundle
 *
 * A bundle of links between two (oriented) supers, characterized by
 * separation, standard deviation, weight (ie number of links) and
 * score (as defined in CombineNormalDistributions).
 */
struct CLinkBundle { 

public:
  
  CLinkBundle( );
  
  CLinkBundle( int sep, int dev, int weight, float score,
	       pair<double,double> win1,
	       pair<double,double> win2 );
  
  void Set( int sep, int dev, int weight, float score,
	    pair<double,double> win1,
	    pair<double,double> win2 );
  
  int Sep( ) const { return sep_; }
  int Dev( ) const { return dev_; }
  int Weight( ) const { return weight_; }
  float Score( ) const { return score_; }
  pair<double,double> Win1( ) const { return win1_; }
  pair<double,double> Win2( ) const { return win2_; }
  
  String AsString( bool brief = false ) const;

  // Returns a combined score of Score and Weight (higher is better, see .cc).
  double CombinedScore( ) const;

  friend bool operator< ( const CLinkBundle &left, const CLinkBundle &right ) {
    if ( left.Sep( ) != right.Sep( ) )
      return ( left.Sep( ) < right.Sep( ) );
    if ( left.Dev( ) != right.Dev( ) )
      return ( left.Dev( ) < right.Dev( ) );
    if ( left.Weight( ) != right.Weight( ) )
      return ( left.Weight( ) < right.Weight( ) );
    return ( left.Score( ) < right.Score( ) );
    // Do not bother with win1_, win2_ (should do for completeness).
  }
  
  
public:
  
  int sep_;
  int dev_;
  int weight_;
  float score_;
  pair<double,double> win1_;
  pair<double,double> win2_;
  
};

/**
 * CLinkBundle_order_best
 *
 * Order bundles based on "best bundles go first" criteria. Sort by:
 *    combined_score (higher score go first)
 *    sep
 *    dev
 */
struct CLinkBundle_order_best
  : public binary_function<const CLinkBundle&, const CLinkBundle&, bool>
{
  bool operator( ) ( const CLinkBundle &left, const CLinkBundle &right ) const {
    if ( left.CombinedScore( ) != right.CombinedScore( ) )
      return ( left.CombinedScore( ) > right.CombinedScore( ) );
    if ( left.Sep( ) != right.Sep( ) )
      return ( left.Sep( ) < right.Sep( ) );
    return ( left.Dev( ) < right.Dev( ) );
    // Do not bother with win1_, win2_ (should do for completeness).
  }
};

struct pCLinkBundle_order_best
  : public binary_function<const CLinkBundle*, const CLinkBundle*, bool>
{
  bool operator( ) ( const CLinkBundle *left, const CLinkBundle *right ) const {
    CLinkBundle_order_best sorter;
    return sorter( *left, *right );
  }
};

/**
 * BinaryWrite
 */
void BinaryWrite( int fd, const CLinkBundle &bundle );

/**
 * BinaryRead
 */
void BinaryRead( int fd, CLinkBundle &bundle );

#endif

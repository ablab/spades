///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef C_OFFSET_H
#define C_OFFSET_H

#include "PairsManager.h"
#include "String.h"
#include "VecTemplate.h"
#include "math/Functions.h"
#include "math/HoInterval.h"

/**
 * struct SLink
 *
 * Core info for a link between two supers.
 */
struct SLink { 

public:
  
  SLink( ) :
    offset_ ( 0 ),
    win1_ ( ho_interval( 0, 0 ) ),
    win2_ ( ho_interval( 0, 0 ) ),
    pair_id_ ( -1 ) { }

  SLink( const normal_distribution &offset,
	 const ho_interval &win1,
	 const ho_interval &win2,
	 const longlong &pair_id ) :
    offset_ ( offset ),
    win1_ ( win1 ),
    win2_ ( win2 ),
    pair_id_ ( pair_id ) { }
  
  bool OverlapsWith( const SLink &other, const float stretch ) const {
    float other_rad = ( stretch * other.offset_.sigma_ );
    float this_rad = ( stretch * offset_.sigma_ );
    float left = Max( offset_.mu_ - this_rad, other.offset_.mu_ - other_rad );
    float right = Min ( offset_.mu_ + this_rad, other.offset_.mu_ + other_rad );
    return ( right - left > 0 );
  }

  friend bool operator< ( const SLink &left, const SLink&right ) {
    if ( left.offset_.mu_ != right.offset_.mu_ )
      return ( left.offset_.mu_ < right.offset_.mu_ );
    if ( left.offset_.sigma_ != right.offset_.sigma_ )
      return ( left.offset_.sigma_ < right.offset_.sigma_ );
    if ( left.win1_ != right.win1_ )
      return ( left.win1_ < right.win1_ );
    if ( left.win2_ != right.win2_ )
      return ( left.win2_ < right.win2_ );
    return ( left.pair_id_ < right.pair_id_ );
  }

  friend bool operator==( SLink const& sl1, SLink const& sl2 )
  { return sl1.offset_ == sl2.offset_ &&
           sl1.win1_ == sl2.win1_ &&
           sl1.win2_ == sl2.win2_ &&
           sl1.pair_id_ == sl2.pair_id_; }

public:
  
  normal_distribution offset_;   // offset and stdev for this link
  ho_interval win1_;             // windows of read on first super
  ho_interval win2_;             // windows of reads on second super
  longlong pair_id_;             // id of pair
  
};

/**
 * class COffset
 *
 * Offset between two supers implied by a single links. The first
 * super is always oriented fw, the second may be rc.
 */
class COffset {
  
public:
  
  COffset( );
  
  COffset( int s1, int s2, bool rc2, int slen1, int slen2 );

  void SetSupers( int s1, int s2, bool rc2, int slen1, int slen2 );
  
  void AddLink( const SLink &link );
  
  // Constant accessors.
  int Super1( ) const { return super1_; }
  int Super2( ) const { return super2_ < 0 ? - super2_ - 1 : super2_ ; }
  int Slen1( ) const { return slen1_; }
  int Slen2( ) const { return slen2_; }
  int Rc2( ) const { return super2_ < 0; }
  const vec<SLink>& Links( int cluster_id ) const { return links_[cluster_id]; }
  
  size_t NClusters( ) const;
  size_t NLinksTotal( ) const;
  size_t NLinks( int cluster_id ) const;
  float Score( int cluster_id ) const;
  double CombinedScore( int cluster_id ) const;
  normal_distribution Offset ( int cluster_id ) const;
  normal_distribution Offset ( int cluster_id,
			       const PairsManager &pairs,
			       vec<normal_distribution> &lib_nds ) const;

  // Estimate the error rate in the measurement of the gap.
  int DevErrorEstimate( int cluster_id, const PairsManager &pairs ) const;

  // Windows are normalized to the supers lengths (result in percent).
  pair<double,double> SpreadWin1( int cluster_id ) const;
  pair<double,double> SpreadWin2( int cluster_id ) const;
  
  // As above, but return windows in base coordinates.
  pair<int,int> SpreadWinBases1( int cluster_id ) const;
  pair<int,int> SpreadWinBases2( int cluster_id ) const;

  // Return separation rather than offset.
  normal_distribution Separation( int cluster_id ) const;

  // Check if supers (and orientation) match those of other.
  bool MatchesSupersWith( const COffset &other ) const;    
  
  // Run validity tests on given cluster (see .cc for details).
  bool IsValid( int cluster_id,
		const size_t *MIN_LINKS = 0,
	        const int *MAX_SPREAD = 0,
		const pair<int,int> *MIN_MAX_GAP = 0 ) const;
  
  // Print info (verbose multi-lines if pairs != 0).
  void Print( ostream &out,
	      const PairsManager *pairs = 0,
	      const bool smart = true ) const;
  
  // Cluster links into consistent (sorted) sets.
  void ClusterLinks( const float stretch = 1.0 ) const;

  // operator<
  friend bool operator< ( const COffset &left, const COffset &right )
  { return left.Super1() < right.Super1() ||
           left.Super1() == right.Super1() &&
               (left.Super2() < right.Super2() ||
                left.Super2() == right.Super2() && left.Rc2() < right.Rc2()); }

  friend bool operator==( COffset const& co1, COffset const& co2 )
  { return co1.super1_ == co2.super1_ &&
           co1.super2_ == co2.super2_ &&
           co1.slen1_ == co2.slen1_ &&
           co1.slen2_ == co2.slen2_ &&
           co1.scores_ == co2.scores_ &&
           co1.offsets_ == co2.offsets_ &&
           co1.links_ == co2.links_; }

private:
  
  // Remove outliers to compute spread. See the WARNING in the cc!
  pair<int,int> CoreSpread( const vec<int> &starts ) const;
  
  // Methods to help logging info.
  String ToStringSupers( ) const;
  String ToStringCluster( const int cluster_id ) const;
  String ToStringND( const normal_distribution &nd ) const;
  
  
private:
  
  int super1_;  // first super
  int super2_;  // second super (signed to capture orientation)
  int slen1_;   // length of super1_
  int slen2_;   // length of super2_
  
  mutable vec<float> scores_;                 // scores of combined nds
  mutable vec<normal_distribution> offsets_;  // combined nd of each cluster
  mutable vec< vec<SLink> > links_;           // clusters of consistent links
  
};

#endif

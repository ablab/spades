///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef C_SUBLINK_H
#define C_SUBLINK_H

#include "PairsManager.h"
#include "SupersHandler.h"
#include "paths/Alignlet.h"

/**
 * class CSublink
 *
 * A class to encapsulate a simplified version of link between
 * scaffolds, one of which is strictly contained in the other (based
 * only on the linking information).
 */
class CSublink {
  
public:

  CSublink( const int bid = -1,
	    const int sid = -1,
	    const bool rc = false,
	    const int start = 0,
	    const int dev = 0,
	    const int weight = 0 );
  
  // Return false if neither scaffold is contained in the other.
  bool Set( const int &pair_sep,
	    const int &pair_stdev,
	    const shandler &supers,
	    const alignlet &hit1,
	    const alignlet &hit2,
	    const int weight = 1 );
  
  void SetStart( const int start ) { start_ = start; }

  int BigId( ) const { return big_id_; }

  int SmallId( ) const { return small_id_; }

  Bool SmallRc( ) const { return small_rc_; }

  int Start( ) const { return start_; }

  int Dev( ) const { return dev_; }

  int Weight( ) const { return weight_; }
  
  // Check if this link is consistent with the given other link.
  bool IsConsistentWith( const CSublink &other,
			 const double max_stretch ) const;
  
  // Implied stretch by using this start instead of start_.
  double Stretch( const int start ) const;

  // Returns -1 on error, or else the id of the closest gap (on big).
  int ClosestGap( const vec<superb> &scaffolds ) const;
  
  // Print some info for this link.
  void PrintInfo( const vec<superb> &scaffolds, ostream &out ) const;

  // operator<
  friend bool operator< ( const CSublink &left, const CSublink &right );
  
  
private:
  
  int big_id_;    // id of large scaffold (oriented fw, it contains small)
  int small_id_;  // id of small scaffold (contained in big)
  bool small_rc_; // orientation of small implied by link
  int start_;     // start of small on big implied by link
  int dev_;       // deviation implied by link
  int weight_;    // weight for this link

};

#endif

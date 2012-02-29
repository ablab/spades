///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "PairsManager.h"
#include "SupersHandler.h"
#include "lookup/LookAlign.h"
#include "paths/Alignlet.h"
#include "paths/CSublink.h"

/**
 * CSublink
 * Constructor
 */
CSublink::CSublink( const int bid,
		    const int sid,
		    const bool rc,
		    const int start,
		    const int dev,
		    const int weight ) :
  big_id_ ( bid ),
  small_id_ ( sid ),
  small_rc_ ( rc ),
  start_ ( start ),
  dev_ ( dev ),
  weight_ ( weight )
{ }

/**
 * CSublink
 * Set
 */
bool CSublink::Set( const int &pair_sep,
		    const int &pair_stdev,
		    const shandler &supers,
		    const alignlet &hit1,
		    const alignlet &hit2,
		    const int weight )
{
  // Detect which scaffold is big and which is small. Return false on failure!
  const alignlet *big_hit = 0;
  const alignlet *small_hit = 0;
  int sup_id1 = supers.ToSuper( hit1.TargetId( ) );
  int sup_id2 = supers.ToSuper( hit2.TargetId( ) );
  if ( sup_id1 < 0 || sup_id2 < 0 ) return false;

  int len1 = supers.TrueLength( sup_id1 );
  int len2 = supers.TrueLength( sup_id2 );
  if ( len1 < len2 ) {
    small_hit = &hit1;
    big_hit = &hit2;
  }
  else if ( len2 < len1 ) {
    small_hit = &hit2;
    big_hit = &hit1;
  }
  else
    return false;
  
  // Set ids and orientation.
  big_id_ = supers.ToSuper( big_hit->TargetId( ) );
  small_id_ = supers.ToSuper( small_hit->TargetId( ) );
  small_rc_ = ( big_hit->Fw1( ) == small_hit->Fw1( ) );
  
  // Set weight.
  weight_ = weight;
  
  // Implied start of small with respect to big.
  int big_contig = big_hit->TargetId( );
  int big_super = supers.ToSuper( big_contig );
  int big_superlen = supers.TrueLength( big_super );
  int big_begin = big_hit->Pos2( ) + supers.StartOnSuper( big_contig );
  int big_len = big_hit->Pos2( ) - big_hit->pos2( );
  int big_end = big_begin + big_len;
  int small_contig = small_hit->TargetId( );
  int small_super = supers.ToSuper( small_contig );
  int small_superlen = supers.TrueLength( small_super );
  int small_begin = small_hit->Pos2( ) + supers.StartOnSuper( small_contig );
  int small_len = small_hit->Pos2( ) - small_hit->pos2( );
  int small_end = small_begin + small_len;
  
  dev_ = pair_stdev;
  if ( big_hit->Fw1( ) ) {
    start_ = big_end + pair_sep;
    start_ -= small_hit->Fw1( ) ? small_superlen - small_end : small_begin;
  }
  else {
    start_ = big_begin - pair_sep;
    start_ -= small_hit->Fw1( ) ? small_end : small_superlen  - small_begin;
  }
  
  // Is small truly embedded?
  return ( start_ >= 0 & start_ + small_superlen <= big_superlen );
  
}

/**
 * CSublink
 * IsConsistentWith
 */
bool CSublink::IsConsistentWith( const CSublink &other,
				 const double max_stretch ) const
{
  if ( this->SmallId( ) != other.SmallId( ) ) return false;
  if ( this->BigId( ) != other.BigId( ) ) return false;
  if ( this->SmallRc( ) != other.SmallRc( ) ) return false;
  
  double dev = Min( this->Dev( ), other.Dev( ) );
  double delta = Abs( this->Start( ) - other.Start( ) );
  return ( delta <= dev * max_stretch );
}

/**
 * CSublink
 * Stretch
 */
double CSublink::Stretch( const int start ) const
{
  return ( double( start - start_ ) / double( dev_ ) );
}

/**
 * CSublink
 * ClosestGap
 *
 * Returns -1 on error, or else the id of the closest gap (on big).
 */
int CSublink::ClosestGap( const vec<superb> &scaffolds ) const
{
  const superb &sup = scaffolds[big_id_];
  if ( sup.Ntigs( ) < 2 ) return -1;
  
  int onsuper_pos = sup.Len( 0 );
  int distance = Abs( start_ - onsuper_pos );
  for (int ii=1; ii<sup.Ntigs( )-1; ii++) {
    onsuper_pos += sup.Gap( ii-1 ) + sup.Len( ii );
    if ( Abs( start_ - onsuper_pos ) > distance ) return ( ii - 1 );
    else distance = Abs( start_ - onsuper_pos );
  }
  
  return sup.Ntigs( ) - 2;
}

/**
 * CSublink
 * PrintInfo
 */
void CSublink::PrintInfo( const vec<superb> &scaffolds, ostream &out ) const
{
  out << "s" << small_id_
      << " (" << scaffolds[small_id_].TrueLength( )
      << " bp, " << scaffolds[small_id_].Ntigs( )
      << " edges) -> s" << big_id_
      << " (" << scaffolds[big_id_].TrueLength( )
      << " bp, " << scaffolds[big_id_].Ntigs( )
      << " edges), " << ( small_rc_ ? "rc" : "fw" )
      << " on [" << start_
      << ", " << start_ + scaffolds[small_id_].TrueLength( )
      << ") +/- " << dev_
      << ", weight=" << weight_
      << "\n";
}

/**
 * CSublink
 * operator<
 */
bool operator< ( const CSublink &left, const CSublink &right )
{
  if ( left.BigId( ) < right.BigId( ) ) return true;
  if ( left.BigId( ) > right.BigId( ) ) return false;
  if ( left.SmallId( ) < right.SmallId( ) ) return true;
  if ( left.SmallId( ) > right.SmallId( ) ) return false;
  if ( left.Start( ) < right.Start( ) ) return true;
  if ( left.Start( ) > right.Start( ) ) return false;
  if ( left.Dev( ) < right.Dev( ) ) return true;
  if ( left.Dev( ) > right.Dev( ) ) return false;
  return ( left.SmallRc( ) < right.SmallRc( ) );
}

///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "SupersHandler.h"
#include "paths/SInsertion.h"

/**
 * SInsertion
 * Constructor
 */
SInsertion::SInsertion( ) :
  small_id_ ( -1 ),
  big_id_ ( -1 ),
  small_rc_ ( false ),
  pos_in_big_ ( -1 ),
  gap_before_ ( make_pair( 0, 0 ) ),
  gap_after_ ( make_pair( 0, 0 ) )
{ }

/**
 * SInsertion
 * Constructor
 */
SInsertion::SInsertion( int small_id,
			int big_id,
			bool small_rc,
			int pos_in_big,
			pair<int,int> gap_before,
			pair<int,int> gap_after) :
  small_id_ ( small_id ),
  big_id_ ( big_id ),
  small_rc_ ( small_rc ),
  pos_in_big_ ( pos_in_big ),
  gap_before_ ( gap_before ),
  gap_after_ ( gap_after )
{ }

/**
 * SInsertion
 * Constructor
 */
SInsertion::SInsertion( const vec<superb> &scaffolds,
			const CSublink &link,
			const int gap_id,
			const int small_start,
			const int gap_len )
{
  small_id_ = link.SmallId( );
  big_id_ = link.BigId( );
  small_rc_ = link.SmallRc( );
  pos_in_big_ = gap_id;
  
  const superb &big_sup = scaffolds[big_id_];
  const superb &small_sup = scaffolds[small_id_];
  
  double linkdev = link.Dev( ) * link.Dev( );
  double scafdev = big_sup.Dev( gap_id ) * big_sup.Dev( gap_id );
  int common_dev = sqrt( linkdev + scafdev );
  
  int gap_end = 0;
  for (int ii=0; ii<=gap_id; ii++)
    gap_end += big_sup.Len( ii ) + big_sup.Gap( ii );
  int gap_begin = gap_end - big_sup.Gap( gap_id );
  int gap_size_before = small_start - gap_begin;
  gap_before_ = make_pair( gap_size_before, common_dev );
  
  int gap_size_after = gap_len - small_sup.TrueLength( ) - gap_size_before;
  gap_after_ = make_pair( gap_size_after, common_dev );
}

/**
 * SInsertion
 * PrintInfo
 */
void SInsertion::PrintInfo( const vec<superb> &scaffolds, ostream &out ) const
{
  const superb& ssup = scaffolds[small_id_];
  const superb& bsup = scaffolds[big_id_];
  
  int gap_end = 0;
  for (int ii=0; ii<=pos_in_big_; ii++)
    gap_end += bsup.Len( ii ) + bsup.Gap( ii );
  int gap_begin = gap_end - bsup.Gap( pos_in_big_ );
  
  out << "inserting s" << small_id_
      << " in s" << big_id_
      << " (" << bsup.TrueLength( )
      << " bp, " << bsup.Ntigs( )
      << " tigs) at gap " << pos_in_big_
      << " (" << gap_end - gap_begin
      << " +/- " << bsup.Dev( pos_in_big_ )
      << ") [" << gap_begin
      << ", " << gap_end
      << "): <" << gap_before_.first
      << " +/- " << gap_before_.second
      << "> <s" << small_id_ << ( small_rc_ ? "-" : "+" )
      << " (" << ssup.TrueLength( )
      << " bp, " << ssup.Ntigs( )
      << " tigs)> <" << gap_after_.first
      << " +/- " << gap_after_.second
      << ">";
  
  int gap_len_1 = gap_end - gap_begin;
  int gap_len_2 = gap_before_.first + ssup.TrueLength( ) + gap_after_.first;
  ForceAssertGe( gap_len_2, gap_len_1 );
  if ( gap_len_2 > gap_len_1 )
    out << " (gap was enlarged of " << gap_len_2 - gap_len_1
	<< " bases)";
  
  out << "\n";
}

bool operator< ( const SInsertion &left, const SInsertion &right )
{
  if ( left.small_id_ < right.small_id_ ) return true;
  if ( left.small_id_ > right.small_id_ ) return false;
  if ( left.big_id_ < right.big_id_ ) return true;
  if ( left.big_id_ > right.big_id_ ) return false;
  return ( left.pos_in_big_ < right.pos_in_big_ );
}


///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Vec.h"
#include "PrettyPrintTable.h"
#include "SupersHandler.h"
#include "lookup/LookAlign.h"
#include "paths/reporting/CGapStats.h"
#include "paths/reporting/EvalScaffoldsUtils.h"

/**
 * struct SPlacement
 * Constructors
 */
SPlacement::SPlacement( ) :
  target_ ( -1 ), query_ ( -1 ), rc_ ( false )
{
  twin_ = make_pair( 0, 0 );
  spos_ = make_pair( 0, 0 );
}

SPlacement::SPlacement( int query, int pos1, int pos2 ) :
  target_ ( -1 ), query_ ( query ), rc_ ( false )
{
  twin_ = make_pair( 0, 0 );
  spos_ = make_pair( pos1, pos2 );
}

SPlacement::SPlacement( const look_align_plus &hit1,
			const look_align_plus &hit2,
			int super_id, int pos1, int pos2 )
{
  this->SetFromLookAligns( hit1, hit2, super_id, pos1, pos2 );
}

SPlacement::SPlacement( int target, int query, bool rc,
			pair<int,int> twin, pair<int,int> spos ) :
  target_ ( target ), query_ ( query ), rc_ ( rc ),
  twin_ ( twin ), spos_ ( spos )
{ }
  
/**
 * struct SPlacement
 * SetFromLookAligns
 */
void SPlacement::SetFromLookAligns( const look_align_plus &hit1,
				    const look_align_plus &hit2,
				    const int super_id,
				    const int pos1,
				    const int pos2 )
{
  target_ = hit1.target_id;
  query_ = super_id;
  rc_ = hit1.Rc1( );
  
  // twin_.
  if ( hit1.Fw1( ) ) twin_ = make_pair( hit1.pos2( ), hit2.Pos2( ) );
  else twin_ = make_pair( hit2.pos2( ), hit1.Pos2( ) );
  
  // spos_.
  spos_ = make_pair( pos1, pos2 );
}

/**
 * struct SPlacement
 * MergesWith
 */
bool SPlacement::MergesWith( const double discrepancy, 
			     const SPlacement &other,
			     const shandler &supers )
{
  ForceAssert( query_ == query_ );
  if ( target_ != other.target_ ) return false;
  if ( rc_ != other.rc_ ) return false;
  
  // Can only merge adjacent blocks.
  if ( this->spos_.second != other.spos_.first ) return false;
  
  // Length of super from one cgpos to the other.
  double super_len = 0;
  const superb &sup = supers[query_];
  for (int ii=spos_.first; ii<other.spos_.second; ii++) {
    super_len += double( sup.Len( ii ) );
    if ( ii < other.spos_.second - 1 )
      super_len += double( sup.Gap( ii ) );
  }
  
  // Excessive compression.
  double correction = discrepancy * super_len;
  double min_slen = super_len - correction;
  double max_slen = super_len + correction;
  double observed_len
    = rc_
    ? double( twin_.second - other.twin_.first )
    : double( other.twin_.second - twin_.first );
  
  if ( observed_len < min_slen || max_slen < observed_len ) return false;
  
  // Ok! Merge.
  if ( rc_ ) twin_.first = other.twin_.first;
  else twin_.second = other.twin_.second;
  spos_.second = other.spos_.second;
  
  return true;
}
  
/**
 * struct SPlacement
 * PrintLegend
 */
void SPlacement::PrintLegend( ostream &out ) const
{
  out << "LEGEND:\n"
      << " w (weight): number of contigs in this stretch\n"
      << " sl (super length): gapped length of portion of super\n"
      << " rl (reference length): length of aligned part on reference\n"
      << " c (compression): ratio between sl and rl\n"
      << "\n";
}

/**
 * struct SPlacement
 * PrintInfo
 */
void SPlacement::PrintInfo( const shandler &supers, vec<String> &info ) const
{
  info.clear( );
  const bool failed = ( target_ < 0 );
  
  const superb &sup = supers[query_];
  longlong super_len = 0;
  for (int ii=spos_.first; ii<spos_.second; ii++) {
    super_len += double( sup.Len( ii ) );
    if ( ii < spos_.second - 1 )
      super_len += double( sup.Gap( ii ) );
  }
  longlong onref_len = twin_.second - twin_.first;
  double compression = double( onref_len ) / double( super_len );

  String str_entry
    = "s" + ToString( query_ )
    + "." + ToString( spos_.first ) + "-" + ToString(  spos_.second - 1 )
    + "/" + ToString( sup.Ntigs( ) - 1 );
  info.push_back( str_entry );
    
  str_entry
    = ( rc_ ? "- on t" : "+ on t" ) + ToString( target_ )
    + " [" + ToString( twin_.first ) +  "," + ToString( twin_.second ) + ")";
  if ( failed ) info.push_back( "na" );
  else info.push_back( str_entry );

  str_entry = "w " + ToString( spos_.second - spos_.first );
  info.push_back( str_entry );
    
  str_entry = "sl " + ToString( super_len );
  info.push_back( str_entry );
	     
  str_entry = "rl " + ToString( onref_len );
  if ( failed ) info.push_back( "na" );
  else info.push_back( str_entry );
      
  str_entry = "c " + ToString( compression, 1 );
  if ( failed ) info.push_back( "na" );
  else info.push_back( str_entry );

  bool tag = false;
  if ( failed ) tag = true;
  else if ( spos_.first > 0 ) tag = true;
  else if ( spos_.second < sup.Ntigs( ) ) tag = true;
  if ( tag ) str_entry += "\t#####";
  info.push_back( str_entry );
    
}

/**
 * struct SPlacement
 * operator<
 */
bool operator< ( const SPlacement &left, const SPlacement &right )
{
  if ( left.target_ < right.target_ ) return true;
  if ( left.target_ > right.target_ ) return false;
  if ( left.twin_ < right.twin_ ) return true;
  if ( left.twin_ > right.twin_ ) return false;
  if ( left.query_ < right.query_ ) return true;
  if ( left.query_ > right.query_ ) return false;
  if ( left.spos_ < right.spos_ ) return true;
  if ( left.spos_ > right.spos_ ) return false;
  return ( left.rc_ < right.rc_ );
}

/**
 * AreConsistent
 */
bool AreConsistent( const look_align_plus &hit1,
		    const look_align_plus &hit2,
		    const shandler &supers,
		    String &info,
		    bool weak,
		    CGapStats *gap_stats )
{
  // HEURISTICS
  const double max_stretch = 5.0;   // tag gaps with excessive stretch

  // On different supers.
  int cg1 = hit1.query_id;
  int cg2 = hit2.query_id;
  int s1 = supers.ToSuper( cg1 );
  int s2 = supers.ToSuper( cg2 );
  if ( s1 < 0 || s1 != s2 ) {
    info = "on_different_supers";
    return false;
  }
  
  // Wrong orientation.
  if ( hit1.Fw1( ) != hit2.Fw1( ) ) {
    info = "orient_on_super";
    return false;
  }
  
  // This could be improved: hits could be consistent even if non-adjacent.
  int super_pos1 = supers.PosOnSuper( cg1 );
  int super_pos2 = supers.PosOnSuper( cg2 );
  if ( Abs( super_pos2 - super_pos1 ) != 1 ) {
    info = "non_adjacent_on_super";
    return false;
  }

  // Relative mapping is wrong.
  bool is_ok = true;
  int beg1 = hit1.pos2( );
  int beg2 = hit2.pos2( );
  if ( hit1.Fw1( ) ) {
    if ( beg1 < beg2 && super_pos2 < super_pos1 ) is_ok = false;
    else if ( beg1 > beg2 && super_pos2 > super_pos1 ) is_ok = false;
  }
  else {
    if ( beg1 < beg2 && super_pos1 < super_pos2 ) is_ok = false;
    else if ( beg1 > beg2 && super_pos1 > super_pos2 ) is_ok = false;
  }
  if ( ! is_ok ) {
    info = "relative_placement";
    return false;
  }
  
  // Gap size and deviation.
  int end1 = hit1.Pos2( );
  int end2 = hit2.Pos2( );
  int observed = ( beg1 < beg2 ) ? beg2 - end1 : beg1 - end2;
  int gap_pos = ( super_pos1 < super_pos2 ) ? super_pos1 : super_pos2;
  int gap_size = supers[s1].Gap( gap_pos );
  int gap_dev = supers[s1].Dev( gap_pos );
  double stretch = ( double( observed - gap_size ) / double( gap_dev ) ) ;
  bool tag_bad = Abs( stretch ) > max_stretch;

  // Fill gap_stats.
  if ( gap_stats ) gap_stats->Add( observed, gap_size, gap_dev );
  
  // Enough for the weak test.
  if ( weak ) return true;

  // Fill info and return.
  info
    = "<" + ToString( gap_size ) + " +/- "  + ToString( gap_dev ) + ">"
    + "   real=" + ToString( observed )
    + "   stretch=" + ToString( stretch, 1 )
    + ( tag_bad ? "   #####" : "" );

  return ( ! tag_bad );

}

/**
 * SuperChain
 */
void SuperChain( const int super_id,
		 const shandler &supers,
		 const vec<look_align_plus> &hits,
		 const vec< vec<int> > &edge2hits,
		 vec<int> &chain )
{
  const superb &sup = supers[super_id];

  // First run to build a chain (a list of hits for the given super).
  chain.clear( );
  chain.resize( sup.Ntigs( ), -1 );
  for (int cgpos=0; cgpos<sup.Ntigs( ); cgpos++) {
    int contig_id = sup.Tig( cgpos );
    const vec<int> &ids = edge2hits[contig_id];
    if ( ids.size( ) == 1 ) chain[cgpos] = ids[0];
    else if ( ids.size( ) > 1 ) chain[cgpos] = -2;
  }
  
  // Incrementally select more chain events for multiply placed contigs.
  String info;
  while ( 1 ) {
    int n_events = 0;
    for (int cgpos=0; cgpos<sup.Ntigs( ); cgpos++) {
      if ( chain[cgpos] != -2 ) continue;
      int contig_id = sup.Tig( cgpos );
      const vec<int> &ids = edge2hits[contig_id];
      
      // Try twice to select a compatible hit from multiply placed contigs.
      for (int pass=0; pass<2; pass++) {
	bool weak = ( pass == 0 ? false : true );

	// Try to anchor to the previous contig.
	if ( cgpos > 1 && chain[cgpos-1] > -1 ) {
	  for (int jj=0; jj<ids.isize( ); jj++) {
	    const look_align_plus &hit1 = hits[ chain[cgpos-1] ];
	    const look_align_plus &hit2 = hits[ ids[jj] ];
	    if ( AreConsistent( hit1, hit2, supers, info, weak ) ) {
	      chain[cgpos] = ids[jj];
	      break;
	    }
	  }
	}
	if ( chain[cgpos] > -1 ) break;
	
	// Try to anchor to the next contig.
	if ( cgpos < sup.Ntigs( ) - 1  && chain[cgpos+1] > -1 ) {
	  for (int jj=0; jj<ids.isize( ); jj++) {
	    const look_align_plus &hit1 = hits[ ids[jj] ];
	    const look_align_plus &hit2 = hits[ chain[cgpos+1] ];
	    if ( AreConsistent( hit1, hit2, supers, info, weak ) ) {
	      chain[cgpos] = ids[jj];
	      break;
	    }
	  }
	}
	if ( chain[cgpos] > -1 ) break;
	
      }
      
      // Have we found a valid chain?
      if ( chain[cgpos] > -1 ) n_events++;
    }

    // If no chains were added, we can leave the while loop.
    if ( n_events < 1 ) break;
    
  }

}

/**
 * DigestChain
 */
void DigestChain( const shandler &supers,
		  const vec<look_align_plus> &hits,
		  const vec<int> &chain,
		  vec<SPlacement> &placs )
{
  placs.clear( );

  // HEURISTIC - Merge placements if this does not cause excessive stretch.
  const double discrepancy = 0.25;

  // Early exit for empty chains.
  int sel_id = -1;
  for (int ii=0; ii<chain.isize( ); ii++) {
    if ( chain[ii] > -1 ) {
      sel_id = ii;
      break;
    }
  }
  if ( sel_id < 0 ) return;

  // Super id.
  const int super_id = supers.ToSuper( hits[ chain[sel_id] ].query_id );
  const superb &sup = supers[super_id];

  // Initial setup: all contigs are accounted for.
  int post = 0;
  String info;
  while ( post < chain.isize( ) ) {
    
    // Interval of unaligned contigs.
    if ( chain[post] < 0 ) {
      int end = post + 1;
      while( end < chain.isize( ) && chain[end] < 0 ) end++;
      placs.push_back( SPlacement( super_id, post, end ) );
      post = end;
      continue;
    }
    
    // Interval of consistent contigs.
    int end = post + 1;
    while ( end < chain.isize( ) ) {
      if ( chain[end] < 0 ) break;
      const look_align_plus &hit1 = hits[ chain[end-1] ];
      const look_align_plus &hit2 = hits[ chain[end] ];
      if ( ! AreConsistent( hit1, hit2, supers, info ) ) break;
      end++;
    }
    const look_align_plus &hit1 = hits[ chain[post] ];
    const look_align_plus &hit2 = hits[ chain[end-1] ];
    placs.push_back( SPlacement( hit1, hit2, super_id, post, end ) );
    post = end;
    
  }

  // Merge consecutive placements, if possible (heurstics here).
  for (int ii=1; ii<placs.isize( ); ii++) {
    if ( placs[ii-1].MergesWith( discrepancy, placs[ii], supers ) ) {
      placs.erase( placs.begin( ) + ii );
      ii--;
    }
  }
  
}

/**
 * PrintChain
 */
void PrintChain( const shandler &supers,
		 const vec<look_align_plus> &hits,
		 const vec<int> &chain,
		 int &n_local_inversions,
		 ostream &out )
{
  n_local_inversions = 0;

  // Early exit for unaligned supers.
  int sel_id = -1;
  for (int ii=0; ii<chain.isize( ); ii++) {
    if ( chain[ii] > -1 ) {
      sel_id = ii;
      break;
    }
  }
  if ( sel_id < 0 ) return;

  // Super id.
  const int super_id = supers.ToSuper( hits[ chain[sel_id] ].query_id );
  const superb &sup = supers[super_id];

  // Generate printable table.
  vec<String> line;
  vec< vec<String> > table;
  const String str_sid = ToString( super_id );
  const String str_ntigs = ToString( sup.Ntigs( ) );
  const String str_pounds = "#####";
  String str_info;
  
  // Loop over all contigs in super.
  for (int cgpos=0; cgpos<sup.Ntigs( ); cgpos++) {
    const int contig_id = sup.Tig( cgpos );
    const int hit_id = chain[cgpos];
    
    // Contig id.
    line.clear( );
    line.push_back( "c" + ToString( contig_id ) );

    // Position of contig in super.
    line.push_back( "s" + str_sid + "." + ToString( cgpos ) + "/" + str_ntigs );
      
    // Contig length.
    line.push_back( ToString( sup.Len( cgpos ) ) + " bp" );

    // Mapping on reference.
    String str_map;
    if ( hit_id < 0 )
      str_map = "no align found";
    else {
      const look_align_plus &hit = hits[hit_id];
      str_map
	= ( hit.Fw1( ) ? "+ on t" : "- on t" ) + ToString( hit.target_id )
	+ " [" + ToString(  hit.pos2() ) + "," + ToString(  hit.Pos2() ) + ")";
    }
    line.push_back( str_map );
    
    // Consistency with previous hit.
    String str_gap = "";
    if ( hit_id > -1 && cgpos < chain.isize( )-1 ) {
      const look_align_plus &this_hit = hits[hit_id];
      const int next_hit_id = chain[cgpos+1];
      if ( next_hit_id < 0 ) 
	str_gap = str_pounds + " unaligned";
      else {
	const look_align_plus &next_hit = hits[next_hit_id];
	bool ok = AreConsistent( this_hit, next_hit, supers, str_info );
	if ( str_info.Contains( "relative_placement" ) ) n_local_inversions++;
	str_gap += str_info;
      }
    }
    line.push_back( str_gap );
    
    // Add line to table.
    table.push_back( line );
  }

  // Print table.
  const String separator = "   ";
  vec<Bool> rjust( 5, 0 );
  rjust[0] = true;
  rjust[1] = true;
  rjust[2] = true;
  
  String str_true_len = ToString( supers.TrueLength( super_id ) );
  String str_reduced_len = ToString( supers[super_id].ReducedLength( ) );
  out << "super_" << ToString( super_id )
      << separator << "gapped length = " << str_true_len << " bp"
      << separator << "ungapped length = " << str_reduced_len << " bp"
      << "\n";
  BeautifyAndPrintTable( table, out, separator, &rjust );
  out << endl;

}

/**
 * ChainGaps
 */
void ChainGaps( const shandler &supers,
		const vec<look_align_plus> &hits,
		const vec<int> &chain,
		CGapStats &gap_stats )
{
  // Early exit for unaligned supers.
  int sel_id = -1;
  for (int ii=0; ii<chain.isize( ); ii++) {
    if ( chain[ii] > -1 ) {
      sel_id = ii;
      break;
    }
  }
  if ( sel_id < 0 ) return;

  // Super id.
  const int super_id = supers.ToSuper( hits[ chain[sel_id] ].query_id );
  const superb &sup = supers[super_id];
  
  // Loop over all contigs in super.
  for (int cgpos=0; cgpos<sup.Ntigs( ); cgpos++) {
    const int hit_id = chain[cgpos];
    if ( hit_id > -1 && cgpos > 0 ) {
      const look_align_plus &this_hit = hits[hit_id];
      const int prev_hit_id = chain[cgpos-1];
      if ( prev_hit_id < 0 ) continue;
      
      const look_align_plus &prev_hit = hits[prev_hit_id];
      String str_bogus;
      AreConsistent( prev_hit, this_hit, supers, str_bogus, true, &gap_stats);
    }
  }
  
}


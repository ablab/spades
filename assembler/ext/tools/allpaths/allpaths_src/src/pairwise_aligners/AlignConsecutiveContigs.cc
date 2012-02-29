/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "Alignment.h"
#include "Basevector.h"
#include "PrintAlignment.h"
#include "ScoreAlignment.h"
#include "Superb.h"
#include "pairwise_aligners/AlignTwoBasevectors.h"
#include "pairwise_aligners/AlignConsecutiveContigs.h"

/**
 * AlignConsecutiveContigs
 */
bool AlignConsecutiveContigs( const int super_id,
			      const int cgpos,
			      vecbvec &contigs,
			      vec<superb> &supers,
			      ostream *log,
			      float MAX_ERROR_RATE )
{
  // Consts (could be turned into extra args).
  int max_stretch = 6;   // used to chop relevant portion of contigs for align
  int min_cglen = 2000;  // min size of chopped portion
  
  // Log stream.
  ofstream devnull( "/dev/null" );
  ostream &out = ( log ) ? *log : devnull;
  
  // Core data.
  const superb &sup = supers[super_id];
  const bvec &contig1 = contigs[sup.Tig( cgpos )];
  const bvec &contig2 = contigs[sup.Tig( cgpos + 1 )];
  const int gap_size = sup.Gap( cgpos );
  const int gap_stdev = sup.Dev( cgpos );
  const int len1 = contig1.size( );
  const int len2 = contig2.size( );
  
  // Find window of overlap on contig1.
  int wbeg1 = 0;
  int wend1 = 0;
  int wsize1 = 0;
  if ( len1 < min_cglen ) {
    wbeg1 = 0;
    wend1 = len1;
    wsize1 = len1;
  }
  else {
    wbeg1 = Max( 0, len1 + gap_size - max_stretch * gap_stdev );
    wend1 = len1;
    wsize1 = wend1 - wbeg1;
    if ( wsize1 < min_cglen ) {
      wbeg1 = Max( 0, wend1 - min_cglen );
      wsize1 = wend1 - wbeg1;
    }
  }

  // Find window of overlap on contig2.
  int wbeg2 = 0;
  int wend2 = 0;
  int wsize2 = 0;
  if ( len2 < min_cglen ) {
    wbeg2 = 0;
    wend2 = len2;
    wsize2 = len2;
  }
  else {
    wbeg2 = 0;
    wend2 = Min( len2, - gap_size + max_stretch * gap_stdev );
    wsize2 = wend2 - wbeg2;
    if ( wsize2 < min_cglen ) {
      wend2 = Min( len2, min_cglen );
      wsize2 = wend2 - wbeg2;
    }
  }

  // Gap size/stdev do not seem to imply an overlap.
  if ( wsize1 < 24 || wsize2 < 24 )
    return false;

  // Build basevectors and qualvectors from the original contigs.
  basevector mtig1;
  basevector mtig2;
  qualvector qual1( 0 );
  qualvector qual2( 0 );
  
  mtig1.SetToSubOf( contig1, wbeg1, wsize1 );
  mtig2.SetToSubOf( contig2, wbeg2, wsize2 );
  
  // Perform alignment.
  out << " aligning [" << sup.Tig( cgpos )
      << ", " << sup.Tig( cgpos + 1 )
      << "] gap=(" << gap_size
      << ", " << gap_stdev
      << ") aligned_lens=(" << mtig1.size( )
      << ", " << mtig2.size( )
      << "): " << flush;
  
  align al;
  int Rc = 0;
  int base_bandwidth = 100;
  int t1_bandwidth = mtig1.size( ) / 2;
  int t2_bandwidth = mtig2.size( ) / 2;
  int bandwidth = Min( base_bandwidth, t1_bandwidth, t2_bandwidth );
  int mlen1 = mtig1.size( );
  int mlen2 = mtig2.size( );
  
  int overlap = AlignTwoBasevectors( mtig1, mtig2, al, 0, 10000000, 1.0,
				     0, Rc, 8, 24, 2, 1, qual1, qual2,
				     10000000.0, 10000, devnull, False, 1000,
				     10000, 1000, 50, True, bandwidth );

  // An overlap value of less than 10 is not trusted.  In that case,
  //  or if no overlap was found at all, try again with different
  //  parameters. This second stage tends to find short overlaps
  //  missed by the first stage.
  if ( overlap < 10 ) {
    out << " align not found. Trying again:";
    overlap = AlignTwoBasevectors( mtig1, mtig2, al, 0, 10000000, 1.0,
				   0, Rc, 8, 24, 2, 1, qual1, qual2,
				   10000000.0, 10000, devnull, False, 1000,
				   10000, 1000, 50 );
  }
  
  // Adjust alignment.
  al.Setpos1( al.pos1( ) + wbeg1 ); 

  // Remove improper aligns.
  if ( overlap > 0 ) {
    if ( al.pos1( ) > 0 && al.pos2( ) > 0 ) overlap = -1;
    if ( al.Pos1( ) < len1 && al.Pos2( ) < len2 ) overlap = -1;
    float err = al.Errors( mtig1, mtig2 );
    float al_len = al.Pos1( ) - al.pos1( );
    if ( ( err / al_len ) > MAX_ERROR_RATE ) overlap = -1;
  }
  
  // No alignment found.
  if ( overlap < 0 ) {
    out << " no align found.\n";
    return false;
  }
  
  // Print alignment and return.
  out << " found align! len=" << al.Pos1( ) - al.pos1( )
      << ", n_errors=" << ToString( al.Errors( contig1, contig2 ) )
      << "\n"
      << " on1: " << al.pos1( ) << "-" << al.Pos1( ) << "_" << len1 << "\n"
      << " on2: " << al.pos2( ) << "-" << al.Pos2( ) << "_" << len2 << "\n"
      << endl;
  
  basevector merged;
  qualvector qmerged;
  qualvector q1( contig1.size( ), 50 );
  qualvector q2( contig2.size( ), 50 );
  PrintVisualAlignment( True, out, contig1, contig2, al, q1, q2 );
  out << flush;

  // Merge contigs and update supers.
  MergeTwoBaseVectors( contig1, contig2, al, merged, q1, q2, qmerged );
  contigs[sup.Tig( cgpos )] = merged;
  contigs[sup.Tig( cgpos + 1 )].resize( 0 );

  if ( cgpos == sup.Ntigs( ) - 2 ) {
    supers[super_id].RemoveTigByPos( cgpos + 1 );
    supers[super_id].SetLen( cgpos, (int)contig1.size( ) );
  }
  else {
    int orig_next_gap = sup.Gap( cgpos + 1 );
    int orig_next_dev = sup.Dev( cgpos + 1 );
    supers[super_id].RemoveTigByPos( cgpos + 1 );
    supers[super_id].SetLen( cgpos, (int)contig1.size( ) );
    
    // Check if contigs (cgpos+1) is embedded in cgpos.
    if ( sup.Gap( cgpos ) + sup.Len( cgpos + 1 ) > 0 ) {
      supers[super_id].SetGap( cgpos, orig_next_gap );
      supers[super_id].SetDev( cgpos, orig_next_dev );
    }
  }

  return true;
  
}

/**
 * order_supers_length
 * 
 * order super by TrueLength (longer supers come first).
 */
struct order_supers_length
  : public binary_function<const superb&, const superb&, bool>
{
public:
  bool operator( ) ( const superb &left, const superb &right ) {
    return ( left.TrueLength( ) > right.TrueLength( ) );
  }
};

/**
 * CompactifyContigs
 */
void CompactifyContigs( vecbvec &contigs,
			vec<superb> &supers )
{
  // Sort supers by length.
  order_supers_length sorter;
  sort( supers.begin( ), supers.end( ), sorter );
  
  // Adjust contigs and supers.
  vecbvec new_contigs;
  new_contigs.reserve( contigs.size( ) );
  
  vec<int> to_old;
  to_old.reserve( contigs.size( ) );
  for (int ii=0; ii<supers.isize( ); ii++) {
    for (int cpos=0; cpos<supers[ii].Ntigs( ); cpos++) {
      int new_id = to_old.size( );
      to_old.push_back( supers[ii].Tig( cpos ) );
      supers[ii].SetTig( cpos, new_id );
    }
  }

  for (int ii=0; ii<to_old.isize( ); ii++) {
    ForceAssert( contigs[ to_old[ii] ].size( ) > 0 );
    new_contigs.push_back( contigs[ to_old[ii] ] );
  }

  swap( contigs, new_contigs );

  // Check.
  for (int ii=0; ii<supers.isize( ); ii++) {
    for (int cpos=0; cpos<supers[ii].Ntigs( ); cpos++) {
      int slen = supers[ii].Len( cpos );
      int clen = contigs[ supers[ii].Tig( cpos ) ].size( );
      ForceAssertEq( slen, clen );
    }
  }
  
}

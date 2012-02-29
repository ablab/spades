///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "CoreTools.h"
#include "MergeReads2.h"
#include "PackAlign.h"
#include "PairsManager.h"
#include "PrintAlignment.h"
#include "Qualvector.h"
#include "VecUtilities.h"
#include "math/HoInterval.h"
#include "kmers/KmerRecord.h"
#include "math/Functions.h"
#include "pairwise_aligners/SmithWatBandedA.h"
#include "paths/AssemblyEdit.h"
#include "paths/LongReadTools.h"

void CorrectPatch( const basevector& LEFT, const basevector& RIGHT,
     const vecbasevector& fbases, const vecqualvector& fquals, 
     const PairsManager& fpairs, const vec< kmer<20> >& fheads, 
     const vec<int64_t>& fids, assembly_edit& e, ostringstream& out,
     const Bool verbose )
{
     // Heuristics.

     const int F = 20;
     const int flank = 200;
     const int bandwidth = 10;
     const int eflank = 10;
     const int min_edits = 5;
     const int min_edits_plus = 5;
     const int min_mult = 5;
     const int min_mult2 = 3;

     int npasses = 1;
     for ( int pass = 1; pass <= npasses; pass++ )
     {

     // Set up.

     int start1 = e.Start1( ), stop2 = e.Stop2( );
     ForceAssertEq( e.Nreps( ), 1 );
     basevector& patch = e.Rep(0); 
     basevector left, right;
     left.SetToSubOf( LEFT, 0, start1 );
     right.SetToSubOf( RIGHT, stop2, (int) RIGHT.size( ) - stop2 );

     // Define extended patch.

     int left_ext = Min( flank, left.isize( ) );
     int right_ext = Min( flank, right.isize( ) );
     basevector epatch = Cat( basevector( left, left.isize( ) - left_ext, left_ext ),
          patch, basevector( right, 0, right_ext ) );
     if (verbose)
     {    out << "\n";
          epatch.Print( out, "epatch" );    }

     // Find hits.

     vec< pair<int,int64_t> > fw_hits, rc_hits;
     vec<ho_interval> cov;
     kmer<F> x;
     for ( int pass = 1; pass <= 2; pass++ )
     {    for ( int j = 0; j <= epatch.isize( ) - F; j++ )
          {    x.SetToSubOf( epatch, j );
               if ( pass == 2 ) x.ReverseComplement( );
               int64_t low = LowerBound( fheads, x ), high = UpperBound( fheads, x );
               for ( int64_t i = low; i < high; i++ )
                    ( pass == 1 ? fw_hits : rc_hits ).push( j, fids[i] );    }    }

     // Find alignments.

     vec< pair<int,int> > start_stop;
     vec< pair<int64_t,int64_t> > id1_id2;
     vec< pair<align,align> > aligns;
     for ( int i1 = 0; i1 < fw_hits.isize( ); i1++ )
     for ( int i2 = 0; i2 < rc_hits.isize( ); i2++ )
     {    int64_t id1 = fw_hits[i1].second, id2 = rc_hits[i2].second;
          if ( fpairs.getPartnerID(id1) != id2 ) continue;
          int start = fw_hits[i1].first, stop = rc_hits[i2].first + F;
          if ( stop < start ) continue;
          start_stop.push( start, stop );
          id1_id2.push( id1, id2 );
          align a1, a2;
          int errors;
          SmithWatBandedA( 
               fbases[id1], epatch, -start, bandwidth, a1, errors, 0, 1, 1 );
          CenterMobileGaps( a1, fbases[id1], epatch, false, out );
          basevector b = fbases[id2];
          b.ReverseComplement( );
          SmithWatBandedA( 
               b, epatch, -( stop - b.isize( ) ), bandwidth, a2, errors, 0, 1, 1 );
          CenterMobileGaps( a2, b, epatch, false, out );
          aligns.push( a1, a2 );    }

     // Find edits.

     vec< pair<ho_interval,basevector> > edits;
     vec<int> xids;
     for ( int z = 0; z < aligns.isize( ); z++ )
     {    int start = start_stop[z].first, stop = start_stop[z].second;
          int64_t id1 = id1_id2[z].first, id2 = id1_id2[z].second;
          const align &a1 = aligns[z].first, &a2 = aligns[z].second;
          if (verbose)
          {    out << "\n";
               PRINT4_TO( out, start, stop, id1, id2 );
               PrintVisualAlignment(
                    True, out, fbases[id1], epatch, a1, fquals[id1] );    }
          vec<ho_interval> perf1, perf2;
          a1.PerfectIntervals1( fbases[id1], epatch, perf1 );
          a1.PerfectIntervals2( fbases[id1], epatch, perf2 );
          for ( int j = 0; j < perf2.isize( ); j++ )
               if ( perf2[j].Length( ) >= 2 * eflank ) cov.push_back( perf2[j] );
          for ( int j = 0; j < perf1.isize( ) - 1; j++ )
          {    if ( perf1[j].Length( ) < eflank ) continue;
               int k;
               for ( k = j + 1; k < perf1.isize( ); k++ )
                    if ( perf1[k].Length( ) >= eflank ) break;
               if ( k == perf1.isize( ) ) continue;
               int estart = perf2[j].Stop( ), estop = perf2[k].Start( );
               basevector s( fbases[id1], perf1[j].Stop( ), 
                    perf1[k].Start( ) - perf1[j].Stop( ) );
               edits.push( ho_interval(estart,estop), s );
               xids.push_back(z);
               if (verbose)
               {    out << "edit " << estart << "-" << estop << " " << s.ToString( )
                         << "\n";    }
               j = k - 1;    }
          basevector b = fbases[id2];
          b.ReverseComplement( );
          qualvector q = fquals[id2];
          q.ReverseMe( );
          if (verbose) PrintVisualAlignment( True, out, b, epatch, a2, q );    
          a2.PerfectIntervals1( b, epatch, perf1 );
          a2.PerfectIntervals2( b, epatch, perf2 );
          for ( int j = 0; j < perf2.isize( ); j++ )
               if ( perf2[j].Length( ) >= 2 * eflank ) cov.push_back( perf2[j] );
          for ( int j = 0; j < perf1.isize( ) - 1; j++ )
          {    if ( perf1[j].Length( ) < eflank ) continue;
               int k;
               for ( k = j + 1; k < perf1.isize( ); k++ )
                    if ( perf1[k].Length( ) >= eflank ) break;
               if ( k == perf1.isize( ) ) continue;
               int estart = perf2[j].Stop( ), estop = perf2[k].Start( );
               basevector s( b, perf1[j].Stop( ), 
                    perf1[k].Start( ) - perf1[j].Stop( ) );
               edits.push( ho_interval(estart,estop), s );
               xids.push_back(z);
               if (verbose)
               {    out << "edit " << estart << "-" << estop << " " << s.ToString( )
                         << "\n";    }
               j = k - 1;    }    }

     // Identify and apply poison.

     SortSync( edits, xids );
     vec<int> poison;
     for ( int i = 0; i < edits.isize( ); i++ )
     {    int j = edits.NextDiff(i);
          int ecount = j - i, refcount = 0;
          if ( ecount >= min_edits )
          {    for ( int l = 0; l < cov.isize( ); l++ )
               {    if ( cov[l].Start( ) <= edits[i].first.Start( ) - eflank
                         && cov[l].Stop( ) >= edits[i].first.Stop( ) + eflank )
                    {    refcount++;    }    }
               if ( ecount >= min_edits_plus && refcount >= min_edits_plus 
                    && ecount < min_mult * refcount )
               {    if (verbose)
                    {    out << "poison: " << ecount << " EDIT " 
                              << edits[i].first.Start( )
                              << "-" << edits[i].first.Stop( ) << " " 
                              << edits[i].second.ToString( ) 
                              << " refcount=" << refcount << "\n";    }
                    for ( int k = i; k < j; k++ )
                         poison.push_back( xids[k] );    }    }
          i = j - 1;    }
     UniqueSort(poison);
     vec<Bool> to_delete( edits.size( ), False );
     for ( int i = 0; i < edits.isize( ); i++ )
          if ( BinMember( poison, xids[i] ) ) to_delete[i] = True;
     EraseIf( edits, to_delete );

     // Tally edits and decide which ones pass.

     if (verbose) out << "\n";
     vec< pair<ho_interval,basevector> > edits_pass;
     vec<int> edits_pass_count;
     for ( int i = 0; i < edits.isize( ); i++ )
     {    int j = edits.NextDiff(i);
          int ecount = j - i, refcount = 0;
          if ( ecount >= min_edits )
          {    for ( int l = 0; l < cov.isize( ); l++ )
               {    if ( cov[l].Start( ) <= edits[i].first.Start( ) - eflank
                         && cov[l].Stop( ) >= edits[i].first.Stop( ) + eflank )
                    {    refcount++;    }    }
               if ( ecount >= min_mult * refcount )
               {    if (verbose)
                    {    out << ecount << " EDIT " << edits[i].first.Start( ) << "-" 
                              << edits[i].first.Stop( ) << " " 
                              << edits[i].second.ToString( ) 
                              << " refcount=" << refcount << "\n";    }
                    edits_pass.push_back( edits[i] );    
                    edits_pass_count.push_back( j - i );    }    }
          i = j - 1;    }
     if (verbose) out << "\n";

     // Identify nonoverlapping edits.

     vec< pair<ho_interval,basevector> > edits_accept;
     for ( int i1 = 0; i1 < edits_pass.isize( ); i1++ )
     {    Bool overlap = False;
          for ( int i2 = 0; i2 < edits_pass.isize( ); i2++ )
          {    if ( i2 == i1 ) continue;
               if ( Distance( edits_pass[i1].first, edits_pass[i2].first ) == 0
                    && edits_pass_count[i1] < min_mult2 * edits_pass_count[i2] )
               {    
                    // If the first edit is a subset of the second one, we're OK.

                    int start1 = edits_pass[i1].first.Start( );
                    int stop1 = edits_pass[i1].first.Stop( );
                    int start2 = edits_pass[i2].first.Start( );
                    int stop2 = edits_pass[i2].first.Stop( );
                    int start = Min( start1, start2 ), stop = Max( stop1, stop2 );
                    basevector s1 = edits_pass[i1].second;
                    basevector s2 = edits_pass[i2].second;
                    basevector b, s;
                    if ( start1 < start2 )
                    {    b.SetToSubOf( epatch, start1, start2 - start1 );
                         s2 = Cat( b, s2 );    }
                    if ( start2 < start1 )
                    {    b.SetToSubOf( epatch, start2, start1 - start2 );
                         s1 = Cat( b, s1 );    }
                    if ( stop1 < stop2 )
                    {    b.SetToSubOf( epatch, stop1, stop2 - stop1 );
                         s1 = Cat( s1, b );    }
                    if ( stop2 < stop1 )
                    {    b.SetToSubOf( epatch, stop2, stop1 - stop2 );
                         s2 = Cat( s2, b );    }
                    s.SetToSubOf( epatch, start, stop - start );
                    align a;
                    int e1, e2, e12;
                    if ( s.size( ) > 0 && s1.size( ) > 0 )
                         e1 = SmithWatFreeSym( s, s1, a, True, True, 1, 1 );
                    else e1 = Max( s.size( ), s1.size( ) );
                    if ( s.size( ) > 0 && s2.size( ) > 0 )
                         e2 = SmithWatFreeSym( s, s2, a, True, True, 1, 1 );
                    else e2 = Max( s.size( ), s2.size( ) );
                    if ( s1.size( ) > 0 && s2.size( ) > 0 )
                         e12 = SmithWatFreeSym( s1, s2, a, True, True, 1, 1 );
                    else e12 = Max( s1.size( ), s2.size( ) );
                    if ( e2 != e1 + e12 ) overlap = True;    }    }
          if ( !overlap ) edits_accept.push_back( edits_pass[i1] );    }

     // Make edits.
     
     int npatch = patch.size( );
     for ( int i = edits_accept.isize( ) - 1; i >= 0; i-- ) 
     {    ho_interval h = edits_accept[i].first - left_ext;
          const basevector& s = edits_accept[i].second;
          if ( h.Start( ) < 0 || h.Stop( ) > npatch ) continue;
          basevector b1( patch, 0, h.Start( ) );
          basevector b2( patch, h.Stop( ), patch.isize( ) - h.Stop( ) );
          patch = Cat( b1, s, b2 );    }    }    }

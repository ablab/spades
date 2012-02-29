///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Close gaps using Mixmers.  Totally experimental.

// int UID = 265; // 1 closures, covering 71500-71662
// int UID = 364; // 2 closures, covering 71294-71481
// int UID = 255; // 6 closures, covering 71175-71526

// int UID = 291; // 2 closures, covering nothing ?? (no, rc problem)
// int UID = 430; // 13 closures, covering nothing

// int UID = 223; // 0 closures
// int UID = 175; // 4 closures, covering 90782-90902, 91329-91449
// int UID = 469; // 3 closures, covering 91066-91181

// gap [3] 68659-68652
// int start1 = 68400, start2 = 68700;
// not closed, in fact underlying unibases don't cover
// oh, this is probably because of error in reference....

// gap [4] 70243-70245
// int start1 = 71050, start2 = 71300; 
// 6 closures, covers

// gap [5] 71421-71432
// int start1 = 71245, start2 = 71386; 
// 3 closures, covers

// gap [6] 90679-90652
// int start1 = 90450, start2 = 90750;
// 27 closures, covering:
// 90485-90735
// 90730-90905

// gap [9] 91197-91199
// int start1 = 91000, start2 = 91250;
// 24 closures, covering:
// 91092-91216
// 91156-91431

// gap [11] 91840-91851
// int start1 = 91600, start2 = 91900;
// 16 closures, covering:
// 91762-92041

bool allow_filled = false;

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "Basevector.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "Qualvector.h"
#include "graph/Digraph.h"
#include "graph/DigraphTemplate.h"
#include "math/HoInterval.h"
#include "paths/GetNexts.h"
#include "paths/HyperBasevector.h"
#include "paths/PairedPair.h"
#include "paths/ReadsToPathsCoreX.h"
#include "paths/Unipath.h"
#include "util/ReadTracker.h"
#include "util/SearchFastb2Core.h"

// CorrectErrors.

void CorrectErrors( vecbasevector& bases, vecqualvector& quals,
     const int flank, const int qmult )
{    int K = 2*flank + 1;
     for ( int left_flank = 0; left_flank < K; left_flank++ )
     {    // if ( left_flank != flank ) continue; // this neuters the for loop!
          vec<basevector> kmer;
          vec<int> id, pos, q;
          for ( size_t i = 0; i < bases.size( ); i++ )
          {    for ( int j = 0; j <= bases[i].isize( ) - K; j++ )
               {    basevector b( bases[i], j, K );
                    kmer.push_back(b), id.push_back(i), pos.push_back(j);
                    q.push_back( quals[i][j+left_flank] );    }    }
          SortSync( kmer, id, pos, q );
          vec<basevector> kmer_fixed(kmer);
          for ( int i = 0; i < kmer.isize( ); i++ )
          {    int j = kmer.NextDiff(i);
               vec<int> qsum( 4, 0 ), aid( 4, vec<int>::IDENTITY );
               for ( int r = i; r < j; r++ )
                    qsum[0] += q[r];
               for ( int add = 1; add <= 3; add++ )
               {    basevector b( kmer[i] );
                    b.Set( left_flank, ( b[left_flank] + add ) % 4 );
                    int s1 = lower_bound( 
                         kmer.begin( ), kmer.end( ), b ) - kmer.begin( );
                    int s2 = upper_bound( 
                         kmer.begin( ), kmer.end( ), b ) - kmer.begin( );
                    for ( int r = s1; r < s2; r++ )
                         qsum[add] += q[r];    }
               ReverseSortSync( qsum, aid );
               if ( aid[0] != 0 && qsum[0] > qmult * qsum[1] && qsum[1] < 30 )
               {    for ( int r = i; r < j; r++ )
                    {    kmer_fixed[r].Set( left_flank, 
                              ( kmer[i][left_flank] + aid[0] ) % 4 );    }    }
               i = j - 1;    }
          for ( int i = 0; i < kmer.isize( ); i++ )
          {    if ( kmer_fixed[i] != kmer[i] )
               {    bases[ id[i] ].Set( pos[i] + left_flank, 
                         kmer_fixed[i][left_flank] );
                    quals[ id[i] ][ pos[i] + left_flank ] = 0;    }    }    }    }

// IsPeriodic: determine if a basevector is periodic with periodicity between
// 1 and top_period.  For example,
// AAAAAA is periodic with period 1
// ATATATA is periodic with period 2.

Bool IsPeriodic( const basevector& x, const int top_period )
{    for ( int p = 1; p <= top_period; p++ )
     {    Bool fail = False;
          for ( int j = 0; j < p; j++ )
          {    if (fail) break;
               for ( int k = j + p; k < x.isize( ); k += p )
               {    if ( x[k] != x[j] ) 
                    {    fail = True;
                         break;    }    }    }
          if ( !fail ) return True;    }
     return False;    }

Bool IsPeriodicEx( const basevector& x, const int top_period, const int max_per )
{    int mpp = max_per + 1;
     for ( int j = 0; j <= x.isize(  ) - mpp; j++ )
     {    basevector b( x, j, mpp );
          if ( IsPeriodic( b, top_period ) ) return True;    }
     return False;    }

// MakeAperiodicMers: given a basevector K, scan it from left to right, breaking it
// up into a series of sub-basevectors, each of length at least K, and for which
// neither the left nor the right (K-1)-mer is periodic with period <= top_period.
// This code appends.

void MakeAperiodicMers( const basevector& x, const int K, const int top_period,
     vec<basevector>& mers )
{    basevector left, right, b;
     for ( int i = 0; i <= x.isize( ) - K; i++ )
     {    left.SetToSubOf( x, i, K-1 );
          if ( IsPeriodicEx( left, top_period, K/2 ) ) continue;
          for ( int len = K; len <= x.isize( ) - i; len++ )
          {    right.SetToSubOf( x, i + len - (K-1), K-1 );
               if ( IsPeriodicEx( right, top_period, K/2 ) ) continue;
               b.SetToSubOf( x, i, len );
               mers.push_back(b);
               break;    }    }    }

// BuildUnibasesFromSeqs: Build "unibases" from a set of sequences.
// Strictly speaking the outputs are not unibases.

void BuildUnibasesFromSeqs( const int K, vec<basevector> seqs,
     vecbasevector& unibases )
{
     // Build graph data.

     UniqueSort(seqs);
     vec< vec<int> > nexts, backs;
     if ( K == 10 ) GetNexts<10>( seqs, nexts );
     else if ( K == 20 ) GetNexts<20>( seqs, nexts );
     else
     {    cout << "That value of K is unsupported." << endl;
          exit(1);    }
     vec<basevector> seqs_rc(seqs);
     for ( size_t i = 0; i < seqs_rc.size( ); i++ )
          seqs_rc[i].ReverseComplement( );
     if ( K == 10 ) GetNexts<10>( seqs_rc, backs );
     else if ( K == 20 ) GetNexts<20>( seqs_rc, backs );

     // Build unibases.

     vec<Bool> used( seqs.size( ), False );
     unibases.clear( );
     for ( size_t i = 0; i < seqs.size( ); i++ )
     {    if ( used[i] ) continue;
          vec<int> chain;
          chain.push_back(i);

          // Go backwards.

          size_t j = i;
          Bool cycle = False;
          while(1)
          {    if ( backs[j].size( ) != 1 ) break;
               else
               {    int l = backs[j][0];
                    if ( nexts[l].size( ) != 1 ) break;
                    j = l;
                    if ( j == i )
                    {    cycle = True;
                         break;    } 
                    chain.push_back(j);    }    }
          chain.ReverseMe( );

          // Go forwards.

          if ( !cycle )
          {    j = chain.back( );
               while(1)
               {    if ( nexts[j].size( ) != 1 ) break;
                    else
                    {    int l = nexts[j][0];
                         if ( backs[l].size( ) != 1 ) break;
                         j = l;
                         chain.push_back(j);    }    }    }

          // Build answer.

          basevector b = seqs[ chain[0] ];
          for ( int r = 1; r < chain.isize( ); r++ )
          {    b.resize( b.isize( ) - (K-1) );
               b = Cat( b, seqs[ chain[r] ] );    }
          for ( int r = 0; r < chain.isize( ); r++ )
               used[ chain[r] ] = True;
          unibases.push_back(b);    }    }

void CleanEnds( HyperBasevector& h, const vec<basevector>& valid_kmers_1, 
     const vec<basevector>& valid_kmers_2, const Bool verbose )
{    int K = h.K( );
     for ( int count = 1; count <= 4; count++ )
     {    vec<int> first_good( h.EdgeObjectCount( ), 1000000000 );
          vec<int> last_good( h.EdgeObjectCount( ), -1 );
          vec<int> edges_to_delete;
          for ( int z = 0; z < h.N( ); z++ )
          {    if ( h.Source(z) )
               {    for ( int r = 0; r < h.From(z).isize( ); r++ )
                    {    basevector& e = h.EdgeObjectByIndexFromMutable( z, r );
                         int ei = h.EdgeObjectIndexByIndexFrom( z, r );
                         Bool has_valid_start = False;
                         for ( int s = 0; s <= e.isize( ) - K; s++ )
                         {    basevector c;
                              c.SetToSubOf( e, s, K );
                              if ( BinMember( valid_kmers_1, c ) )
                              {    has_valid_start = True;
                                   first_good[ei] = s;
                                   break;    }    }
                         if ( !has_valid_start ) edges_to_delete.push_back(ei);
                         else if ( first_good[ei] > 0 )
                         {    e.SetToSubOf( e, first_good[ei],
                                   e.isize( ) - first_good[ei] );    }    }    }
               if ( h.Sink(z) )
               {    for ( int r = 0; r < h.To(z).isize( ); r++ )
                    {    basevector& e = h.EdgeObjectByIndexToMutable( z, r );
                         int ei = h.EdgeObjectIndexByIndexTo( z, r );
                         Bool has_valid_stop = False;
                         // for ( int s = 0; s <= e.isize( ) - K; s++ )
                         for ( int s = e.isize( ) - K; s >= 0; s-- )
                         {    basevector c;
                              c.SetToSubOf( e, s, K );
                              if ( BinMember( valid_kmers_2, c ) )
                              {    has_valid_stop = True;
                                   last_good[ei] = s;
                                   break;    }    }
                         if ( !has_valid_stop )
                         {    int ei = h.EdgeObjectIndexByIndexTo( z, r );
                              edges_to_delete.push_back(ei);    }    
                         else if ( last_good[ei] < e.isize( ) - K )
                         {    e.resize( e.isize( ) - 
                                   ( e.isize( ) - K - last_good[ei] ) );    }
                                   }    }    }
          UniqueSort(edges_to_delete);
          if (verbose)
          {    cout << "CleanEnds deleting " << edges_to_delete.size( ) << " edges" 
                    << endl;    }
          h.DeleteEdges(edges_to_delete);
          h.RemoveUnneededVertices( );
          h.RemoveDeadEdgeObjects( );    }    }

Bool cmp2( const triple<int64_t,int64_t,int>& t1,
     const triple<int64_t,int64_t,int>& t2 )
{    if ( t1.second < t2.second ) return True;
     if ( t1.second > t2.second ) return False;
     int xpos1 = t1.third, xpos2 = t2.third;
     if ( xpos1 < 0 ) xpos1 = -xpos1-1;
     if ( xpos2 < 0 ) xpos2 = -xpos2-1;
     if ( xpos1 < xpos2 ) return True;
     if ( xpos1 > xpos2 ) return False;
     if ( t1.first < t2.first ) return True;
     return False;    }

Bool cmp3( const triple<int64_t,int64_t,int>& t1,
     const triple<int64_t,int64_t,int>& t2 )
{    int xpos1 = t1.third, xpos2 = t2.third;
     // if ( xpos1 < 0 ) xpos1 = -xpos1-1;
     // if ( xpos2 < 0 ) xpos2 = -xpos2-1;
     if ( xpos1 < xpos2 ) return True;
     if ( xpos1 > xpos2 ) return False;
     if ( t1.second < t2.second ) return True;
     if ( t1.second > t2.second ) return False;
     if ( t1.first < t2.first ) return True;
     return False;    }

void FindLocs( const vec< vec< pair<int,ho_interval> > >& truth,
     const pp_closure& c, vec<ho_interval>& locs )
{    locs.clear( );
     for ( int i = 0; i < truth.isize( ); i++ )
     {    for ( int j = 0; j <= truth[i].isize( ) - c.isize( ); j++ )
          {    Bool mismatch = False;
               for ( int k = 0; k < c.isize( ); k++ )
               {    if ( truth[i][j+k].first != c[k] )
                    {    mismatch = True;
                         break;    }    }
               if (mismatch) continue;
               locs.push( truth[i][j].second.Start( ), 
                    truth[i][ j + c.isize( ) - 1 ].second.Stop( ) );    }    }    }

void MakeReads( const vec< pair<int,int> > & IDS, const vecbasevector& bases,
     const vecqualvector& quals, vecbasevector& bases1, vecbasevector& bases2, 
     vecqualvector& quals1, vecqualvector& quals2, const Bool verbose )
{    for ( int i = 0; i < IDS.isize( ); i++ )
     {    int id1 = IDS[i].first, id2 = IDS[i].second;
          bases1.push_back_reserve( bases[id1] );
          basevector b( bases[id2] );
          b.ReverseComplement( );
          bases2.push_back_reserve(b);
          quals1.push_back_reserve( quals[id1] );
          qualvector q( quals[id2] );
          q.ReverseMe( );
          quals2.push_back_reserve(q);    }
     {    vecbasevector bases12(bases1);
          bases12.Append(bases2);
          vecqualvector quals12(quals1);
          quals12.Append(quals2);
          const int correct_errors_flank = 6;
          const int correct_errors_qmult = 10;
          CorrectErrors( bases12, quals12, 
               correct_errors_flank, correct_errors_qmult );
          for ( size_t i = 0; i < bases1.size( ); i++ )
          {    bases1[i] = bases12[i];
               bases2[i] = bases12[ i + bases1.size( ) ];
               quals1[i] = quals12[i];
               quals2[i] = quals12[ i + bases1.size( ) ];    }
          if (verbose)
          {    bases1.WriteAll( "xxx1.fastb" );
               bases2.WriteAll( "xxx2.fastb" );    }    }    }

void GetSoloKmers( const vecbasevector& bases1, const vecbasevector& bases2,
     const int K, vec<basevector>& solos )
{    vec<basevector> kmers;
     for ( size_t i = 0; i < bases1.size( ); i++ )
     {    for ( int j = 0; j <= bases1[i].isize( ) - K; j++ )
          {    basevector b( bases1[i], j, K );
               kmers.push_back(b);    }    }
     for ( size_t i = 0; i < bases2.size( ); i++ )
     {    for ( int j = 0; j <= bases2[i].isize( ) - K; j++ )
          {    basevector b( bases2[i], j, K );
               kmers.push_back(b);    }    }
     Sort(kmers);
     for ( int i = 0; i < kmers.isize( ); i++ )
     {    int j = kmers.NextDiff(i);
          if ( j - i == 1 ) solos.push_back( kmers[i] );
          i = j - 1;    }    }

// turned off for the moment, not necessarily the right thing to do
void GetDirectClosures( const vec<pp_pair>& ppairs, const vec<int>& L,
     const vec< vec< pair<int,ho_interval> > >& truth,
     vec<pp_closure>& all_closures )
{    for ( int i = 0; i < ppairs.isize( ); i++ )
     {    if ( ppairs[i].LeftSize( ) == 0 || ppairs[i].RightSize( ) == 0 ) continue;
          const double dmult = 3.0;
          cout << "\nPAIR " << i << " = " << ppairs[i] << endl;
          vec<pp_closure> closures;
          GetClosures( ppairs[i], L, dmult, closures );
          all_closures.append(closures);
          for ( int j = 0; j < closures.isize( ); j++ )
          {    cout << "closure " << j << " = " << closures[j];
               vec<ho_interval> locs;
               FindLocs( truth, closures[j], locs );
               if ( locs.nonempty( ) )
               {    cout << " (";
                    for ( int l = 0; l < locs.isize( ); l++ )
                    {    if ( l > 0 ) cout << ", ";
                         cout << locs[l];    }
                    cout << ")";    }
               cout << "\n";    }    }    }

void GetSingleReadPatches( const vec<pp_pair>& ppairs, const HyperBasevector& h,
     const vecbasevector& bases1, const vecbasevector& bases2,
     const vec<int>& L, vec<pp_closure>& all_closures, const Bool verbose )
{    int K = h.K( );
     for ( int i1 = 0; i1 < ppairs.isize( ); i1++ )
     {    if ( bases1[i1].isize( ) < K || bases2[i1].isize( ) < K ) continue;
          const pp_pair& p1 = ppairs[i1];
          if ( p1.LeftSize( ) == 0 || p1.RightSize( ) == 0 ) continue;
          vec<pp_closure> this_closures;
          for ( int i2 = 0; i2 < ppairs.isize( ); i2++ )
          {    if ( i2 == i1 ) continue;
               if ( bases1[i2].isize( ) < K || bases2[i2].isize( ) < K ) continue;
               const pp_pair& p2 = ppairs[i2];
               if ( p2.LeftSize( ) == 0 || p2.RightSize( ) == 0 ) continue;
               vec<pp_mpair> pnew;
               const double dmult = 3.0;
               FindJoins( p1, p2, L, dmult, pnew );
               for ( int j = 0; j < pnew.isize( ); j++ )
               {    const pp_pair& p = pnew[j];
                    if (verbose)
                    {    cout << "\njoin of pairs " << i1 << " and " << i2
                              << " yielding " << pnew[j] << endl;    }
                    vec<pp_closure> closures;
                    GetClosures( p, L, dmult, closures );
                    for ( int z = 0; z < closures.isize( ); z++ )
                    {    
                         basevector c = h.EdgeObject( closures[z][0] );
                         for ( int r = 1; r < closures[z].isize( ); r++ )
                         {    c.resize( c.isize( ) - (K-1) );
                              c = Cat( c, h.EdgeObject( closures[z][r] ) );    }
                         
                         if (verbose)
                         {    cout << "closure " << z << " = " 
                                   << closures[z] << endl;    }
                         vec<int> offsets_L1, offsets_R1, offsets_L2, offsets_R2;
                         GetOverlaps( p1.Left( ), closures[z], offsets_L1, True );
                         GetOverlaps( p1.Right( ), closures[z], offsets_R1, True );
                         GetOverlaps( p2.Left( ), closures[z], offsets_L2, True );
                         GetOverlaps( p2.Right( ), closures[z], offsets_R2, True );
                         for ( int L1 = 0; L1 < offsets_L1.isize( ); L1++ )
                         for ( int R1 = 0; R1 < offsets_R1.isize( ); R1++ )
                         for ( int L2 = 0; L2 < offsets_L2.isize( ); L2++ )
                         for ( int R2 = 0; R2 < offsets_R2.isize( ); R2++ )
                         {    if (verbose)
                              {    cout << "offsets = " << offsets_L1[L1] << ","
                                        << offsets_R1[R1] << ","
                                        << offsets_L2[L2] << ","
                                        << offsets_R2[R2] << endl;    }
                              int start = -offsets_L1[L1];
                              int stop = -offsets_R1[R1] + p1.RightSize( );
                              if ( start < stop )
                              {    pp_closure ce;
                                   for ( int u = start; u < stop; u++ )
                                        ce.push_back( closures[z][u] );
                                   this_closures.push_back(ce);    }
                              
                              // all_closures.push_back( closures[z] );

                              if (verbose)
                              {    int offL1 = 0, offL2 = 0, offR1 = 0, offR2 = 0;
                                   for ( int x = 0; x < -offsets_L1[L1]; x++ )
                                        offL1 += L[ closures[z][x] ];
                                   for ( int x = 0; x < -offsets_L2[L2]; x++ )
                                        offL2 += L[ closures[z][x] ];
                                   basevector b1(bases1[i1], 0, K); 
                                   basevector b2(bases1[i2], 0, K);
                                   int o1 = -offsets_L1[L1], o2 = -offsets_L2[L2];
                                   const basevector& 
                                        e1 = h.EdgeObject( closures[z][o1] );
                                   const basevector& 
                                        e2 = h.EdgeObject( closures[z][o2] );
                                   int start1L = e1.Find(b1), start2L = e2.Find(b2);
                                   offL1 += start1L, offL2 += start2L;
                                   int ox1 = -offsets_R1[R1] + p1.RightSize( );
                                   for ( int x = ox1; x < closures[z].isize( ); x++ )
                                        offR1 += L[ closures[z][x] ];
                                   int ox2 = -offsets_R2[R2] + p2.RightSize( );
                                   for ( int x = ox2; x < closures[z].isize( ); x++ )
                                        offR2 += L[ closures[z][x] ];
                                   basevector 
                                        d1( bases2[i1], bases2[i1].isize() - K, K );
                                   basevector 
                                        d2( bases2[i2], bases2[i2].isize() - K, K );
                                   const basevector& f1 
                                        = h.EdgeObject( closures[z][ox1-1] );
                                   const basevector& f2 
                                        = h.EdgeObject( closures[z][ox2-1] );
                                   int start1R = f1.Find(d1), start2R = f2.Find(d2);
                                   offR1 += ( f1.isize( ) - K - start1R );
                                   offR2 += ( f2.isize( ) - K - start2R );
                                   offR1 = c.isize( ) - bases2[i1].isize( ) - offR1;
                                   offR2 = c.isize( ) - bases2[i2].isize( ) - offR2;
                                   int m = Min( offL1, offR1, offL2, offR2 );
                                   if ( m < 0 ) m = -m;
                                   else m = 0;
                                   cout << String( m, ' ' ) << c.ToString( ) << endl;
                                   cout << String( m + offL1, ' ' ) 
                                        << bases1[i1].ToString( ) << endl;
                                   cout << String( m + offR1, ' ' ) 
                                        << bases2[i1].ToString( ) << endl;
                                   cout << String( m + offL2, ' ' ) 
                                        << bases1[i2].ToString( ) << endl;
                                   cout << String( m + offR2, ' ' ) 
                                        << bases2[i2].ToString( ) << endl;    }

                                   }    }    }    }

          UniqueSort(this_closures);
          if ( this_closures.solo( ) ) all_closures.append(this_closures);
          else if ( this_closures.size( ) > 1 && verbose ) 
          {    cout << "> 1 closure, rejecting\n";    }    }    }

void PrintUnibases( const HyperBasevector& h )
{    cout << "\nunibases:\n";
     for ( int j = 0; j < h.EdgeObjectCount( ); j++ )
          cout << BaseAlpha(j) << " = " << h.EdgeObject(j).ToString( ) << endl;
     cout << "\nconnections:\n";
     vec< pair<int,int> > edges;
     for ( int i = 0; i < h.N( ); i++ )
     {    for ( int j1 = 0; j1 < h.To(i).isize( ); j1++ )
          {    for ( int j2 = 0; j2 < h.From(i).isize( ); j2++ )
               {    int e1 = h.EdgeObjectIndexByIndexTo( i, j1 );
                    int e2 = h.EdgeObjectIndexByIndexFrom( i, j2 );
                    edges.push( e1, e2 );    }    }    }
     Sort(edges);
     for ( int i = 0; i < edges.isize( ); i++ )
     {    cout << BaseAlpha( edges[i].first ) << " --> "
               << BaseAlpha( edges[i].second ) << "\n";    }    }

void FindPlacements( const HyperBasevector& h, const String& data_dir,
     vec< vec< pair<int,ho_interval> > >& truth )
{    cout << "\nplacements:\n";
     {    vecbasevector edgesb;
          for ( int i = 0; i < h.EdgeObjectCount( ); i++ )
               edgesb.push_back_reserve( h.EdgeObject(i) );
          edgesb.WriteAll( "ccc.fastb" );    }
     vec< triple<int64_t,int64_t,int> > UALIGNS;
     SearchFastb2( "ccc.fastb", data_dir + "/genome.fastb", 20, &UALIGNS, 0,
          -1, 0.90, False );
     sort( UALIGNS.begin( ), UALIGNS.end( ), cmp3 );
     int last_stop = -1;
     vec< pair<int,ho_interval> > truthx;
     for ( int i = 0; i < UALIGNS.isize( ); i++ )
     {    int id = UALIGNS[i].first, start = UALIGNS[i].third;
          // if ( start < 0 ) start = -start-1;
          int stop = start + h.EdgeObject(id).isize( );
          if ( last_stop >= 0 && last_stop - start != h.K( ) - 1 )
          {    cout << "----------------------------\n";
               truth.push_back(truthx);
               truthx.clear( );    }
          last_stop = stop;
          cout << BaseAlpha(id) << " --> " << start << "-" << stop << " perfectly\n";
          truthx.push( id, ho_interval( start, stop ) );    }
     truth.push_back(truthx);
     flush(cout);    }

void ReportClosures( const vec<pp_closure>& all_closures,
     vec< vec< pair<int,ho_interval> > >& truth )
{    cout << "\n";
     int unmapped = 0;
     vec<ho_interval> all_locs;
     for ( int r = 0; r < all_closures.isize( ); r++ )
     {    cout << "closure = " << all_closures[r];
          vec<ho_interval> locs;
          FindLocs( truth, all_closures[r], locs );
          all_locs.append(locs);
          if ( locs.nonempty( ) )
          {    cout << " (";
               for ( int l = 0; l < locs.isize( ); l++ )
               {    if ( l > 0 ) cout << ", ";
                    cout << locs[l];    }
               cout << ")";    }
          else unmapped++;
          cout << "\n";    }
     UniqueSort(all_locs);
     cout << "\n// " << all_closures.size( ) 
          << " closures (" << unmapped << " unmapped), covering:" << endl;
     vec<Bool> locs_to_delete( all_locs.size( ), False );
     for ( int i1 = 0; i1 < all_locs.isize( ); i1++ )
     {    for ( int i2 = 0; i2 < all_locs.isize( ); i2++ )
          {    if ( all_locs[i1] != all_locs[i2] 
                    && Subset( all_locs[i1], all_locs[i2] ) )
               {    locs_to_delete[i1] = True;    }    }    }
     EraseIf( all_locs, locs_to_delete );
     int KK = 96;
     for ( int i = 0; i < all_locs.isize( ); i++ )
     {    int j;
          for ( j = i + 1; j < all_locs.isize( ); j++ )
               if ( all_locs[j-1].Stop( ) - all_locs[j].Start( ) < KK - 1 ) break;
          cout << "// " << all_locs[i].Start( ) << "-" 
               << all_locs[j-1].Stop( ) << "\n";
          i = j - 1;    }    }

void MapReadsToGraph( const vec< pair<int,int> >& IDS, const vecbasevector& bases1, 
     const vecbasevector& bases2, const vecqualvector& quals1, 
     const vecqualvector& quals2, const PairsManager& pairs, 
     const vec<int>& left_trim, const vec<int>& right_trim, HyperBasevector& h, 
     const vec<int>& L, vec<pp_pair>& ppairs, const Bool verbose )
{    int npairs = bases1.size( ), K = h.K( );
     vecbasevector bases2rc(bases2);
     ReverseComplement(bases2rc);
     vec<String> decomp1, decomp2;
     vec<int> to_left, to_right;
     h.ToLeft(to_left), h.ToRight(to_right);
     vec<pp_read> reads1, reads2;
     vec<int> left_extend( npairs, 0 ), right_extend( npairs, 0 );
     for ( int pass = 1; pass <= 2; pass++ )
     {    const vecbasevector& bases = ( pass == 1 ? bases1 : bases2rc );
          vec<basevector> kmrs;
          vec<int> vs, js, poss;
          for ( int v = 0; v < h.N( ); v++ )
          {    for ( int j = 0; j < h.From(v).isize( ); j++ )
               {    const basevector& e = h.EdgeObjectByIndexFrom( v, j );
                    for ( int pos = 0; pos <= e.isize( ) - K; pos++ )
                    {    basevector b( e, pos, K );
                         kmrs.push_back(b), vs.push_back(v), js.push_back(j),
                         poss.push_back(pos);    }    }    }
          SortSync( kmrs, vs, js, poss );
          for ( size_t i = 0; i < bases.size( ); i++ )
          {    vec<String> units;
               pp_read pp;
               basevector r = bases[i];
               basevector b( r, 0, K );
               int p1 = lower_bound( kmrs.begin( ), kmrs.end( ), b ) - kmrs.begin( );
               int p2 = upper_bound( kmrs.begin( ), kmrs.end( ), b ) - kmrs.begin( );
               if ( p2 - p1 != 1 ) 
                    units.push_back( "[" + ToString(p2-p1) + " starts]" );
               else 
               {    int v = vs[p1], j = js[p1], pos = poss[p1];
                    int e = h.EdgeObjectIndexByIndexFrom( v, j );
                    units.push_back( BaseAlpha(e) );
                    pp.push_back(e);
                    if ( pass == 1 ) left_extend[i] = pos;
                    else right_extend[i] = pos;
                    if ( pos + r.isize( ) > h.EdgeObject(e).isize( ) )
                    {    int start = h.EdgeObject(e).isize( ) - pos - (K-1);
                         if ( start < r.isize( ) )
                         {    r = basevector( r, start, r.isize( ) - start );
                              while( r.isize( ) >= K )
                              {    int w = h.From(v)[j];
                                   int nedges = h.From(w).size( );

                                   // For each edge leading out from w, count the 
                                   // number of agreeing bases before there is a 
                                   // disagreement.

                                   vec<int> agrees, ids(nedges, vec<int>::IDENTITY);
                                   for ( int m = 0; m < nedges; m++ )
                                   {    const basevector& e 
                                             = h.EdgeObjectByIndexFrom( w, m );
                                        int z;
                                        for ( z = 0; 
                                             z < Min( r.isize( ), e.isize( ) ); z++ )
                                        {    if ( r[z] != e[z] ) break;    }
                                        agrees.push_back(z);    }
                                   ReverseSortSync( agrees, ids );
                                   int nexts = 0;
                                   if ( nedges > 0 )
                                   {    for ( nexts = 1; nexts < nedges; nexts++ )
                                        {    if ( agrees[nexts] != agrees[0] ) 
                                                  break;    }    }

                                   if ( nexts != 1 )
                                   {    if (verbose)
                                        {    String u = "[" + ToString(nexts)
                                                  + " starts]";
                                             if ( pass == 1 ) u += r.ToString( );
                                             else
                                             {    basevector rrc(r);
                                                  rrc.ReverseComplement( );
                                                  u = rrc.ToString( ) + u;    }
                                             units.push_back(u);    }
                                        break;    }
                                   int e = h.EdgeObjectIndexByIndexFrom( w, ids[0] );
                                   v = w;
                                   j = ids[0];
                                   if (verbose) units.push_back( BaseAlpha(e) );
                                   pp.push_back(e);
                                   int pos = h.EdgeObject(e).isize( ) - (K-1);
                                   int len = r.isize( ) - pos;
                                   if ( len <= 0 ) r.resize(0);
                                   else r.SetToSubOf( r, pos, len );    
                                        }    }    }    }
               if ( pass == 2 ) 
               {    if (verbose) units.ReverseMe( );
                    pp.ReverseMe( );    }
               if ( pass == 1 ) reads1.push_back(pp);
               else reads2.push_back(pp);
               if (verbose)
               {    String d;
                    for ( int r = 0; r < units.isize( ); r++ )
                    {    if ( r > 0 ) d += ".";
                         d += units[r];    }
                    if ( pass == 1 ) decomp1.push_back(d);
                    else decomp2.push_back(d);    }    }
          h.Reverse( );    }
     for ( int i = 0; i < npairs; i++ )
     {    int pid = pairs.getPairID( IDS[i].first );
          int rl = bases1[i].size( ) + bases2[i].size( );
          int len = reads1[i].Length(L) + reads2[i].Length(L) + 2*(K-1);
          int sep = pairs.sep(pid) + rl - len - left_trim[i] - right_trim[i]
               + left_extend[i] + right_extend[i];
          ppairs.push( reads1[i], reads2[i], sep, pairs.sd(pid) );    }
     if (verbose)
     {    cout << "\nmapping reads back:\n";
          for ( size_t i = 0; i < bases1.size( ); i++ )
          {    ostringstream s;
               s << i;
               cout << "[" << i << "]: " << decomp1[i] << "\n";
               cout << " " << String( s.str( ).size( ), ' ' ) 
                    << "   " << decomp2[i] << "\n";    }    }    }

void BuildHyper( const vecbasevector& munibases, const int K,
     const vecbasevector& bases1, const vecbasevector& bases2,
     const int top_period, HyperBasevector& h, const vec<basevector>& solos,
     const Bool verbose )
{
     // Build graph.

     vec< vec<int> > from, to( munibases.size( ) );
     GetNexts( K, munibases, from );
     for ( size_t i = 0; i < munibases.size( ); i++ )
     {    Sort( from[i] );
          for ( int j = 0; j < from[i].isize( ); j++ )
               to[ from[i][j] ].push_back(i);    }
     digraph A( from, to );
     BuildUnibaseAdjacencyHyperBasevector( K, A, munibases, h );

     // Clean ends.

     vec<basevector> valid_kmers_1, valid_kmers_2;
     for ( size_t z = 0; z < bases1.size( ); z++ )
     {    basevector b1 = bases1[z], b2;
          if ( b1.isize( ) >= K )
          {    b1.resize(K);
               valid_kmers_1.push_back(b1);    }
          if ( bases2[z].isize( ) >= K )
          {    b2.SetToSubOf( bases2[z], bases2[z].isize( ) - K, K );
               valid_kmers_2.push_back(b2);    }    }
     UniqueSort(valid_kmers_1), UniqueSort(valid_kmers_2);
     CleanEnds( h, valid_kmers_1, valid_kmers_2, verbose );

     // Delete edges with periodic ends.  Not sure this makes sense.

     vec<int> to_delete;
     for ( int i = 0; i < h.EdgeObjectCount( ); i++ )
     {    const basevector& e = h.EdgeObject(i);
          basevector e1( e, 0, K-1 ), e2( e, e.isize( ) - (K-1), K-1 );
          if ( IsPeriodicEx( e1, top_period, K/2 ) 
               || IsPeriodicEx( e2, top_period, K/2 ) )
          {    to_delete.push_back(i);    }    }
     h.DeleteEdges(to_delete);
     h.RemoveUnneededVertices( );
     h.RemoveDeadEdgeObjects( );

     // Identify edges having kmers that appear only once in the reads.

     if (verbose) cout << "\n";
     to_delete.clear( );
     for ( int i = 0; i < h.EdgeObjectCount( ); i++ )
     {    int uniques = 0;
          for ( int j = 0; j <= h.EdgeObject(i).isize( ) - K; j++ )
          {    basevector b( h.EdgeObject(i), j, K );
               if ( BinMember( solos, b ) ) uniques++;    }
          if ( uniques > 0 )
          {    if (verbose)
               {    cout << uniques << " unique kmers in edge " 
                         << BaseAlpha(i) << "\n";    }
               to_delete.push_back(i);    }    }
     h.DeleteEdges(to_delete);
     h.RemoveUnneededVertices( );
     h.RemoveDeadEdgeObjects( );

     // Announce number of edges.

     if (verbose) cout << "\n#edges = " << h.EdgeObjectCount( ) << "\n";    }

void BuildClosuresFromReads( const vec< pair<int,int> >& IDS,
     const vecbasevector& bases, const vecqualvector& quals, 
     const PairsManager& pairs, const int K, const int top_period,
     const String& data_dir, vec<basevector>& CLOSURES, const Bool verbose )
{
     // Define reads.

     vecbasevector bases1, bases2;
     vecqualvector quals1, quals2;
     MakeReads( IDS, bases, quals, bases1, bases2, quals1, quals2, verbose );

     // Find kmers appearing only once in the reads.

     vec<basevector> solos;
     GetSoloKmers( bases1, bases2, K, solos );

     // Remove unique kmers from the beginnings of reads.

     int npairs = bases1.size( );
     vec<int> left_trim, right_trim;
     for ( int i = 0; i < npairs; i++ )
     {    int j;
          for ( j = 0; j <= bases1[i].isize( ) - K; j++ )
          {    basevector b( bases1[i], j, K );
               if ( !BinMember( solos, b ) ) break;    }
          bases1[i].SetToSubOf( bases1[i], j, bases1[i].isize( ) - j );
          left_trim.push_back(j);
          for ( j = bases2[i].isize( ) - K; j >= 0; j-- )
          {    basevector b( bases2[i], j, K );
               if ( !BinMember( solos, b ) ) break;    }
          bases2[i].SetToSubOf( bases2[i], 0, K + j );    
          right_trim.push_back( bases2[i].isize( ) - (K+j) );    }

     // Build mers.

     vec<basevector> mers;
     for ( int i = 0; i < npairs; i++ )
     {    MakeAperiodicMers( bases1[i], K, top_period, mers );
          MakeAperiodicMers( bases2[i], K, top_period, mers );    }
     /*
     cout << "mers:\n";
     for ( size_t i = 0; i < mers.size( ); i++ )
          cout << mers[i].ToString( ) << "\n"; 
     */

     // Build hyper.

     vecbasevector munibases;
     BuildUnibasesFromSeqs( K, mers, munibases );
     HyperBasevector h;
     BuildHyper( munibases, K, bases1, bases2, top_period, h, solos, verbose );

     // Map the reads back to the graph.

     vec<int> L;
     for ( int i = 0; i < h.EdgeObjectCount( ); i++ )
          L.push_back( h.EdgeObject(i).size( ) - K + 1 );
     vec<pp_pair> ppairs;
     MapReadsToGraph( IDS, bases1, bases2, quals1, quals2, pairs, left_trim, 
          right_trim, h, L, ppairs, verbose );

     // Print stuff.

     if (verbose)
     {    Ofstream( dout, "xxx.dot" );
          h.PrintSummaryDOT0w( dout, True, False, True, NULL, True );
          PrintUnibases(h);    }

     // Find placements of unibases on the genome.

     vec< vec< pair<int,ho_interval> > > truth;
     if (verbose) FindPlacements( h, data_dir, truth );

     // Try to close pairs.

     vec<pp_closure> all_closures;
     // GetDirectClosures( ppairs, L, truth, all_closures );
     GetSingleReadPatches( ppairs, h, bases1, bases2, L, all_closures, verbose );
     UniqueSort(all_closures);
     if (verbose) ReportClosures( all_closures, truth );
     for ( int i = 0; i < all_closures.isize( ); i++ )
     {    const pp_closure& p = all_closures[i];
          basevector b = h.EdgeObject( p[0] );
          for ( int j = 1; j < p.isize( ); j++ )
          {    b.resize( b.isize( ) - (K-1) );
               b = Cat( b, h.EdgeObject( p[j] ) );    }
          CLOSURES.push_back(b);    }    }

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String( PRE );
     CommandArgument_String( DATA );
     CommandArgument_String( RUN );
     CommandArgument_Int(K);
     CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
       "Number of threads to use (use all available processors if set to 0)");
     CommandArgument_String_OrDefault( UNIBASES_IN, "all_reads");
     CommandArgument_String_OrDefault( FRAG_READS_IN, "frag_reads_filt");
     CommandArgument_String_OrDefault( FRAG_READS_CORR_IN, "frag_reads_corr");
     CommandArgument_String_OrDefault( UNIBASES_OUT, "all_reads.added" );
     CommandArgument_Int_OrDefault( UID, -1 );
     CommandArgument_Int_OrDefault( START1, -1 );
     CommandArgument_Int_OrDefault( START2, -1 );
     CommandArgument_Bool_OrDefault( VERBOSE, False );
     EndCommandArguments;

     // Thread control
  
     NUM_THREADS = configNumThreads(NUM_THREADS);
     omp_set_num_threads( NUM_THREADS );

     String KS = ToString(K);

     // Define directories.

     String data_dir = PRE + "/" + DATA;
     String run_dir = data_dir + "/" + RUN;

     // Define heuristics.

     int evil_K = 20;   // was original just K, but that has a well established meaning in our code
     const int top_period = 5;

     // Load read pairs.

     PairsManager pairs( run_dir + "/frag_reads_filt.pairs" );
     PairsManager pairs_corr( run_dir + "/frag_reads_corr.pairs" );

     // Get number of reads.

     size_t nreads = MastervecFileObjectCount( 
          run_dir + "/frag_reads_filt.fastb" );
     size_t nreads_corr 
          = MastervecFileObjectCount( run_dir + "/frag_reads_corr.fastb" );

     // Generate map of "corr" reads to reads.

     vec<int> id_filt(nreads_corr);
     {    ReadTracker rt;
          rt.Load( run_dir + "/frag_reads_corr" );
          for ( size_t r = 0; r < nreads_corr; r++ )
               id_filt[r] = rt.GetReadIndex(r);    }
               
     // Figure out which pairs were filled.

     vec<Bool> filled( pairs.nPairs( ), False );
     {    String line;
          fast_ifstream fin( run_dir + "/filled_reads.info" );
          while(1)
          {    getline( fin, line );
               if ( fin.fail( ) ) break;
               int cpid = line.Before( " " ).Int( ); // .corr pair id
               int cid1 = pairs_corr.ID1(cpid);
               int id1 = id_filt[cid1];
               filled[ pairs.getPairID(id1) ] = True;    }    }

     // Align the unfilled pairs to the unibases.

     const int K0 = 26;
     const int MAX_PLACEMENTS = 100;
     vec< triple<int64_t,int64_t,int> > XALIGNS;
     {    String bases_tmp = run_dir + "/frag_reads_unfilled.tmp.fastb";
          {    vecbasevector bases( 
                    run_dir + "/frag_reads_filt.fastb" );
               for ( size_t id = 0; id < bases.size( ); id++ )
               {    int pid = pairs.getPairID(id);
                    if ( pid < 0 || filled[pid] ) bases[id].resize(0);
                    else bases[id].resize(K0);    }
               bases.WriteAll(bases_tmp);    }
          SearchFastb2( bases_tmp, run_dir +  "/" + UNIBASES_IN + ".unibases.k" + KS, K0,
               &XALIGNS, 0, MAX_PLACEMENTS );
          Remove(bases_tmp);
          Sort( XALIGNS, cmp2 );    }

     // Load the reads and unibases.

     vecbasevector bases( run_dir + "/frag_reads_filt.fastb" );
     vecqualvector quals( run_dir + "/frag_reads_filt.qualb" );
     vecbasevector unibases( run_dir + "/" + UNIBASES_IN + ".unibases.k" + KS );

     // Remove placed pairs and placements that are not close enough to the 
     // end of a unipath.  Also keep only forward facing alignments.

     size_t naligns = XALIGNS.size( );
     const int dev_mult = 4;
     vec<Bool> to_remove( XALIGNS.size( ), False );
     for ( size_t i = 0; i < XALIGNS.size( ); i++ )
     {    size_t j;
          for ( j = i + 1; j < XALIGNS.size( ); j++ )
               if ( XALIGNS[j].second != XALIGNS[i].second ) break;
          vec< triple<int,int,int> > pid_ind_pos;
          for ( size_t k = i; k < j; k++ )
          {    int rid = XALIGNS[k].first;
               pid_ind_pos.push( pairs.getPairID(rid), k, XALIGNS[k].third );    }
          Sort(pid_ind_pos);
          for ( int l = 0; l < pid_ind_pos.isize( ); l++ )
          {    int m;
               for ( m = l + 1; m < pid_ind_pos.isize( ); m++ )
                    if ( pid_ind_pos[m].first != pid_ind_pos[l].first ) break;
               if ( m - l == 2 )
               {    int pid = pid_ind_pos[l].first;
                    int id1 = pairs.ID1(pid), id2 = pairs.ID2(pid);
                    int sep = pairs.sep(pid), dev = pairs.sd(pid);
                    int n1 = bases[id1].size( ), n2 = bases[id2].size( );
                    int ind1 = pid_ind_pos[l].second, ind2 = pid_ind_pos[l+1].second;
                    int pos1 = pid_ind_pos[l].third, pos2 = pid_ind_pos[l+1].third;
                    if ( pos1 < 0 ) swap( pos1, pos2 );
                    if ( pos1 >= 0 && pos2 < 0 )
                    {    pos2 = -pos2-1 + K0;
                         int fraglen = pos2 - pos1;
                         int low = n1 + sep - dev_mult * dev + n2;
                         int high = n1 + sep + dev_mult * dev + n2;
                         if ( low <= fraglen && fraglen <= high )
                              to_remove[ind1] = to_remove[ind2] = True;    }    }
               l = m - 1;    }
          i = j - 1;    }
     for ( size_t i = 0; i < XALIGNS.size( ); i++ )
     {    int id1 = XALIGNS[i].first, u = XALIGNS[i].second, pos = XALIGNS[i].third;
          int id2 = pairs.getPartnerID(id1);
          int n1 = bases[id1].size( ), n2 = bases[id2].size( );
          int pid = pairs.getPairID(id1);
          int sep = pairs.sep(pid), dev = pairs.sd(pid);
          if ( pos >= 0 )
          {    if ( pos + n1 + sep + dev_mult * dev + n2 <= unibases[u].isize( ) )
                    to_remove[i] = True;    }
          else to_remove[i] = True;    }
     EraseIf( XALIGNS, to_remove );

     // Set up for closures.

     vec<basevector> CLOSURES;

     // Handle the cheating case where range of reads was provided.  Note that this
     // won't work correctly unless REMOVE_DODGY_READS_FRAG=False, because the 
     // numbering changes between 'orig' and 'filt' reads.

     if ( START1 >= 0 )
     {    vec< pair<int,int> > IDS;
          fast_ifstream in( run_dir + "/frag_reads_orig.info" );
          String line;
          // cout << "using reads:\n";
          while(1)
          {    getline( in, line );
               if ( in.fail( ) ) break;
               if ( line.Contains( "filled" ) && !allow_filled ) continue;
               if ( line.Contains( "no" ) ) continue;
               istrstream iline( line.c_str( ) );
               String read1, junk, loc, read2;
               iline >> read1 >> junk >> loc >> junk >> read2;
               int id1 = read1.Before( "fw" ).Int( );
               read2 = read2.GlobalReplaceBy( "rc", "" );
               int id2 = read2.Int( );
               int start = loc.Between( ".", "-" ).Int( );
               if ( start < START1 || start > START2 ) continue;
               if (VERBOSE) 
                    cout << "[" << IDS.size( ) << "] " << id1 << " " << id2 << "\n";
               IDS.push( id1, id2 );    }
          BuildClosuresFromReads( IDS, bases, quals, pairs, evil_K, 
               top_period, data_dir, CLOSURES, VERBOSE );    }

     // Handle the case where read clusters are derived from alignments.

     else
     {    // Report alignments.

          vec<size_t> XALIGNS_index( unibases.size( ) + 1 );
          {    size_t i = 0;
               for ( size_t u = 0; u <= unibases.size( ); u++ )
               {    while( i < XALIGNS.size( ) && XALIGNS[i].second < (int64_t) u ) 
                         ++i;
                    XALIGNS_index[u] = i;    }    }
          vec< vec<basevector> > CLOSURES_u( unibases.size( ) );
          #pragma omp parallel for
          for ( size_t u = 0; u < unibases.size( ); u++ )
          {    size_t i = XALIGNS_index[u], j = XALIGNS_index[u+1];
               if ( UID >= 0 && (int) u != UID ) { continue; }
               if (VERBOSE)
               {    cout << "\nunfilled reads pointing off end of unibase " 
                         << u << endl;    }
               for ( size_t k = i; k < j; k++ )
               {    int id1 = XALIGNS[k].first;
                    if (VERBOSE)
                    {    cout << id1 << " at " 
                              << XALIGNS[k].third - unibases[u].isize( ) 
                              << ", partner = " 
                              << pairs.getPartnerID(id1) << "\n";    }    }
               vec<int> ids;
               for ( size_t k = i; k < j; k++ )
                    ids.push_back( XALIGNS[k].first );
               Sort(ids);
               for ( int r = 0; r < ids.isize( ); r++ )
               {    int s = ids.NextDiff(r);
                    if ( s - r > 1 && VERBOSE )
                    {    cout << "Warning: read " << ids[r] << " placed "
                              << s - r << " times\n";    }
                    r = s - 1;    }
               UniqueSort(ids);
               vec< pair<int,int> > IDS;
               if (VERBOSE) cout << "\nusing:\n";
               for ( int x = 0; x < ids.isize( ); x++ )
               {    if (VERBOSE)
                    {    cout << "[" << x << "] " << ids[x] << " "
                              << pairs.getPartnerID( ids[x] ) << "\n";    }
                    IDS.push( ids[x], pairs.getPartnerID( ids[x] ) );    }
               if (VERBOSE) cout << "total reads placed = " << ids.size( ) << "\n";
               BuildClosuresFromReads( IDS, bases, quals, pairs, evil_K, 
                    top_period, data_dir, CLOSURES_u[u], VERBOSE );    }
          for ( size_t u = 0; u < unibases.size( ); u++ )
               CLOSURES.append( CLOSURES_u[u] );    }

     // Trim closures.  To keep an end of a closure, it must end in a (K-1)-mer
     // present in the unibases or another closure.

     UniqueSort(CLOSURES);
     int evil_KM = evil_K - 1;
     vec<basevector> uni;
     for ( size_t i = 0; i < unibases.size( ); i++ )
     {    const basevector& u = unibases[i];
          for ( int j = 0; j <= u.isize( ) - evil_KM; j++ )
               uni.push( u, j, evil_KM );    }
     UniqueSort(uni);
     vec<basevector> CK;
     for ( size_t i = 0; i < CLOSURES.size( ); i++ )
     {    const basevector& u = CLOSURES[i];
          for ( int j = 0; j <= u.isize( ) - evil_K; j++ )
               CK.push( u, j, evil_K );    }
     UniqueSort(CK);
     vec< vec<int> > NEXTS, BACKS;
     {    GetNexts( evil_K, CK, NEXTS );
          vec<basevector> CKR(CK);
          for ( int i = 0; i < CKR.isize( ); i++ )
               CKR[i].ReverseComplement( );
          GetNexts( evil_K, CKR, BACKS );    }
     vec<Bool> del_right( CK.size( ), False ), del_left( CK.size( ), False );
     while(1)
     {    int dels = 0;
          for ( int i = 0; i < CK.isize( ); i++ )
          {    if ( del_right[i] ) continue;
               Bool OK = False;
               for ( int j = 0; j < NEXTS[i].isize( ); j++ )
                    if ( !del_right[ NEXTS[i][j] ] ) OK = True;
               if ( !OK )
               {    basevector b( CK[i], 1, evil_K-1 );
                    if ( !BinMember( uni, b ) )
                    {    del_right[i] = True;
                         dels++;    }    }    }
          if ( dels == 0 ) break;    }
     while(1)
     {    int dels = 0;
          for ( int i = 0; i < CK.isize( ); i++ )
          {    if ( del_left[i] ) continue;
               Bool OK = False;
               for ( int j = 0; j < BACKS[i].isize( ); j++ )
                    if ( !del_left[ BACKS[i][j] ] ) OK = True;
               if ( !OK )
               {    basevector b( CK[i], 0, evil_K-1 );
                    if ( !BinMember( uni, b ) )
                    {    del_left[i] = True;
                         dels++;    }    }    }
          if ( dels == 0 ) break;    }
     for ( size_t i = 0; i < CLOSURES.size( ); i++ )
     {    basevector& c = CLOSURES[i];
          int j;
          for ( j = 0; j <= c.isize( ) - evil_K; j++ )
          {    basevector d( c, j, evil_K );
               if ( !del_left[ BinPosition( CK, d ) ] ) break;    }
          c.SetToSubOf( c, j, c.isize( ) - j );
          for ( j = c.isize( ) - evil_K; j >= 0; j-- )
          {    basevector d( c, j, evil_K );
               if ( !del_right[ BinPosition( CK, d ) ] ) break;    }
          c.resize( c.isize( ) - ( c.isize( ) - evil_K - j ) );
          if ( c.isize( ) < evil_K ) c.resize(0);    }
     UniqueSort(CLOSURES);

     // Report results about all closures and rebuild unipaths.

     cout << "\ntotal closures = " << CLOSURES.size( ) << endl;
     cout << Date( ) << ": " << unibases.size( ) << " old unibases" << endl;
     int nuni = unibases.size( );
     vec< vec<int> > nexts;
     GetNexts( K, unibases, nexts );
     for ( int u = 0; u < nuni; u++ )
     {    for ( int j = 0; j < nexts[u].isize( ); j++ )
          {    basevector b1 = unibases[u], b2 = unibases[ nexts[u][j] ];
               b1.SetToSubOf( b1, b1.isize( ) - K, K );
               b2.SetToSubOf( b2, K-1, 1 );
               basevector b = Cat( b1, b2 );
               unibases.push_back_reserve(b);    }    }
     size_t N = unibases.size( );
     for ( size_t j = 0; j < CLOSURES.size( ); j++ )
     {    if ( CLOSURES[j].size( ) == 0 ) continue;
          unibases.push_back_reserve( CLOSURES[j] );
          CLOSURES[j].ReverseComplement( );
          unibases.push_back_reserve( CLOSURES[j] );    }
     vecKmerPath paths, pathsrc, unipaths;
     vec<tagged_rpint> pathsdb, unipathsdb;
     const int threads_to_use = 2 * omp_get_max_threads( ) / 3;
     Mkdir777( run_dir + "/tmp" );
     ReadsToPathsCoreY( unibases, K, paths, pathsrc, 
          pathsdb, run_dir + "/tmp", threads_to_use );
     Unipath( paths, pathsrc, pathsdb, unipaths, unipathsdb );
     digraph A;
     BuildUnipathAdjacencyGraph( paths, pathsrc, pathsdb, unipaths, unipathsdb, A );
     HyperKmerPath h;
     BuildUnipathAdjacencyHyperKmerPath( K, A, unipaths, h );
     KmerBaseBroker kbb( K, paths, pathsrc, pathsdb, unibases );
     unibases.clear( );
     for ( int i = 0; i < h.EdgeObjectCount( ); i++ )
          unibases.push_back( kbb.Seq( h.EdgeObject(i) ) );
     cout << Date( ) << ": " << unibases.size( ) << " new unibases" << endl;

     // Write new files.

     unibases.WriteAll( run_dir + "/" + UNIBASES_OUT + ".unibases.k" + KS );
     paths.WriteAll( run_dir + "/" + UNIBASES_OUT + ".paths.k" + KS );
     pathsrc.WriteAll( run_dir + "/" + UNIBASES_OUT + ".paths_rc.k" + KS );
     unipaths.WriteAll( run_dir + "/" + UNIBASES_OUT + ".unipaths.k" + KS );
     BinaryWrite3( run_dir + "/" + UNIBASES_OUT + ".pathsdb.k" + KS, pathsdb );
     BinaryWrite3( run_dir + "/" + UNIBASES_OUT + ".unipathsdb.k" + KS,
          unipathsdb );    }

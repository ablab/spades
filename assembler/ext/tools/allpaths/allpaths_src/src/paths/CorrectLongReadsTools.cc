///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "Basevector.h"
#include "CoreTools.h"
#include "PackAlign.h"
#include "ParallelVecUtilities.h"
#include "PrintAlignment.h"
#include "Qualvector.h"
#include "graph/DigraphTemplate.h"
#include "math/HoInterval.h"
#include "pairwise_aligners/SmithWatBandedA.h"
#include "paths/CorrectLongReadsTools.h"
#include "paths/CorrectPatch.h"
#include "paths/CorrectPatch3.h"
#include "paths/KmerAlignSet.h"
#include "paths/LongReadPatchOptimizer.h"
#include "paths/LongReadTools.h"

// FollowVotes.
// - ppaths and votes are in bijective correspondence
// - low and high define a range of entries in ppaths

void FollowVotes( int i, int low, int high,
     const int np, const vec< vec<int> >& ppaths, const vec<double>& votes,
     const vec< pair<int,int> >& verts, const vecbasevector& unibases, const int K, 
     const heuristics& heur, ostream& hout,
     vec< vec<int> >& answer, vec<double>& answer_votes )
{    
     // Tabulate calls.

     int orig_low = low;
     hout << "\n";
     PRINT3_TO( hout, i, low, high );
     double vote_total = 0.0;
     for ( int j = low; j < high; j++ )
          vote_total += votes[j];
     if ( vote_total < heur.min_votes ) return;
     while( low < high && i >= ppaths[low].isize( ) ) low++;
     vec< triple<int,double,int> > calls;
     vec< triple<double,int,int> > calls2;
     for ( int j = low; j < high; j++ )
     {    if ( i < ppaths[j].isize( ) )
               calls.push( ppaths[j][i], votes[j], j );    }
     Sort(calls);
     for ( int j = 0; j < calls.isize( ); j++ )
     {    int k;
          for ( k = j + 1; k < calls.isize( ); k++ )
               if ( calls[k].first != calls[j].first ) break;
          double sum = 0.0;
          for ( int r = j; r < k; r++ )
               sum += calls[r].second;
          calls2.push( sum, calls[j].first, calls[j].third );
          j = k - 1;    }
     ReverseSort(calls2);

     // If the winning call defines a full-length overlap, delete any calls 
     // that do not.

     if ( i > 0 && calls2.nonempty( ) )
     {    int id0 = ppaths[low][i-1], id1 = calls2[0].second;
          int stop0 = verts[id0].first + unibases[ verts[id0].second ].isize( );
          vec<Bool> to_delete( calls2.size( ), False );
          if ( stop0 - verts[id1].first == K - 1 )
          {    for ( int r = 1; r < calls2.isize( ); r++ )
               {    int id2 = calls2[r].second;
                    if ( stop0 - verts[id2].first < K - 1 )
                         to_delete[r] = True;    }    }
          EraseIf( calls2, to_delete );    }

     // Announce results.

     hout << "position " << i+1 << endl;
     for ( int r = 0; r < calls2.isize( ); r++ )
     {    if ( calls2[r].first < calls2[0].first / 4.0 ) break;
          hout << "[" << calls2[r].second << "], votes = " << calls2[r].first;
          if ( r == 0 ) hout << " (WINNER)";
          hout << "\n";    }

     // Explore the alternatives.

     Bool found = False;
     for ( int j = 0; j < calls2.isize( ); j++ )
     {    Bool Protected = ( heur.min_votes_to_protect >= 0 &&
               calls2[j].first >= heur.min_votes_to_protect );
          if ( !Protected && ( calls2[j].first < heur.min_votes
               || calls2[0].first >= heur.min_win_ratio * calls2[j].first ) ) 
          {    break;    }
          int winner = calls2[j].second, LOW(low), HIGH(high);
          while(1)
          {    if ( i >= ppaths[LOW].isize( ) || ppaths[LOW][i] != winner ) LOW++;
               else break;    }
          while(1)
          {    if ( ppaths[HIGH-1][i] != winner ) HIGH--;
               else break;    }
          if ( i == np - 1 ) break;
          found = True;
          FollowVotes( i+1, LOW, HIGH, np, ppaths, votes, verts, unibases, K, 
               heur, hout, answer, answer_votes );    }
     if (found) return;

     // Save.

     int w = ( low < high ? low : orig_low );
     if ( calls2.nonempty( ) ) w = calls2[0].third;
     vec<int> p = ppaths[w];

     if ( i == np - 1 && low < high ) p.resize(i+1);
     else if ( i > 0 ) p.resize(i);
     else p.resize(1);

     answer.push_back(p);
     answer_votes.push_back(vote_total);
     hout << "SAVING ";
     for ( int l = 0; l < p.isize( ); l++ )
     {    if ( l > 0 ) hout << " --> ";
          hout << "[" << p[l] << "]";    }
     hout << "; votes = " << vote_total << "\n";    }

void ComputeMatches( 
     
     // inputs:

     const vec<int>& p,                            // a vertex path in the graph G
     const vecbasevector& unibases,
     const vec< pair<int,int> >& verts,            // vertices for the graph G
                                                   // (position p, unipath v)
     const vec<int>& unis,                         // unique-sorted verts.second
     const vec< vec< pair<int,int> > >& umatches,  // matches of read to unibases
     
     // logging:

     ostream& out, const Bool PRINT_MATCHES, 

     // outputs:

     vec< pair<int,int> >& matches )
{
     // Compute matches.  

     vec<int> starts( p.size( ) );
     basevector tx = unibases[ verts[ p[0] ].second ];
     starts[0] = 0;
     for ( int l = 1; l < p.isize( ); l++ )
     {    int o = verts[ p[l-1] ].first + unibases[ verts[ p[l-1] ].second ].isize( )
               - verts[ p[l] ].first;
          tx.resize( tx.isize( ) - o );
          starts[l] = tx.size( );
          tx = Cat( tx, unibases[ verts[ p[l] ].second ] );    }
     for ( int l = 0; l < p.isize( ); l++ )
     {    int uj = BinPosition( unis, verts[ p[l] ].second );
          vec< pair<int,int> > m = umatches[uj];
          for ( int j = 0; j < m.isize( ); j++ )
               m[j].second += starts[l];
          matches.append(m);    }
     UniqueSort(matches);

     // Filter matches.   First do a linear regression (based on a reduced match set
     // after removing outliers and points for which the read position appears more 
     // than once).  Then filter again based on the results of the linear regression.

     const int max_delta_offset = 200;
     const double max_delta_tpos = 200.0;
     vec< pair<int,int> > matches2(matches);
     int N2 = matches2.size( );
     vec<int> offset;
     for ( int j = 0; j < N2; j++ )
          offset.push_back( matches2[j].second - matches2[j].first );
     Sort(offset);
     if ( offset.empty( ) ) return; // *********************************************
     int median_offset = offset[ offset.size( ) / 2 ];
     vec<Bool> to_delete( N2, False );
     for ( int j = 0; j < N2; j++ )
     {    int offset = matches2[j].second - matches2[j].first;
          if ( Abs( offset - median_offset ) > max_delta_offset )
               to_delete[j] = True;
          if ( j > 0 && matches2[j].first == matches2[j-1].first )
               to_delete[j] = True;
          if ( j < N2 - 1 && matches2[j].first == matches2[j+1].first )
               to_delete[j] = True;    }
     EraseIf( matches2, to_delete );
     N2 = matches2.size( );

     // Compute the linear regression.

     double rpos_mean = 0.0, tpos_mean = 0.0;
     for ( int j = 0; j < N2; j++ )
     {    rpos_mean += matches2[j].first, tpos_mean += matches2[j].second;    }
     rpos_mean /= double(N2), tpos_mean /= double(N2);
     double rv = 0.0, mv = 0.0;
     for ( int j = 0; j < N2; j++ )
     {    rv += (matches2[j].first - rpos_mean) * (matches2[j].first - rpos_mean);
          mv += matches2[j].second * (matches2[j].first - rpos_mean);    }
     double slope = mv/rv;
     int N = matches.size( );
     vec<double> delta_tpos(N);
     for ( int j = 0; j < N; j++ )
     {    double tpos_predicted = tpos_mean 
               + slope * ( matches[j].first - rpos_mean );
          delta_tpos[j] = matches[j].second - tpos_predicted;    }

     // Filter versus the linear regression.

     to_delete.resize_and_set( N, False );
     for ( int j = 0; j < N; j++ )
     {    if ( Abs( delta_tpos[j] ) > max_delta_tpos )
               to_delete[j] = True;    }
     EraseIf( matches, to_delete );
     EraseIf( delta_tpos, to_delete );

     // In preparation for more filtering, form the rolling mean over windows.  For 
     // this purpose, again exclude points for which the read position appears more 
     // than once.

     const int roll_window = 100;
     vec< pair<int,int> > matches3(matches);
     to_delete.resize_and_set( matches3.size( ), False );
     for ( int j = 0; j < matches3.isize( ); j++ )
     {    if ( j > 0 && matches3[j].first == matches3[j-1].first )
               to_delete[j] = True;
          if ( j < matches3.isize( ) - 1 && matches3[j].first 
               == matches3[j+1].first )
          {    to_delete[j] = True;    }    }
     EraseIf( matches3, to_delete );
     vec<double> delta_tpos3(delta_tpos);
     EraseIf( delta_tpos3, to_delete );

     // Get the rolling mean.

     int N3 = matches3.size( );
     vec<double> roll_mean3( N3, 0.0 );
     const int wpoints_goal = 20;
     int wpoints = Min( N3, wpoints_goal );
     double sum = 0.0, count = 0.0;
     for ( int j = 0; j < wpoints; j++ )
     {    count++;
          sum += delta_tpos3[j];    }
     vec<Bool> computed( matches.size( ), False );
     vec<double> roll_mean( matches.size( ), 0.0 );
     if ( count > 0 ) 
     {    roll_mean3[0] = sum / count;
          for ( int j = wpoints; j < N3; j++ )
          {    sum += delta_tpos3[j] - delta_tpos3[j-wpoints];
               roll_mean3[j-wpoints+1] = sum/count;    }
          for ( int j = N3 - wpoints + 1; j < N3; j++ )
               roll_mean3[j] = roll_mean3[ N3 - wpoints ];
          int mpos = 0;
          for ( int j = 0; j < matches3.isize( ); j++ )
          {    while( matches[mpos] != matches3[j] ) mpos++;
               roll_mean[mpos] = roll_mean3[j];
               computed[mpos] = True;    }
          // The following loop is unnecessarily quadratic.
          for ( int j = 0; j < matches.isize( ); j++ )
          {    int k = j;
               while( k < matches.isize( ) && !computed[k] ) k++;
               if ( k < matches.isize( ) ) 
                    roll_mean[j] = roll_mean[k];
               else 
               {    if ( j == 0 ) // not sure what to do here
                         roll_mean[j] = 0.0;
                    else roll_mean[j] = roll_mean[j-1];    }    }

          // Now use the rolling mean to filter.

          const int max_dist = 20;
          to_delete.resize_and_set( matches.size( ), False );
          for ( int j = 0; j < matches.isize( ); j++ )
          {    if ( Abs( delta_tpos[j] - roll_mean[j] ) > max_dist )
                    to_delete[j] = True;    }
          EraseIf( matches, to_delete );
          EraseIf( delta_tpos, to_delete );
          EraseIf( roll_mean, to_delete );     }

     // Filter out mult-hits.

     to_delete.resize_and_set( matches.size( ), False );
     for ( int j = 0; j < matches.isize( ); j++ )
     {    if ( j > 0 && matches[j].first == matches[j-1].first ) to_delete[j] = True;
          if ( j < matches.isize( ) - 1 && matches[j].first == matches[j+1].first )
               to_delete[j] = True;    }
     EraseIf( matches, to_delete );
     EraseIf( delta_tpos, to_delete );
     EraseIf( roll_mean, to_delete );
               
     // Print matches.

     if (PRINT_MATCHES)
     {    out << "\n\n";
          for ( int j = 0; j < matches.isize( ); j++ )
          {    out << "new match: " << matches[j].first << " --> " 
                    << matches[j].second << "   " << setiosflags(ios::fixed) 
                    << setprecision(1) << delta_tpos[j] << "   " << roll_mean[j]
                    << resetiosflags(ios::fixed) << "\n";    } 
          out << "\n";    }    }

void AlignReadToPath( 

     // inputs:

     const basevector& r, const int id, const Bool fw, const vec<int>& p, 
     const int path_id, const vecbasevector& unibases,
     const vec< pair<int,int> >& verts, const vec<int>& unis,
     const vec< vec< pair<int,int> > >& umatches, const vec< vec<int> >& ppaths,
     vec<Bool>& ppath_seen, const int L, const Bool PRINT_MATCHES,
     const Bool DUMP_LOCAL, const Bool PRINT_LM1, 

     // outputs:

     double& ERRS, int& LMATCHES,
     triple<int,int,String>& EREPORTS, int& LASTPOS, int& LAST, Bool& SKIPPED,
     map< pair<ho_interval,basevector>, pair<align,double> >& A )
{    
     ostringstream out;
     out << "read " << ( !fw ? "-" : "" ) << id << " vs path " << path_id << " = ";
     for ( int l = 0; l < p.isize( ); l++ )
     {    if ( l > 0 ) out << " --> ";
          out << "[" << p[l] << "]";    }
     out << "; errs = ";

     basevector t = unibases[ verts[ p[0] ].second ];
     vec<int> stops( p.size( ) );
     stops[0] = t.size( );
     for ( int l = 1; l < p.isize( ); l++ )
     {    int o = verts[ p[l-1] ].first + unibases[ verts[ p[l-1] ].second ].isize( )
               - verts[ p[l] ].first;
          t.resize( t.isize( ) - o );
          t = Cat( t, unibases[ verts[ p[l] ].second ] );
          stops[l] = t.size( );    }

     vec< pair<int,int> > matches;
     vec<ho_interval> LM1, LM2;
     ComputeMatches(p, unibases, verts, unis, umatches, out, PRINT_MATCHES, matches);
     const int infinity = 1000000000;
     if ( matches.empty( ) ) // shouldn't ever happen
     {    LAST = -1;
          ERRS = infinity;
          LMATCHES = -infinity;    }
     else
     {    int tlast = matches.back( ).second;
          LAST = 0;
          while( LAST < p.isize( ) - 1 && tlast > stops[LAST] ) LAST++;

          // A given partial path should only be processed once.

          vec<int> px(p);
          px.resize( LAST + 1 );
          int pp = BinPosition( ppaths, px );
          ForceAssertGe( pp, 0 );
          if ( ppath_seen[pp] ) 
          {    SKIPPED = True;
               return;    }
          ppath_seen[pp] = True;

          // Note a key problem here: we're going to assess path p based on it
          // mapping to a certain region of the read, and different ps will map to
          // different regions, making the comparison in its current form unfair.

          basevector rx, tx;
          align a;
          for ( int j = 0; j < matches.isize( ); j++ )
          {    int k, rpos1, tpos1, rpos2 = -1, tpos2 = -1;
               rpos1 = matches[j].first, tpos1 = matches[j].second;
               LM1.push( rpos1, rpos1 + L ), LM2.push( tpos1, tpos1 + L );
               for ( k = j + 1; k < matches.isize( ); k++ )
               {    rpos2 = matches[k].first, tpos2 = matches[k].second;
                    if ( rpos2 == rpos1 || tpos2 <= tpos1 ) continue;
                    if ( rpos2 <= rpos1 + L || tpos2 <= tpos1 + L ) continue;
                    break;    }
               if ( k == matches.isize( ) ) 
               {    for ( int z = j + 1; z < matches.isize( ); z++ )
                    {    rpos2 = matches[z].first, tpos2 = matches[z].second;
                         if ( rpos2 - rpos1 == tpos2 - tpos1 )
                         {    LM1.push( rpos2, rpos2 + L );
                              LM2.push( tpos2, tpos2 + L );    }    }
                    break;    }
               rx.SetToSubOf( r, rpos1 + L, rpos2 - (rpos1 + L ) );
               tx.SetToSubOf( t, tpos1 + L, tpos2 - (tpos1 + L ) );

               // Align, or fetch cached copy.
                              
               ho_interval h( rpos1 + L, rpos2 );
               map< pair<ho_interval,basevector>, pair<align,double> > 
                    ::iterator Ax = A.find( make_pair( h, tx ) );
               if ( Ax == A.end( ) )
               {    double err = SmithWatFreeSym( rx, tx, a, True, True, 1, 1 );    
                    ERRS += err;
                    A[ make_pair( h, tx ) ] = make_pair( a, err );    }
               else
               {    a = Ax->second.first;
                    ERRS += Ax->second.second;    }

               // Print alignment.

               if (DUMP_LOCAL)
               {    out << "\nlocal alignment:\n";
                    PRINT4_TO( out, rpos1+L, rpos2, tpos1+L, tpos2 );
                    PrintVisualAlignment( False, out, rx, tx, a );    }

               // Add to LM1, LM2.

               vec<ho_interval> perfs1, perfs2;
               a.PerfectIntervals1( rx, tx, perfs1 );
               a.PerfectIntervals2( rx, tx, perfs2 );
               for ( int l = 0; l < perfs1.isize( ); l++ )
               {    perfs1[l].Shift(rpos1+L), perfs2[l].Shift(tpos1+L);     }
               LM1.append(perfs1), LM2.append(perfs2);
               LM1.push( rpos2, rpos2 + L ), LM2.push( tpos2, tpos2 + L );

               LASTPOS = rpos2 + L;
               j = k - 1;    }    }

     // Condense LM1, LM2 and print.

     if ( LM1.nonempty( ) )
     {    vec<ho_interval> LM1_new, LM2_new;
          LM1_new.push_back( LM1[0] ), LM2_new.push_back( LM2[0] );
          for ( int l = 1; l < LM1.isize( ); l++ )
          {    if ( LM1[l].Start( ) <= LM1_new.back( ).Stop( )
                    && LM2[l].Start( ) <= LM2_new.back( ).Stop( )
                    && LM1_new.back( ).Stop( ) - LM1[l].Start( )
                    == LM2_new.back( ).Stop( ) - LM2[l].Start( ) )
               {    LM1_new.back( ).SetStop( LM1[l].Stop( ) );
                    LM2_new.back( ).SetStop( LM2[l].Stop( ) );    }
               else
               {    LM1_new.push_back( LM1[l] );
                    LM2_new.push_back( LM2[l] );    }    }
          LM1 = LM1_new, LM2 = LM2_new;
          LM1_new.clear( ), LM2_new.clear( );
          for ( int l = 0; l < LM1.isize( ); l++ )
          {    if ( LM1[l].Length( ) >= L ) 
               {    LM1_new.push_back( LM1[l] );
                    LM2_new.push_back( LM2[l] );    }    }
          LM1 = LM1_new, LM2 = LM2_new;    }
     if (PRINT_LM1)
     {    out << "\nmatching intervals:\n";
          for ( int l = 0; l < LM1.isize( ); l++ )
          {    if ( LM1[l].Length( ) >= L )
                    out << LM1[l] << " --> " << LM2[l] << "\n";    }    }

     LMATCHES = 0;
     for ( int l = 0; l < LM1.isize( ); l++ )
          if ( LM1[l].Length( ) >= L ) LMATCHES += LM1[l].Length( ) - L + 1;

     out << setprecision(5) << ERRS << ", Lmatches = " << LMATCHES 
          << ", lastpos = " << LASTPOS << ", last = " << LAST << "\n";    

     EREPORTS = make_triple( LASTPOS, /* ERRS */ LMATCHES, out.str( ) );    }

// Partition reads into similar orbits, discarding those that do not clearly lie
// in a single orbit.  We proceed by left trimming the uniseqs (removing stuff 
// before u), then finding the top 40 by length, trimming them to the minimum
// of these lengths, then aligning them to each other using a banded Smith-Waterman
// with bandwidth 10.

void PartitionReads( const int u, const vec< triple<int,Bool,uniseq> >& rights,
     vec< vec<int> >& orbits, ostringstream& hout )
{
     double eclock = WallClockTime( );

     // Heuristics.

     const int max_align = 40;
     const int bandwidth = 10;
     const int max_dist = 10;
     const float score_mult = 3.0;
     const float score_diff = 5.0;

     // Do left trimming.

     vec<basevector> rbases;
     int nr = rights.size( );
     vec<int> rsizes, ids( nr, vec<int>::IDENTITY );
     for ( int i = 0; i < nr; i++ )
     {    uniseq q = rights[i].third;
          int n;
          for ( n = 0; n < q.N( ); n++ )
               if ( q.U(n) == u ) break;
          q.TrimLeft(n);
          rbases.push_back( q.Bases( ) );
          rsizes.push_back( rbases.back( ).size( ) );    }
     ReverseSortSync( rsizes, rbases, ids );
     int ns = Min( max_align, nr );

     // Trim all sequences to have the same length, then sort again, then determine
     // which are equal.

     vec<basevector> rbasesx(ns);
     vec<int> idsx(ns);
     for ( int i = 0; i < ns; i++ )
     {    rbasesx[i] = rbases[i];
          rbasesx[i].resize( rbases[ns-1].size( ) );
          idsx[i] = ids[i];    }
     SortSync( rbasesx, idsx );
     equiv_rel e(ns);
     for ( int i = 0; i < ns; i++ )
     {    int j = rbasesx.NextDiff(i);
          for ( int k = i + 1; k < j; k++ )
               e.Join( k, i );
          j = i - 1;    }

     // Now align.  Use the known equality to avoid excess aligning.

     vec< vec<int> > score( ns, vec<int>(ns, 0) );
     align a;
     int errors;
     for ( int i1 = 0; i1 < ns; i1++ )
     for ( int i2 = 0; i2 < ns; i2++ )
     {    if ( !e.Representative(i1) || !e.Representative(i2) ) continue;
          float s = 2.0 * SmithWatBandedA( 
               rbasesx[i1], rbasesx[i2], 0, bandwidth, a, errors, 0, 1, 1 );
          vec<int> o1, o2;
          e.Orbit( i1, o1 );
          e.Orbit( i2, o2 );
          for ( int j1 = 0; j1 < o1.isize( ); j1++ )
          for ( int j2 = 0; j2 < o2.isize( ); j2++ )
          {    score[ o1[j1] ][ o2[j2] ] = int(round(s));    }    }

     // Join close things.
          
     for ( int i1 = 0; i1 < ns; i1++ )
     for ( int i2 = 0; i2 < ns; i2++ )
     {    if ( score[i1][i2] <= max_dist ) e.Join( i1, i2 );    }
     vec<int> reps;
     e.OrbitRepsAlt(reps);
     hout << "\nlengths:";
     for ( int i = 0; i < ns; i++ )
          hout << " " << rsizes[i];
     hout << "\n\n";
     for ( int i = 0; i < reps.isize( ); i++ )
     {    hout << "orbit " << i << ":";
          vec<int> o;
          e.Orbit( reps[i], o );
          for ( int j = 0; j < o.isize( ); j++ )
               hout << " " << o[j];
          hout << "\n";    }
     hout << "\n";
     for ( int i1 = 0; i1 < reps.isize( ); i1++ )
     for ( int i2 = i1; i2 < reps.isize( ); i2++ )
     {    vec<int> o1, o2;
          e.Orbit( reps[i1], o1 );
          e.Orbit( reps[i2], o2 );
          hout << "distances between orbits " << i1 << " and " << i2 << ":";
          for ( int j1 = 0; j1 < o1.isize( ); j1++ )
          for ( int j2 = 0; j2 < o2.isize( ); j2++ )
          {    int x1 = o1[j1], x2 = o2[j2];
               if ( i1 != i2 || x1 < x2 ) hout << " " << score[x1][x2];    }
          hout << "\n";    }

     // Define initial orbits.

     for ( int j = 0; j < reps.isize( ); j++ )
     {    vec<int> o, ox;
          e.Orbit( reps[j], o );
          for ( int i = 0; i < o.isize( ); i++ )
               ox.push_back( idsx[ o[i] ] );
          orbits.push_back(ox);    }
     hout << "\ntime used so far = " << TimeSince(eclock) << endl;

     // Now figure out what to do with the rest of the reads.

     for ( int i = ns; i < nr; i++ )
     {    if ( orbits.solo( ) ) continue;
          else
          {    hout << "\n";
               vec<float> s;
               for ( int j = 0; j < reps.isize( ); j++ )
               {    hout << "read " << i << " versus orbit " << j << ": ";
                    vec<int> o;
                    e.Orbit( reps[j], o );
                    float errs = 0.0;
                    for ( int k = 0; k < o.isize( ); k++ )
                    {    errs += 2.0 * SmithWatBandedA( rbases[i], rbasesx[ o[k] ], 
                              0, bandwidth, a, errors, 0, 1, 1 );    }
                    errs /= double( o.size( ) );
                    s.push_back(errs);
                    hout << setiosflags(ios::fixed) << setprecision(1) << errs
                         << resetiosflags(ios::fixed) << "\n";    }    
               vec<int> oids( reps.size( ), vec<int>::IDENTITY );
               SortSync( s, oids );
               if ( s.size( ) >= 2 && s[0] * score_mult >= s[1] ) continue;
               if ( s.size( ) >= 2 && s[0] + score_diff >= s[1] ) continue;
               if ( orbits[ oids[0] ].isize( ) == max_align ) continue;
               orbits[ oids[0] ].push_back( ids[i] );
               hout << "adding to orbit " << oids[0] << "\n";    }    }
     hout << "\ntotal reads = " << rights.size( ) << "\n";
     for ( int i = 0; i < orbits.isize( ); i++ )
          hout << "orbit " << i << " has size " << orbits[i].size( ) << "\n";
     hout << "\ntime used = " << TimeSince(eclock) << endl;
     /*
     cout << hout.str( );
     flush(cout);
     _exit(0);
     */
          }

void ScoreGraph( const digraphE<int>& G, const int source,
     const vec< pair< int, vec< pair<int,int> > > >& ALIGNS,
     int& matches, int& errs 
          , vec<int>& path
          , vec<int>& upath // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          )
{
     vec<int> sources, sinks, to_left, to_right;
     G.Sources(sources), G.Sinks(sinks);
     G.ToLeft(to_left), G.ToRight(to_right);
     vec< vec<int> > successors( sources.size( ) );
     for ( int i = 0; i < sources.isize( ); i++ )
     {    G.GetSuccessors1( sources[i], successors[i] );
          Sort( successors[i] );    }
     vec< pair<int,int> > matches_errs;
     vec< vec<int> > paths, upaths;
     for ( int i1 = 0; i1 < sources.isize( ); i1++ )
     for ( int i2 = 0; i2 < sinks.isize( ); i2++ )
     {    if ( sources[i1] != source ) continue;
          if ( !BinMember( successors[i1], sinks[i2] ) ) continue;
          vec<int> epath;
          vec<int> path, upath;
          G.ShortestPath( source, sinks[i2], path );
          int errs = 0;
          for ( int m = 0; m < path.isize( ) - 1; m++ )
          {    int x1 = path[m], x2 = path[m+1];
               errs += G.EdgeObject( G.EdgesBetween( x1, x2 )[0] );    }
          for ( int j = 0; j < path.isize( ) - 1; j++ )
          {    int x = path[j], y = path[j+1];
               int min_edge = 1000000000;
               for ( int l = 0; l < G.From(x).isize( ); l++ )
               {    if ( G.From(x)[l] == y )
                    {    min_edge = Min( min_edge, G.
                              EdgeObjectByIndexFrom( x, l ) );    }    }
               for ( int l = 0; l < G.From(x).isize( ); l++ )
               {    if ( G.From(x)[l] == y )
                    {    if ( G.EdgeObjectByIndexFrom( x, l ) == min_edge )
                         {    epath.push_back( 
                                   G.EdgeObjectIndexByIndexFrom( x, l ) );
                              break;   }    }    }    }
          vec<int> M;
          for ( int r = 0; r < epath.isize( ); r++ )
          {    int e = epath[r];
               int v = to_left[e], w = to_right[e];
               for ( int m = 0; m < ALIGNS[v].second.isize( ); m++ )
                    M.push_back( ALIGNS[v].second[m].first );
                    for ( int m = 0; m < ALIGNS[w].second.isize( ); m++ )
                    M.push_back( ALIGNS[w].second[m].first );     }
          UniqueSort(M);
          for ( int j = 0; j < path.isize( ); j++ ) // XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
               upath.push_back( ALIGNS[ path[j] ].first ); // XXXXXXXXXXXXXXXXXXXXXX
          matches_errs.push( -M.isize( ), errs );
          paths.push_back(path);
          upaths.push_back(upath);    }
     ReverseSortSync( matches_errs, paths, upaths );
     if ( matches_errs.empty( ) ) //  ??????????????????????????????????????????????
     {    matches = -1, errs = -1;
          cout << "matches_errs is empty; weird!\n";
          return;    }
     ForceAssert( matches_errs.nonempty( ) ); // ***********************************
     matches = -matches_errs.back( ).first, errs = matches_errs.back( ).second;
     path = paths.back( );
     upath = upaths.back( );    }

void Phase2( const int u, const vecbasevector& unibases, 
     const vec< vec< pair<int,int> > >& Ulocs, const vec<int>& to_rc,
     const vecbasevector& longreads, const int K, const vec< vec<int> >& hits, 
     const vec<uniseq>& UNISEQ, const vec< vec<int> >& UNISEQ_ID, ostream& hout, 
     const Bool DOT2, const Bool PRINT_MATCHES, const Bool PRINT_LM1, 
     const Bool DUMP_LOCAL, 
     const Bool PRINT_DISCARDS, const Bool VALIDATE1, const String& data_dir, 
     const vecbasevector& genome2, const vec< vec< pair<int,int> > >& Glocs, 
     const int LG, const heuristics& heur, vecbasevector& new_stuff,
     vec< vec<int> >& right_exts,
     const Bool NEW_FILTER, const Bool ANNOUNCE, const int SUPER_VERBOSITY,
     const Bool PRINT_BASES2, const int MIN_SPREAD2 )

{
     hout << "\n=========================================================="
          << "==========================\n";
     double clock = WallClockTime( );
     hout << "\n" << Date( ) << ": LOOKING RIGHT FROM UNIBASE " << u 
          << " (L=" << unibases[u].size( ) << ")\n" << endl;
     vec< triple<int,Bool,uniseq> > rights;

     // Process hits to u.

     for ( int j = 0; j < hits[u].isize( ); j++ )
     {    int id = hits[u][j];
          const uniseq& q = UNISEQ[id];
          if ( q.U( ).back( ) == u ) continue;
          rights.push( id, True, q );    }

     // Process hits to the reverse complement of u.

     int ru = to_rc[u];
     for ( int j = 0; j < hits[ru].isize( ); j++ )
     {    int id = hits[ru][j];
          uniseq q = UNISEQ[id];
          q.ReverseMe(to_rc);
          if ( q.U( ).back( ) == u ) continue;
          rights.push( id, False, q );    }

     // Note read ids.

     vec<int> all_ids;
     for ( int j = 0; j < rights.isize( ); j++ )
          all_ids.push_back( rights[j].first );
     UniqueSort(all_ids);
     hout << "U2=" << u << " READS_TO_USE=\"{";
     for ( int j = 0; j < all_ids.isize( ); j++ )
     {    if ( j > 0 ) hout << ",";
          hout << all_ids[j];    }
     hout << "}\"\n\n";
     if (ANNOUNCE)
     {    
          #pragma omp critical
          {    cout << "\nstarting U2=" << u << " READS_TO_USE=\"{";
               for ( int j = 0; j < all_ids.isize( ); j++ )
               {    if ( j > 0 ) cout << ",";
                    cout << all_ids[j];    }
               cout << "}\"\n" << endl;    }    }

     // Build a graph G.  The vertices are pairs (p, v), corresponding to unipaths v
     // (to the right of u) and their position p, relative to u, so that a given
     // unipath may appear more than once at different positions.  The edges are the
     // edges that appeared in at least one 'winning edge' for a read.

     hout << Date( ) << ": building graph" << endl;
     for ( int j = 0; j < rights.isize( ); j++ )
     {    hout << "<" << j << "> " << ( !rights[j].second ? "-" : "" ) 
               << rights[j].first << "   " << rights[j].third << "\n";    }
     vec< pair<int,int> > verts, edges;
     for ( int j = 0; j < rights.isize( ); j++ )
     {    const uniseq& q = rights[j].third;
          int pos = -unibases[u].isize( ), x = 0;
          while( q.U(x) != u ) x++;
          for ( int i = x; i < q.N( ); i++ )
          {    verts.push( pos, q.U(i) );
               if ( i == q.N( ) - 1 ) break;
               pos += unibases[ q.U(i) ].size( ) - q.Over(i);    }    }

     // Delete rare vertices.

     hout << Date( ) << ": deleting rare verts" << endl;
     const int max_verts = 400;
     Sort(verts);
     int mult_floor = 1;
     while(1)
     {    vec<Bool> verts_to_delete( verts.size( ), False );
          for ( int i = 0; i < verts.isize( ); i++ )
          {    int j = verts.NextDiff(i);
               if ( j - i >= mult_floor )
               {    for ( int k = i + 1; k < j; k++ )
                         verts_to_delete[k] = True;    }
               else      
               {    for ( int k = i; k < j; k++ )
                         verts_to_delete[k] = True;    }
               i = j - 1;    }
          if ( verts.isize( ) - Sum(verts_to_delete) <= max_verts )
          {    EraseIf( verts, verts_to_delete );
               break;    }
          else mult_floor++;    }
     // UniqueSort(verts);

     hout << Date( ) << ": printing verts" << endl;
     hout << "\nverts:\n\n";
     for ( int j = 0; j < verts.isize( ); j++ )
     {    hout << "[" << j << "] " << verts[j].first << "   " << verts[j].second 
               << " (l=" << unibases[ verts[j].second ].size( ) << ")\n";    }
     for ( int j = 0; j < rights.isize( ); j++ )
     {    vec<int> vert_ids;
          const uniseq& q = rights[j].third;
          int pos = -unibases[u].isize( ), x = 0;
          while( q.U(x) != u ) x++;
          for ( int i = x; i < q.N( ); i++ )
          {    vert_ids.push( BinPosition( verts, make_pair( pos, q.U(i) ) ) );
               if ( i == q.N( ) - 1 ) break;
               pos += unibases[ q.U(i) ].size( ) - q.Over(i);    }
          for ( int l = 0; l < vert_ids.isize( ) - 1; l++ )
          {    if ( vert_ids[l] >= 0 && vert_ids[l+1] >= 0 )
                    edges.push( vert_ids[l], vert_ids[l+1] );    }    }
     UniqueSort(edges);
     vec< vec<int> > from( verts.size( ) ), to( verts.size( ) );
     for ( int i = 0; i < edges.isize( ); i++ )
     {    int v = edges[i].first, w = edges[i].second;
          from[v].push_back(w), to[w].push_back(v);    }
     for ( int x = 0; x < verts.isize( ); x++ )
     {    Sort( from[x] ), Sort( to[x] );    }
     digraph G( from, to );

     // Remove transitive edges.

     hout << Date( ) << ": removing trans edges" << endl;
     vec< pair<int,int> > edges_to_remove;
     vec< vec<int> > successors( G.N( ) );
     for ( int v = 0; v < G.N( ); v++ )
          G.GetSuccessors1( v, successors[v] );
     for ( int x = 0; x < G.N( ); x++ )
     for ( int j1 = 0; j1 < G.From(x).isize( ); j1++ )
     for ( int j2 = 0; j2 < G.From(x).isize( ); j2++ )
     {    int y1 = G.From(x)[j1], y2 = G.From(x)[j2];
          if ( y1 != y2 && BinMember( successors[y1], y2 ) )
               edges_to_remove.push( x, y2 );    }
     UniqueSort(edges_to_remove);
     vec< vec<int> > from2( G.N( ) ), to2( G.N( ) );
     for ( int x = 0; x < G.N( ); x++ )
     {    for ( int j = 0; j < G.From(x).isize( ); j++ )
          {    int y = G.From(x)[j];
               if ( !BinMember( edges_to_remove, make_pair( x, y ) ) )
                    from2[x].push_back(y);    }
          for ( int j = 0; j < G.To(x).isize( ); j++ )
          {    int y = G.To(x)[j];
               if ( !BinMember( edges_to_remove, make_pair( y, x ) ) )
                    to2[x].push_back(y);    }    }
     G = digraph( from2, to2 );

     // Print graph.

     hout << Date( ) << ": printing graph" << endl;
     hout << "\nedges:\n\n";
     for ( int v = 0; v < G.N( ); v++ )
     {    for ( int j = 0; j < G.From(v).isize( ); j++ )
          {    int w = G.From(v)[j];
               hout << "[" << v << "] --> [" << w << "]\n";    }    }
     if (DOT2)
     {    Ofstream( out, "xxx.dot" );
          digraph(G).PrettyDOT( out, False, True );    }

     // Compute frequency of short partial paths.

     const int depth = 3;
     const int min_count = 3;
     vec< vec< vec<int> > > dpaths(depth+1);
     Bool depth_verbose = False;
     if (depth_verbose) hout << "\npartial paths of depth <= " << depth << "\n";
     vec< vec<int> > hpaths;
     for ( int v = 0; v < G.N( ); v++ )
     {    vec<int> p;
          p.push_back(v);
          hpaths.push_back(p);    }
     while( hpaths.nonempty( ) )
     {    vec<int> p = hpaths.back( );
          hpaths.pop_back( );
          dpaths[ p.size( ) ].push_back(p);
          if ( p.isize( ) < depth )
          {    int v = p.back( );
               for ( int j = 0; j < G.From(v).isize( ); j++ )
               {    vec<int> q(p);
                    q.push_back( G.From(v)[j] );
                    hpaths.push_back(q);    }    }    }
     for ( int d = 0; d <= depth; d++ )
          Sort( dpaths[d] );
     vec< vec<int> > dcount(depth+1);
     for ( int d = 1; d <= depth; d++ )
          dcount[d].resize( dpaths[d].size( ) );
     for ( int j = 0; j < rights.isize( ); j++ )
     {    const uniseq& q = rights[j].third;
          int pos = -unibases[u].isize( ), x = 0;
          while( q.U(x) != u ) x++;
          vec<int> v;
          for ( int i = x; i < q.N( ); i++ )
          {    v.push_back( BinPosition( verts, make_pair( pos, q.U(i) ) ) );
               if ( i == q.N( ) - 1 ) break;
               pos += unibases[ q.U(i) ].size( ) - q.Over(i);    }
          for ( int d = 1; d <= depth; d++ )
          {    for ( int l = 0; l <= v.isize( ) - d; l++ )
               {    vec<int> w;
                    w.SetToSubOf( v, l, d );
                    int p = BinPosition( dpaths[d], w );
                    // HOW COULD THIS BE NEGATIVE??
                    // PRINT4( j, d, l, p ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                    if ( p >= 0 ) dcount[d][p]++;    }    }    }
     if (depth_verbose)
     {    for ( int m = 0; m < dpaths[depth].isize( ); m++ )
          {    const vec<int>& p = dpaths[depth][m];
               hout << "[" << m << "]";
               for ( int j = 0; j < p.isize( ); j++ )
                    hout << " " << p[j];
               int start = verts[ p.front( ) ].first;
               int nbases = verts[ p.back( ) ].first 
                    + unibases[ verts[ p.back( ) ].second ].isize( ) - start;
               hout << " (bases=" << nbases << ",start=" 
                    << start << ",count=" << dcount[depth][m] << ")\n";    }    }
     vec< vec<Bool> > acceptable( depth + 1 );
     for ( int d = 1; d <= depth; d++ )
     {    acceptable[d].resize( dcount[d].size( ) );
          for ( int j = 0; j < dcount[d].isize( ); j++ )
          {    acceptable[d][j] = ( dcount[d][j] >= min_count );
               if (depth_verbose)
               {    if ( acceptable[d][j] )
                    {    hout << "accepting";
                         for ( int l = 0; l < dpaths[d][j].isize( ); l++ )
                              hout << " " << dpaths[d][j][l];
                         hout << "\n";    }    }    }    }

     // New filtering.

     vec<int> xsuggested_kills;
     vec< vec<int> > vpaths_final;
     int L = heur.L;
     if (NEW_FILTER)
     {

     // Find matches of reads to unibases.

     hout << Date( ) << ": compute matches for all reads" << endl;
     vec<KmerAlignSet> kaligns( rights.size( ) );
     vec<String> kreport( rights.size( ) );
     #pragma omp parallel for
     for ( int j = 0; j < rights.isize( ); j++ )
     {    Bool fw = rights[j].second;
          int id = rights[j].first;
          basevector r = longreads[id];
          if ( !fw ) r.ReverseComplement( );
          KmerAlignSet alx( r, L, Ulocs ), aly;
          const Bool CLUSTER_ALIGNS_NEW_CLEAN = False;
          ClusterAlignsNew( alx, aly, CLUSTER_ALIGNS_NEW_CLEAN, MIN_SPREAD2 );
          KillInferiorClustersNew( aly, unibases );
          if ( SUPER_VERBOSITY >= 1 )
          {    ostringstream rout;
               rout << "\n" << aly.NAligns( ) << " clustered alignments "
                    << "of read(right) " << id << "(" << j << ")\n";
               for ( int i = 0; i < aly.NAligns( ); i++ )
               {    rout << "\n[" << i << "] clustered alignment of read(right) " 
                         << id << "(" << j << ") to unibase " << aly.U(i) 
                         << "[l=" << unibases[ aly.U(i) ].size( ) << "]:\n";
                    for ( int k = 0; k < aly.Count(i); k++ )
                    {    if ( k > 0 && k % 6 == 0 ) rout << "\n";
                         rout << "(" << aly.Rpos(i,k) << "," << aly.Upos(i,k) 
                              << ")";    }
                    rout << "\n";    }    
               kreport[j] = rout.str( );    }
          kaligns[j] = aly;    }
     if ( SUPER_VERBOSITY >= 1 )
     {    for ( int j = 0; j < rights.isize( ); j++ )
               hout << kreport[j];    }

     // For each read, and each edge, find the pairs of alignments that support it.

     hout << Date( ) << ": finding support" << endl;
     vec< vec< vec< vec< pair<int,int> > > > > support( rights.size( ) );
     for ( int r = 0; r < rights.isize( ); r++ )
     {    support[r].resize( G.N( ) );
          const KmerAlignSet& a = kaligns[r];
          for ( int v1 = 0; v1 < G.N( ); v1++ )
          {    support[r][v1].resize( G.From(v1).size( ) );
               for ( int j = 0; j < G.From(v1).isize( ); j++ )
               {    int v2 = G.From(v1)[j];
                    int u1 = verts[v1].second, u2 = verts[v2].second;
                    const basevector &U1 = unibases[u1], &U2 = unibases[u2];
                    int pos1 = verts[v1].first, pos2 = verts[v2].first;
                    for ( int l1 = 0; l1 < a.NAligns( ); l1++ )
                    {    if ( a.U(l1) != u1 ) continue;
                         for ( int l2 = 0; l2 < a.NAligns( ); l2++ )
                         {    if ( a.U(l2) != u2 ) continue;

                              // Check for compatibility of the alignment to u2.
                              // We look for a single kmer match that validates
                              // this compatibility.
                              //
                              // a1 --> a2
                              // pos1 --> pos2
                              // u1 --> u2
                              // u1 starts at pos1
                              // u2 starts at pos2
                              // a1.second = {(rpos,upos)}
                              // a2.second = {(rpos,upos)}
                              //
                              // pos1
                              // -----------------------|----|-----
                              //                   pos2
                              //                   -----|----|--------------------

                              Bool compatible = False;
                              vec<int> o, ox, ocount;
                              for ( int m1 = 0; m1 < a.Count(l1); m1++ )
                              for ( int m2 = 0; m2 < a.Count(l2); m2++ )
                              {    if ( a.Rpos(l1,m1) != a.Rpos(l2,m2) ) continue;
                                   int upos1 = a.Upos(l1,m1), upos2 = a.Upos(l2,m2);
                                   if ( pos1 + upos1 != pos2 + upos2 ) continue;
                                   o.push_back( upos1 - upos2 );    }
                              if ( o.empty( ) ) continue;

                              // Find the most common offset.

                              Sort(o);
                              for ( int i = 0; i < o.isize( ); i++ )
                              {    int j = o.NextDiff(i);
                                   ox.push_back( o[i] );
                                   ocount.push_back( j - i );
                                   i = j - 1;    }
                              ReverseSortSync( ocount, ox );
                              int off = ox[0];
                              
                              // We require that u2 is a proper right
                              // extension of u1.

                              if ( U1.isize( ) - U2.isize( ) >= off ) continue;

                              // Make sure that the unibases align.

                              Bool mismatch = False;
                              for ( int up1 = Max( 0, off );
                                   up1 < Min( U1.isize( ), 
                                        U2.isize( ) + off ); up1++ )
                              {    int up2 = up1 - off;
                                   if ( U1[up1] != U2[up2] )
                                   {    mismatch = True;
                                        break;    }    }
                              if ( !mismatch ) compatible = True;  
                              if (compatible)
                              {    support[r][v1][j].push( 
                                        l1, l2 );    }    }    }    }    }    }

     // Build successive partial paths.

     hout << Date( ) << ": building partial paths" << endl;
     const int iterations_to_use_depth = 100000;
     Bool depth_on = False;
     restart_partial:
     vec< vec<int> > vpaths;
     vpaths_final.clear( );
     for ( int x = 0; x < G.N( ); x++ )
     {    if ( G.Source(x) ) 
          {    vec<int> v;
               v.push_back(x);
               vpaths.push_back(v);    }    }
     int64_t iterations = 0;
     while( vpaths.nonempty( ) )
     {    if ( SUPER_VERBOSITY >= 1 )
          {    if ( iterations % 10000 == 0 )
               {    DPRINT3_TO( hout, iterations, 
                         vpaths.size( ), vpaths_final.size( ) );    }    }
          if ( ++iterations == iterations_to_use_depth )
          {    if ( !depth_on )
               {    depth_on = True;
                    hout << "\nrestarting with depth filter on\n" << endl;
                    goto restart_partial;    }
               else
               {    hout << "too many iterations, giving up" << endl;
                    return;    }    }
          vec<int> p = vpaths.back( );
          if ( SUPER_VERBOSITY >= 2 )
          {    hout << "exploring path";
               for ( int j = 0; j < p.isize( ); j++ )
                    hout << " " << p[j];
               hout << endl;    }
          vpaths.pop_back( );
          int x = p.back( );

          // If a path has no extension, it is declared a final path.
          // If a path has a unique extension, we take it.

          if ( G.From(x).empty( ) ) vpaths_final.push_back(p);
          if ( G.From(x).solo( ) )
          {    int y = G.From(x)[0];
               p.push_back(y);
               vpaths.push_back(p);    }
          if ( G.From(x).size( ) <= 1 ) continue;

          // Now go through each extension of the path p.  The structure 
          // rpos_ppos_all captures what we learn about each extension.

          vec< vec< vec< vec< pair<int,int> > > > > 
               rpos_ppos_all( G.From(x).size( ) );
          for ( int j = 0; j < G.From(x).isize( ); j++ )
          {    int y = G.From(x)[j];
               vec<int> q(p);
               q.push_back(y);
               rpos_ppos_all[j].resize( rights.size( ) );

               // See if q is acceptable.

               if (depth_on)
               {    Bool unacceptable = False;
                    for ( int d = 1; d <= Min( depth, q.isize( ) ); d++ )
                    {    vec<int> z;
                         z.SetToSubOf( q, q.isize( ) - d, d );
                         int bp = BinPosition( dpaths[d], z );
                         if ( bp < 0 || !acceptable[d][bp] ) 
                              unacceptable = True;    }
                    if (unacceptable) continue;    }

               // The path q is the extension of p.  Now we go through all the 
               // reads.

               #pragma omp parallel for
               for ( int r = 0; r < rights.isize( ); r++ )
               {    vec< vec<int> > walks, walks_final;

                    // Find the alignments that could start the read walking.
                    // These define the initial walks.

                    const KmerAlignSet& a = kaligns[r];
                    for ( int l = 0; l < a.NAligns( ); l++ )
                    {    int u = verts[ q[0] ].second;
                         if ( kaligns[r].U(l) != u ) continue;
                         vec<int> w;
                         w.push_back(l);
                         walks.push_back(w);    }

                    // For each walk, iteratively extend along q.

                    while( walks.nonempty( ) )
                    {    vec<int> w = walks.back( );
                         walks.pop_back( );
                         int d = w.size( );
                         if ( d == q.isize( ) )
                         {    walks_final.push_back(w);
                              continue;    }
                         int l1 = w.back( );
                         int u1 = verts[ q[d-1] ].second, u2 = verts[ q[d] ].second;
                         int pos1 = verts[ q[d-1] ].first;
                         int pos2 = verts[ q[d] ].first;

                         // Identify the possible extensions.

                         vec<int> exts;
                         int v1 = q[d-1], v2 = q[d], jj;
                         for ( jj = 0; jj < G.From(v1).isize( ); jj++ )
                              if ( G.From(v1)[jj] == v2 ) break;
                         ForceAssert( jj < G.From(v1).isize( ) );
                         for ( int z = 0; z < support[r][v1][jj].isize( ); z++ )
                         {    if ( support[r][v1][jj][z].first != l1 ) continue;
                              exts.push_back( support[r][v1][jj][z].second );    }

                         // Find the alignments to u2.  These define the possible
                         // extensions.  Don't allow an alignment to occur twice.

                         for ( int m = 0; m < exts.isize( ); m++ )
                         {    int l2 = exts[m];
                              if ( Member( w, l2 ) ) continue;
                              vec<int> z(w);
                              z.push_back(l2);
                              walks.push_back(z);    }    }

                    // Now go through the final walks, and define the matches
                    // corresponding to each.

                    if ( SUPER_VERBOSITY >= 3 )
                    {    
                         #pragma omp critical
                         {    hout << "\nfinal walks for read " << r << ":\n";
                              for ( int i = 0; i < walks_final.isize( ); i++ )
                              {    const vec<int>& w = walks_final[i];
                                   hout << "(" << i << ")";
                                   for ( int j = 0; j < w.isize( ); j++ )
                                        hout << " " << w[j];
                                   hout << "\n";    }    }    }
                    for ( int i = 0; i < walks_final.isize( ); i++ )
                    {    const vec<int>& w = walks_final[i];
                         vec< pair<int,int> > rpos_ppos;
                         for ( int d = 0; d < w.isize( ); d++ )
                         {    int pos = verts[ q[d] ].first;
                              for ( int m = 0; m < a.Count( w[d] ); m++ )
                              {    rpos_ppos.push( a.Rpos( w[d], m ),
                                        pos + a.Upos( w[d], m ) );    }    }
                         UniqueSort(rpos_ppos);
                         rpos_ppos_all[j][r].push_back(rpos_ppos);    }    }    }

          // Now compare branches.

          int nb = G.From(x).size( );
          vec< vec<int> > support( nb, vec<int>( nb, 0 ) );
          for ( int j1 = 0; j1 < nb; j1++ )
          for ( int j2 = j1 + 1; j2 < nb; j2++ )
          {    int y1 = G.From(x)[j1], y2 = G.From(x)[j2];
               int u1 = verts[y1].second, u2 = verts[y2].second;
               int pos1 = verts[y1].first, pos2 = verts[y2].first;
               int mp = Min( pos1 + unibases[u1].isize( ),
                    pos2 + unibases[u2].isize( ) );

               // Make sure that there is a difference to be detected.

               Bool diff = False;
               if ( pos1 <= pos2 )
               {    for ( int p1 = 0; p1 < unibases[u1].isize( ); p1++ )
                    {    int p2 = p1 + pos1 - pos2;
                         if ( p2 < 0 ) continue;
                         if ( p2 >= unibases[u2].isize( ) ) break;
                         if ( unibases[u1][p1] != unibases[u2][p2] ) 
                         {    diff = True;
                              break;    }    }    }
               else
               {    for ( int p2 = 0; p2 < unibases[u2].isize( ); p2++ )
                    {    int p1 = p2 + pos2 - pos1;
                         if ( p1 < 0 ) continue;
                         if ( p1 >= unibases[u1].isize( ) ) break;
                         if ( unibases[u1][p1] != unibases[u2][p2] ) 
                         {    diff = True;
                              break;    }    }    }
               if ( !diff ) continue;

               // Check support.

               int xsupport1 = 0, xsupport2 = 0;
               for ( int r = 0; r < rights.isize( ); r++ )
               {    int count1 = 0, count2 = 0;
                    for ( int k = 0; k < rpos_ppos_all[j1][r].isize( ); k++ )
                    {    int count = 0;
                         for ( int z = 0; z < rpos_ppos_all[j1][r][k].isize( ); z++ )
                         {    int ppos = rpos_ppos_all[j1][r][k][z].second;
                              if ( ppos <= mp ) count++;    }
                         count1 = Max( count1, count );    }
                    for ( int k = 0; k < rpos_ppos_all[j2][r].isize( ); k++ )
                    {    int count = 0;
                         for ( int z = 0; z < rpos_ppos_all[j2][r][k].isize( ); z++ )
                         {    int ppos = rpos_ppos_all[j2][r][k][z].second;
                              if ( ppos <= mp ) count++;    }
                         count2 = Max( count2, count );    }
                    if ( count1 > count2 ) 
                    {    xsupport1++;
                         if ( SUPER_VERBOSITY >= 2 )
                         {    hout << "xread (right) " << r << " supports " << y1
                                   << " over " << y2 << "\n";    }    }
                    if ( count2 > count1 ) 
                    {    xsupport2++;    
                         if ( SUPER_VERBOSITY >= 2 )
                         {    hout << "xread (right) " << r << " supports " << y2
                                   << " over " << y1 << "\n";    }    }    }
               support[j1][j2] = xsupport1;
               support[j2][j1] = xsupport2;
               if ( SUPER_VERBOSITY >= 2 )
                    PRINT4_TO( hout, y1, y2, xsupport1, xsupport2 );
               const int min_support = 5;
               if ( xsupport1 >= min_support 
                    && xsupport2 <= 1 && G.To(y2).solo( ) )
               {    if ( SUPER_VERBOSITY >= 2 )
                         hout << "xsuggest killing " << y2 << "\n";
                    xsuggested_kills.push_back(y2);    }
               if ( xsupport2 >= min_support 
                    && xsupport1 <= 1 && G.To(y1).solo( ) )
               {    if ( SUPER_VERBOSITY >= 2 )
                         hout << "xsuggest killing " << y1 << "\n";
                    xsuggested_kills.push_back(y1);    }    }    

          // Keep branches.  Note minimal gating, should be made more stringent.

          vec<Bool> keep( nb, True );
          if ( nb == 2 )
          {    if ( support[0][1] >= 2 && support[1][0] == 0 ) keep[1] = False;
               if ( support[1][0] >= 2 && support[0][1] == 0 ) keep[0] = False;
               if ( support[0][1] >= 5 && support[1][0] == 1 ) keep[1] = False;
               if ( support[1][0] >= 5 && support[0][1] == 1 ) keep[0] = False;    }
          int keepers = 0;
          for ( int j = 0; j < nb; j++ )
          {    if ( !keep[j] ) continue;
               int y = G.From(x)[j];
               vec<int> q(p);
               q.push_back(y);
               int count = 0;
               for ( int r = 0; r < rights.isize( ); r++ )
                    if ( rpos_ppos_all[j][r].nonempty( ) ) count++;
               if ( count == 0 ) 
               {    if ( SUPER_VERBOSITY >= 2 )
                         hout << "no support at all for " << y << "\n";    }
               else 
               {    keepers++;
                    vpaths.push_back(q);    }    }
          if ( keepers == 0 ) vpaths_final.push_back(p);    }
     UniqueSort(xsuggested_kills);    }

     // Compute matches for all reads.
                    
     hout << Date( ) << ": computing matches for all reads\n";
     vec< vec<int> > unis_all( rights.size( ) );
     vec< vec< vec< pair<int,int> > > > umatches_all( rights.size( ) );
     vec<String> ureport( rights.size( ) );
     const Bool PRINT_UMATCHES = True;
     #pragma omp parallel for
     for ( int xxx = 0; xxx < rights.isize( ); xxx++ )
     {    Bool fw = rights[xxx].second;
          int id = rights[xxx].first;

          // Find matches of read to unibases.

          basevector r = longreads[id];
          if ( !fw ) r.ReverseComplement( );

          KmerAlignSet alx( r, L, Ulocs ), aly;
          ClusterAlignsOld( alx, aly );
          KillInferiorClusters(aly);

          vec<int>& unis = unis_all[xxx];
          for ( int j = 0; j < verts.isize( ); j++ )
               unis.push_back( verts[j].second );
          UniqueSort(unis);
          vec< vec< pair<int,int> > >& umatches = umatches_all[xxx];
          umatches.resize( unis.size( ) );
          ostringstream rout;
          for ( int uj = 0; uj < unis.isize( ); uj++ )
          {    int v = unis[uj];
               const basevector& t = unibases[v];
               for ( int i = 0; i < aly.NAligns( ); i++ )
               {    if ( aly.U(i) == v )
                    {    for ( int j = 0; j < aly.Count(i); j++ )
                              umatches[uj].push( aly.RposUpos(i,j) );    }    }
               Sort( umatches[uj] );
               if (PRINT_UMATCHES)
               {    rout << "\nmatches of read(right) " << id << "(" << xxx << ")"
                         << " to unibase " << v << "\n";
                    for ( int j = 0; j < umatches[uj].isize( ); j++ )
                    {    if ( j > 0 && j % 6 == 0 ) rout << "\n";
                         rout << "(" << umatches[uj][j].first << ","
                              << umatches[uj][j].second << ")";    }
                    rout << "\n";    }    }    
          ureport[xxx] = rout.str( );    }
     if (PRINT_UMATCHES)
     {    for ( int j = 0; j < rights.isize( ); j++ )
               hout << ureport[j];    }

     // Find all paths, then discard ones yielding duplicate sequence.

     hout << Date( ) << ": finding all paths" << endl;
     vec< vec<int> > paths;
     if (NEW_FILTER) paths = vpaths_final;
     else G.AllPaths( -1, -1, paths );
     hout << Date( ) << ": constructing bpaths" << endl;
     vec<basevector> bpaths( paths.size( ) );
     #pragma omp parallel for
     for ( int i = 0; i < paths.isize( ); i++ )
     {    const vec<int>& p = paths[i];
          basevector& t = bpaths[i];
          t = unibases[ verts[ p[0] ].second ];
          for ( int l = 1; l < p.isize( ); l++ )
          {    int o = verts[ p[l-1] ].first
                    + unibases[ verts[ p[l-1] ].second ].isize( )
                    - verts[ p[l] ].first;
               t.resize( t.isize( ) - o );
               t = Cat( t, unibases[ verts[ p[l] ].second ] );    }     }
     ParallelUniqueSortSync( bpaths, paths );

     // Define partial paths.  Way dumb.

     hout << Date( ) << ": constructing ppaths" << endl;
     vec< vec<int> > ppaths;
     for ( int i = 0; i < paths.isize( ); i++ )
     {    const vec<int>& p = paths[i];
          for ( int j = 1; j <= p.isize( ); j++ )
          {    vec<int> q(p);
               q.resize(j);
               ppaths.push_back(q);    }    }
     ParallelUniqueSort(ppaths);
     hout << Date( ) << ": constructing bppaths" << endl;
     vec<basevector> bppaths( ppaths.size( ) );
     #pragma omp parallel for
     for ( int i = 0; i < ppaths.isize( ); i++ )
     {    const vec<int>& p = ppaths[i];
          basevector& t = bppaths[i];
          t = unibases[ verts[ p[0] ].second ];
          for ( int l = 1; l < p.isize( ); l++ )
          {    int o = verts[ p[l-1] ].first
                    + unibases[ verts[ p[l-1] ].second ].isize( )
                    - verts[ p[l] ].first;
               t.resize( t.isize( ) - o );
               t = Cat( t, unibases[ verts[ p[l] ].second ] );    }     }

     // Initial proof of principle - unwind reads onto graph and Smith-Waterman
     // against all paths.

     PRINT2_TO( hout, ppaths.size( ), rights.size( ) );
     vec<double> votes( ppaths.size( ), 0 );
     hout << Date( ) << ": going through the reads" << endl;
     #pragma omp parallel for
     for ( int xxx = 0; xxx < rights.isize( ); xxx++ )
     {    ostringstream rout;
          rout << "\n";
          Bool fw = rights[xxx].second;
          int id = rights[xxx].first;
          vec<Bool> ppath_seen( ppaths.size( ), False );

          // Find matches of read to unibases.

          basevector r = longreads[id];
          if ( !fw ) r.ReverseComplement( );
          const vec<int>& unis = unis_all[xxx];
          const vec< vec< pair<int,int> > >& umatches = umatches_all[xxx];

          // Set up to track alignments of the read to a candidate correction.
          // This is for efficiency.  The map A takes as input an interval on
          // the read and the basevector it is to be aligned to.  As output
          // it returns the alignment and the error count.

          map< pair<ho_interval,basevector>, pair<align,double> > A;

          // Proceed.

          const int infinity = 1000000000;
          vec<double> errs( paths.size( ), 0.0 );
          vec<int> Lmatches( paths.size( ), 0 );
          vec< triple<int,int,String> > ereports( paths.size( ) );
          vec<int> lastpos( paths.size( ) ), last( paths.size( ) );
          vec<Bool> skipped( paths.size( ), False );
          int L = heur.L;
          for ( int i = 0; i < paths.isize( ); i++ )
          {    AlignReadToPath( r, id, fw, paths[i], i, unibases, verts,  unis,
                    umatches, ppaths, ppath_seen, L, PRINT_MATCHES, DUMP_LOCAL, 
                    PRINT_LM1, errs[i], Lmatches[i], ereports[i], lastpos[i], 
                    last[i], skipped[i], A );    }

          vec<Bool> discard( ereports.size( ), False );

          // VERSION 3
          for ( int i1 = 0; i1 < ereports.isize( ); i1++ )
          {    for ( int i2 = 0; i2 < ereports.isize( ); i2++ )
               {    if ( skipped[i1] || skipped[i2] ) continue;
                    if ( 10*Lmatches[i1] - errs[i1] > 10*Lmatches[i2] - errs[i2] )
                         discard[i2] = True;    }    }

          int keepers = discard.isize( ) - Sum(discard) - Sum(skipped);
          for ( int i = 0; i < ereports.isize( ); i++ )
          {    if ( skipped[i] ) continue;
               if ( discard[i] ) 
               {    if (PRINT_DISCARDS) rout << "DISCARD: " << ereports[i].third;
                    continue;    }
               if (PRINT_DISCARDS) rout << "KEEP: ";
               rout << ereports[i].third;
               if ( last[i] < 0 ) continue;
               vec<int> p = paths[i];
               p.resize( last[i] + 1 );
               int pp = BinPosition( ppaths, p );
               ForceAssertGe( pp, 0 ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
               #pragma omp critical
               {    votes[pp] += 1.0 / double(keepers);    }    }    
          #pragma omp critical
          {    hout << rout.str( );    }    }

     // Summarize votes.

     hout << Date( ) << ": summarizing votes" << endl;
     if ( votes.empty( ) ) 
     {    hout << "no votes" << endl;
          hout << Date( ) << ": done with unibase " << u 
               << ", time used = " << TimeSince(clock) << endl;
          return;    }
     hout << "\nvotes for unibase " << u << ":\n\n";
     double maxvotes = Max(votes);
     const double vote_thresh = 4.0;
     Bool bad = False;
     for ( int i = 0; i < ppaths.isize( ); i++ )
     {    if ( votes[i] > maxvotes/vote_thresh )
          {    hout << "{" << i << "} ";
               const vec<int>& p = ppaths[i];
               for ( int l = 0; l < p.isize( ); l++ )
               {    if ( l > 0 ) hout << " -->";
                    hout << " [" << p[l] << "]";    }
               hout << "; votes = " << votes[i];
               if (VALIDATE1)
               {    vec<placementx> places = FindGenomicPlacements( 
                         bppaths[i], LG, genome2, Glocs );
                    if ( places.nonempty( ) ) hout << " [TRUE]";
                    else if ( votes[i] == maxvotes ) bad = True;    }
               hout << endl;    }    }
     if (bad)
     {    hout << "\nWarning: best extension of unibase " << u
               << " is nongenomic." << endl;    }

     // Vote on tree.

     hout << "\nVote tree analysis for unibase " << u << ":\n";
     int np = 0;
     for ( int j = 0; j < ppaths.isize( ); j++ )
          np = Max( np, ppaths[j].isize( ) );
     vec< vec<int> > answer;
     vec<double> answer_votes;
     FollowVotes( 0, 0, ppaths.size( ), np, ppaths, votes, verts, unibases, K, 
          heur, hout, answer, answer_votes );
     hout << "\nFINAL ANSWER FOR UNIBASE " << u << "\n\n";
     vec<Bool> genomic( answer.size( ), False );
     vec<int> ext( answer.size( ) );
     vec<basevector> banswer;
     vec< vec<int> > ianswer;
     for ( int j = 0; j < answer.isize( ); j++ )
     {    const vec<int>& p = answer[j];
          for ( int l = 0; l < p.isize( ); l++ )
          {    if ( BinMember( xsuggested_kills, p[l] ) )
                    hout << "Warning: would have killed " << p[l] << "\n";    }
          vec<int> us;
          for ( int l = 0; l < p.isize( ); l++ )
               us.push_back( verts[ p[l] ].second );
          hout << "{" << j << "} ";
          for ( int l = 0; l < p.isize( ); l++ )
          {    if ( l > 0 ) hout << " --> ";
               hout << "[" << p[l] << "]";    }
          ext[j] = verts[ p.back( ) ].first
               + unibases[ verts[ p.back( ) ].second ].isize( );
          hout << "; votes = " << answer_votes[j] << ", extension = " << ext[j];
          int bpos = BinPosition( ppaths, p );
          ForceAssertGe( bpos, 0 );
          banswer.push_back( bppaths[bpos] );
          if (VALIDATE1)
          {    vec<placementx> places = FindGenomicPlacements( 
                    bppaths[bpos], LG, genome2, Glocs );
               if ( places.nonempty( ) ) 
               {    hout << " [TRUE]";

               Bool unacceptable = False;
               for ( int d = 1; d <= Min( depth, p.isize( ) ); d++ )
               {    vec<int> z;
                    z.SetToSubOf( p, p.isize( ) - d, d );
                    int bp = BinPosition( dpaths[d], z );
                    if ( bp < 0 || !acceptable[d][bp] ) unacceptable = True;    }
               if (unacceptable)
               {    hout << " [UNACCEPTABLE]";    }

                    genomic[j] = True;    }    }
          hout << "\n=";
          for ( int l = 0; l < us.isize( ); l++ )
               hout << " " << us[l];
          hout << endl;
          ianswer.push_back(us);    }
     if (VALIDATE1)
     {    Bool good = answer.solo( ) && genomic[0] && ext[0] >= 1000;
          if ( !good) hout << "\nFUNNY!\n";    }
     if (PRINT_BASES2)
     {    for ( int j = 0; j < banswer.isize( ); j++ )
               banswer[j].Print( hout, j );    }
     #pragma omp critical
     {    for ( int j = 0; j < banswer.isize( ); j++ )
          {    new_stuff.push_back_reserve( banswer[j] );    
               right_exts.push_back( ianswer[j] );    }    }
     hout << Date( ) << ": done with unibase " << u 
          << ", time used = " << TimeSince(clock) << endl;    }

void Validate( vecbasevector all, const vecbasevector& unibases,
     const vecbasevector& genome, const int min_mult, 
     const Bool PRINT_MISSING_KMERS )
{    
     const int M = 400;
     vec< pair< kmer<M>, int > > kmers;
     int owned = 0, total = 0;
     all.reserve( all.size( ) + unibases.size( ) );
     for ( size_t j = 0; j < unibases.size( ); j++ )
          all.push_back( unibases[j] );
     vec<int64_t> starts;
     starts.push_back(0);
     for ( size_t i = 0; i < all.size( ); i++ )
     {    const basevector& u = all[i];
          starts.push_back( starts.back( ) + Max( 0, u.isize( ) - M + 1 ) );    }
     kmers.resize( starts.back( ) );
     #pragma omp parallel for
     for ( size_t i = 0; i < all.size( ); i++ )
     {    const basevector& u = all[i];
          kmer<M> x, xrc;
          for ( int j = 0; j <= u.isize( ) - M; j++ )
          {    int64_t r = starts[i] + j;
               x.SetToSubOf( u, j );
               xrc = x;
               xrc.ReverseComplement( );
               kmers[r] = make_pair( ( x <= xrc ? x : xrc ), 0 );    }    }
     starts.clear( );
     starts.push_back(0);
     for ( size_t i = 0; i < genome.size( ); i++ )
     {    const basevector& u = genome[i];
          starts.push_back( starts.back( ) + Max( 0, u.isize( ) - M + 1 ) );    }
     int64_t nkmers = kmers.size( );
     kmers.resize( nkmers + starts.back( ) );
     for ( size_t i = 0; i < genome.size( ); i++ )
     {    const basevector& u = genome[i];
          #pragma omp parallel for
          for ( int jz = 0; jz <= u.isize( ) - M; jz += 10000 )
          {    kmer<M> x, xrc;
               for ( int j = jz; j <= Min( u.isize( ) - M, jz + 10000 ); j++ )
               {    int64_t r = nkmers + starts[i] + j;
                    x.SetToSubOf( u, j );
                    xrc = x;
                    xrc.ReverseComplement( );
                    kmers[r] = make_pair( ( x <= xrc ? x : xrc ), 1 );    }    }    }
     ParallelSort(kmers);
     int missing = 0, bads = 0;
     for ( int64_t i = 0; i < (int64_t) kmers.size( ); i++ )
     {    int64_t j;
          int mult0 = 0;
          for ( j = i; j < (int64_t) kmers.size( ); j++ )
          {    if ( kmers[j].first != kmers[i].first ) break;
               if ( kmers[j].second == 0 ) mult0++;    }
          if ( kmers[j-1].second == 1 )
          {    if ( mult0 >= min_mult ) owned++;
               else
               {    if ( PRINT_MISSING_KMERS && missing++ % 20 == 0 )
                    {    basevector b(M);
                         kmers[i].first.GetBasevector(b);
                         b.Print( cout, "missing_" + ToString(missing) );    }    }
               total++;    }
          else if ( mult0 >= min_mult ) bads++;
          i = j - 1;    }

     // Report results of validation.

     cout << "missing = " << ToStringAddCommas( total - owned ) << " = "
          << PERCENT_RATIO( 5, total - owned, total ) << endl;
     cout << "bads = " << ToStringAddCommas(bads) << " = " 
          << PERCENT_RATIO( 3, bads, total ) << endl;    }

template<int K> 
void BuildPatches(
        const vecbasevector& unibases, 
        const vec< vec<int> >& nexts,
        const vec<Bool>& cn_plus,
        const int L, 
        const vec< vec< pair<int,int> > >& Ulocs, 
        vec<GapPatcher> & patchers_all, 
        const Bool CORRECT_PATCHES, 
        const Bool CORRECT_PATCHES_NEW, 
        const Bool CORRECT_PATCH_VERBOSE,
        const Bool ANNOUNCE_PATCH_CORRECTION,
        const vecbasevector& genome2, 
        const int PATCH_VERBOSITY, 
        const Bool VALIDATE_PATCHES, 
        const int PATCH_CORRECT_VERBOSITY, 
        const int LG, 
        const vec< vec< pair<int,int> > > & Glocs, 
        const String & data_dir, 
        const String & run_dir, 
        vec<basevector> & bpatches, 
        const uint n_patchers_min, 
        const vecbasevector& fbases,
        const vecqualvector& fquals,
        const PairsManager& fpairs,
        const vec< kmer<20> >& fheads,
        const vec<int64_t>& fids
                    )
{
  const unsigned nb_rad_SW = 10;  // the radius of the local Smith-Waterman
  const unsigned sz_patcher_min = 2 * L;
  const unsigned sz_padding_min = 2 * nb_rad_SW + L;
  const int patch_mode = 2; // 0:shortest size  1:median size  2:median gap

  // Hash unibases.
  
  vec<int64_t> starts;
  starts.push_back(0);
  for (size_t i = 0; i < unibases.size(); i++) {
    const basevector& u = unibases[i];
    starts.push_back(starts.back() + u.isize() - K + 1);
  }
  vec< triple<kmer<K>,int,int> > kmers_plus(starts.back());
  #pragma omp parallel for schedule(dynamic, 1)
  for (size_t i = 0; i < unibases.size(); i++) {
    const basevector& u = unibases[i];
    for (int j = 0; j <= u.isize() - K; j++) {
      int64_t r = starts[i] + j;
      kmers_plus[r].first.SetToSubOf(u, j);
      kmers_plus[r].second = i, kmers_plus[r].third = j;
    }
  }
  ParallelSort(kmers_plus);

  // Form consensus for patchers.  This code is copied from LongReadPostPatcher.

  Sort(patchers_all);
  DPRINT(patchers_all.size());

  // ---- Subdivide patchers by gaps and filter out poor patchers

  vec< vec<GapPatcher> > patchers_gap;
  patchers_gap_collect(unibases, patchers_all, n_patchers_min, sz_padding_min,
                       & patchers_gap);
  const size_t n_gaps = patchers_gap.size();
  cout << Date() << ": finding patches for " << n_gaps << " gaps" << endl;

  // ---- Choose, for each gap, the first-guess patcher 
     
  vec<unsigned> ips_best;
  patchers_first_guesses(unibases, patchers_gap, patch_mode, sz_patcher_min, 
                         & ips_best);

  // ---- Sort gaps according to their predicted performance
  //      with omp use schedule(dynamic, 1) to run largest jobs first 

  vec<unsigned> is_gaps;
  indexes_gaps_for_parallel(patchers_gap, ips_best, & is_gaps);

  // ---- Declare some stuff and optimize the guess patchers

  cout << "\n" << Date() << ": patching" << endl;

  vec<String> yreports(n_gaps);
  vec<bool> patched(unibases.size(), false);
  vec<GapPatcher0> all_paths(unibases.size());


  #pragma omp parallel for schedule(dynamic, 1)
  for (unsigned j_gap = 0; j_gap < n_gaps; j_gap++) {

    const unsigned i_gap = is_gaps[j_gap];
    const vec<GapPatcher> & patchers = patchers_gap[i_gap];
    const unsigned n_patchers = patchers.size();

    const size_t    ibv1 = patchers[0].t1;
    const size_t    ibv2 = patchers[0].t2;

    if ( cn_plus.nonempty( ) && ( cn_plus[ibv1] || cn_plus[ibv2] ) ) continue;

    const unsigned ip_best = ips_best[i_gap];
    if (ip_best >= n_patchers) {  // invalid best
      if (PATCH_VERBOSITY >= 1) 
        #pragma omp critical
        cout << Date() << ": i_gap= " << i_gap << ": invalid best patch." << endl; 
    }
    else {

      GapPatcher0 & p_opt = all_paths[ibv1];
      vec<double> timers(3, 0);
      patcher_optimal(unibases, patchers, ip_best, L, nb_rad_SW, sz_padding_min, i_gap, 
                      & p_opt, PATCH_VERBOSITY, & timers);
      patched[ibv1] = true;





      // Correct and save the patch.
      const BaseVec & bv1  = unibases[ibv1];
      const BaseVec & bv2  = unibases[ibv2];
         
      int start1 = p_opt.upos1;
      int stop2 = p_opt.upos2;
      vec<basevector> reps;
      reps.push_back(p_opt.r);
      assembly_edit edit(assembly_edit::GAP_CLOSER, ibv1, start1, ibv2, stop2, reps);
      Bool DIRECT = False;
      basevector bv_corr;
      ostringstream out_tmp;
        if (ANNOUNCE_PATCH_CORRECTION)
        {
             #pragma omp critical
             cout << Date( ) << ": begin correction for gap " << i_gap << endl;    }
        Bool patch_corrected = 
             CorrectPatch3( unibases[ibv1], unibases[ibv2], fbases, fquals, fpairs, 
             fheads, fids, edit, out_tmp, CORRECT_PATCH_VERBOSE, 
             genome2, LG, Glocs );
        if (ANNOUNCE_PATCH_CORRECTION)
        {
             #pragma omp critical
             cout << Date( ) << ": end correction for gap " << i_gap << endl;    }
        if ( !patch_corrected ) 
        {    if ( PATCH_CORRECT_VERBOSITY >= 1 )
                  out_tmp << "patch not corrected, ignoring patch\n";
             yreports[i_gap] = out_tmp.str();
             continue;    }
        bv_corr = Cat(basevector(unibases[ibv1], 0, edit.Start1()), 
                      edit.Rep(0),
                      basevector(unibases[ibv2], edit.Stop2(),
                                 unibases[ibv2].isize() - edit.Stop2())); 

        // Determine if patch is buildable from existing unibases.

        Bool buildable = False;
        if (VALIDATE_PATCHES)
        {    vec< pair<int,int> > upos;
             upos.push( ibv1, 0 );
             while( upos.nonempty( ) )
             {    int u = upos.back( ).first, pos = upos.back( ).second;
                  upos.pop_back( );
                  Bool match = True;
                  for ( int j = 0; j < unibases[u].isize( ); j++ )
                  {    if ( pos >= bv_corr.isize( ) ) break;
                       if ( unibases[u][j] != bv_corr[pos] )
                       {    match = False;
                            break;    }
                       pos++;    }
                  if ( !match ) break;
                  if ( pos == bv_corr.isize( ) )
                  {    buildable = True;
                       break;    }    
                  for ( int j = 0; j < nexts[u].isize( ); j++ )
                       upos.push( nexts[u][j], pos - (K-1) );    } 
              if (buildable) out_tmp << "buildable!\n";    }

#pragma omp critical
      {
        bpatches.push_back(bv_corr);
      }
         
      // Validate patch.
         
      if ( VALIDATE_PATCHES && !buildable ) {
        vec<placementx> totals = FindGenomicPlacements(bv_corr, LG, genome2, Glocs);
#pragma omp critical
        {
          out_tmp << "matches for patch to gap " << i_gap 
               << " from " << ibv1 << " to " << ibv2 << "\n";
          for (int l = 0; l < totals.isize(); l++) {
            out_tmp << "perfectly matches " << totals[l].g << "."
                    << totals[l].pos << "-" << totals[l].pos + bv_corr.isize() 
                    << "\n";
          }
          if (totals.empty()) {
            {    Ofstream(out, "slobber.fasta");
                 bv_corr.Print(out, "patch_for_gap_" + ToString(j_gap));    }
            out_tmp << AllOfOutput1("QueryLookupTable K=12 MM=12 MC=0.15 "
                                     "SEQS=slobber.fasta L=" + data_dir + "/genome.lookup "
                                     "VISUAL=True NH=True QUIET=True");
          }
        }
      }
         
      yreports[i_gap] = out_tmp.str();

    }
  }

  cout << endl;

  int yc = 1;
  for (int i = 0; i < yreports.isize(); i++) {
    if (yreports[i].size() > 0) {
      cout << "\n[" << yc++ << "] " << yreports[i];
    }
  }
  cout << Date() << ": " << bpatches.size() 
       << " gap patches prepared" << endl;
}

template void BuildPatches<96>(const vecbasevector& unibases, 
                               const vec< vec<int> >& nexts,
                               const vec<Bool>& cn_plus,
                               const int L, 
                               const vec< vec< pair<int,int> > >& Ulocs, 
                               vec<GapPatcher>& patchers, 
                               const Bool CORRECT_PATCHES, 
                               const Bool CORRECT_PATCHES_NEW, 
                               const Bool CORRECT_PATCH_VERBOSE,
                               const Bool ANNOUNCE_PATCH_CORRECTION,
                               const vecbasevector& genome2, 
                               const int PATCH_VERBOSITY, 
                               const Bool VALIDATE_PATCHES, 
                               const int PATCH_CORRECT_VERBOSITY, 
                               const int LG, 
                               const vec< vec< pair<int,int> > >& Glocs, 
                               const String& data_dir, 
                               const String& run_dir, 
                               vec<basevector>& bpatches, 
                               const uint n_patchers_min,
                               const vecbasevector& fbases,
                               const vecqualvector& fquals,
                               const PairsManager& fpairs,
                               const vec< kmer<20> >& fheads,
                               const vec<int64_t>& fids );

template void BuildPatches<640>(const vecbasevector& unibases, 
                                const vec< vec<int> >& nexts,
                                const vec<Bool>& cn_plus,
                                const int L, 
                                const vec< vec< pair<int,int> > >& Ulocs, 
                                vec<GapPatcher>& patchers, 
                                const Bool CORRECT_PATCHES, 
                                const Bool CORRECT_PATCHES_NEW, 
                                const Bool CORRECT_PATCH_VERBOSE,
                               const Bool ANNOUNCE_PATCH_CORRECTION,
                                const vecbasevector& genome2, 
                                const int PATCH_VERBOSITY, 
                                const Bool VALIDATE_PATCHES, 
                                const int PATCH_CORRECT_VERBOSITY, 
                                const int LG, 
                                const vec< vec< pair<int,int> > >& Glocs, 
                                const String& data_dir, 
                                const String& run_dir, 
                                vec<basevector>& bpatches, 
                                const uint n_patchers_min,
                                const vecbasevector& fbases,
                                const vecqualvector& fquals,
                                const PairsManager& fpairs,
                                const vec< kmer<20> >& fheads,
                                const vec<int64_t>& fids );

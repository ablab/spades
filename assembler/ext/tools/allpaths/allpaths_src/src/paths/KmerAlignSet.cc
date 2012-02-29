///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "CoreTools.h"
#include "Equiv.h"
#include "STLExtensions.h"
#include "math/Functions.h"
#include "math/HoInterval.h"
#include "paths/KmerAlignSet.h"
#include "paths/LongReadTools.h"

void ClusterAlignsNew( const KmerAlignSet& in, KmerAlignSet& out, const Bool clean,
     const int min_spread )
{    for ( int i = 0; i < in.NAligns( ); i++ )
     {    int u = in.U(i);
          vec< pair<int,int> > rpos_upos;
          int n = in.Count(i);
          for ( int j = 0; j < n; j++ )
               rpos_upos.push( in.Rpos(i,j), in.Upos(i,j) );
          Sort(rpos_upos);
          vec<int> from( n, -1 ), to( n, -1 );
          vec< triple<int,int,int> > dist;
          const int udist_mult = 5;
          const int max_dist = 1000;
          const int delta_r_max = 350;
          for ( int l1 = 0; l1 < n; l1++ )
          {    int rpos1 = rpos_upos[l1].first, upos1 = rpos_upos[l1].second;
               for ( int l2 = l1 + 1; l2 < n; l2++ )
               {    int rpos2 = rpos_upos[l2].first, upos2 = rpos_upos[l2].second;
                    int delta_r = rpos2 - rpos1, delta_u = upos2 - upos1;
                    if ( delta_r > delta_r_max ) break;
                    if ( delta_r == 0 || delta_u <= 0 ) continue;
                    int offset1 = upos1 - rpos1, offset2 = upos2 - rpos2;
                    int delta_offset = offset2 - offset1;
                    double d = delta_r + udist_mult * Abs(delta_offset);
                    if ( d <= max_dist ) dist.push( d, l1, l2 );    }    }
          Sort(dist);
          for ( int k = 0; k < dist.isize( ); k++ )
          {    int l1 = dist[k].second, l2 = dist[k].third;
               if ( from[l1] >= 0 || to[l2] >= 0 ) continue;
               from[l1] = l2, to[l2] = l1;    }
          equiv_rel e(n);
          for ( int l = 0; l < n; l++ )
               if ( from[l] >= 0 ) e.Join( l, from[l] );
          vec<int> reps;
          e.OrbitRepsAlt(reps);
          for ( int l = 0; l < reps.isize( ); l++ )
          {    vec< pair<int,int> > a;
               vec<int> o;
               e.Orbit( reps[l], o );
               Sort(o);
               if ( rpos_upos[ o.back( ) ].first
                    - rpos_upos[ o.front( ) ].first >= min_spread )
               {    for ( int x = 0; x < o.isize( ); x++ )
                         a.push_back( rpos_upos[ o[x] ] );
                    if (clean)
                    {    const int min_prox = 50;
                         int start, stop;
                         for ( start = 0; start < a.isize( ); start++ )
                         {    int p = BinPosition( rpos_upos, a[start] );
                              Bool too_close = False;
                              if ( p > 0 
                                   && rpos_upos[p].first == rpos_upos[p-1].first
                                   && rpos_upos[p].second - rpos_upos[p-1].second
                                        < min_prox )
                              {    too_close = True;    }
                              if ( p < rpos_upos.isize( ) - 1
                                   && rpos_upos[p+1].first == rpos_upos[p].first
                                   && rpos_upos[p+1].second - rpos_upos[p].second
                                        < min_prox )
                              {    too_close = True;    }
                              if ( !too_close ) break;    }
                         for ( stop = a.isize( ) - 1; stop >= 0; stop-- )
                         {    int p = BinPosition( rpos_upos, a[stop] );
                              Bool too_close = False;
                              if ( p > 0 
                                   && rpos_upos[p].first == rpos_upos[p-1].first
                                   && rpos_upos[p].second - rpos_upos[p-1].second
                                        < min_prox )
                              {    too_close = True;    }
                              if ( p < rpos_upos.isize( ) - 1
                                   && rpos_upos[p+1].first == rpos_upos[p].first
                                   && rpos_upos[p+1].second - rpos_upos[p].second
                                        < min_prox )
                              {    too_close = True;    }
                              if ( !too_close ) break;    }
                         stop++;
                         if ( stop <= start ) a.clear( );
                         else a.SetToSubOf( a, start, stop - start );    }
                    if ( a.nonempty( ) ) out.AddAlign( u, a );    }    }    }    }

// Form an equivalence relation on the matches between the read and a
// given unibase.  This has not been thought through particularly carefully.

void ClusterAlignsOld( const KmerAlignSet& in, KmerAlignSet& out )
{    for ( int i = 0; i < in.NAligns( ); i++ )
     {    int u = in.U(i);
          int n = in.Count(i);
          const int delta_r_max = 350;
          const int delta_sub = 10;
          const double max_delta_ratio = 1.5;
          const int min_spread = 10;
          equiv_rel e(n);
          vec< pair<int,int> > rpos_upos;
          for ( int j = 0; j < n; j++ )
               rpos_upos.push( in.Rpos(i,j), in.Upos(i,j) );
          Sort(rpos_upos);
          for ( int l1 = 0; l1 < n; l1++ )
          {    int rpos1 = rpos_upos[l1].first, upos1 = rpos_upos[l1].second;
               for ( int l2 = l1 + 1; l2 < n; l2++ )
               {    int rpos2 = rpos_upos[l2].first;
                    int upos2 = rpos_upos[l2].second;
                    int delta_r = rpos2 - rpos1, delta_u = upos2 - upos1;
                    if ( delta_r == 0 ) continue;
                    if ( delta_r > delta_r_max ) break;
                    if ( delta_u <= 0 ) continue;
                    if ( double( delta_r - delta_sub ) / double(delta_u)
                         > max_delta_ratio )
                    {    continue;    }
                    if ( double( delta_u - delta_sub ) / double(delta_r)
                         > max_delta_ratio )
                    {    continue;    }
                    e.Join( l1, l2 );    }    }
          vec<int> reps;
          e.OrbitRepsAlt(reps);
          for ( int l = 0; l < reps.isize( ); l++ )
          {    vec< pair<int,int> > a;
               vec<int> o;
               e.Orbit( reps[l], o );
               Sort(o);
               if ( rpos_upos[ o.back( ) ].first
                    - rpos_upos[ o.front( ) ].first >= min_spread )
               {    for ( int x = 0; x < o.isize( ); x++ )
                         a.push_back( rpos_upos[ o[x] ] );
                    out.AddAlign( u, a );    }    }    }    }

void KillInferiorClusters( KmerAlignSet& x )
{    vec<Bool> to_remove( x.NAligns( ), False );
     const double min_ratio_to_kill = 5.0;
     for ( int i1 = 0; i1 < x.NAligns( ); i1++ )
          for ( int i2 = 0; i2 < x.NAligns( ); i2++ )
     {    if ( x.X( )[i1].second.front( ).first > x.X( )[i2].second.front( ).first )
               continue;    
          if ( x.X( )[i1].second.back( ).first < x.X( )[i2].second.back( ).first )
               continue;
          if ( double( x.X( )[i1].second.size( ) )
               < min_ratio_to_kill * double( x.X( )[i2].second.size( ) ) )
          {    continue;    }
          /*
          cout << "killing " << x.X( )[i2].first
               << "." << x.X( )[i2].second.front( ).first
               << "-" << x.X( )[i2].second.back( ).first
               << "[" << x.X( )[i2].second.size( ) << "] using "
               << x.X( )[i1].first
               << "." << x.X( )[i1].second.front( ).first
               << "-" << x.X( )[i1].second.back( ).first
               << "[" << x.X( )[i1].second.size( ) << "]\n";
          */
          to_remove[i2] = True;    }
     EraseIf( x.XMutable( ), to_remove );    }

void KillInferiorClustersNew( KmerAlignSet& x, const vecbasevector& unibases,
     const double min_ratio_to_kill )
{    int N = x.NAligns( );
     vec<Bool> to_remove( N, False );
     vec<ho_interval> read_range(N);
     int top = 0;
     for ( int i = 0; i < N; i++ )
     {    read_range[i] = ho_interval( x.X( )[i].second.front( ).first,
               x.X( )[i].second.back( ).first );
          top = Max( top, read_range[i].Stop( ) );    }

     // For efficiency, create an index so that we can find all the alignments
     // overlapping a given window.

     const int window = 10;
     vec< vec<int> > own( top/window + 1 );
     for ( int i = 0; i < N; i++ )
     {    for ( int j = 0; j < own.isize( ); j++ )
          {    if ( Meets( ho_interval( j*window, (j+1)*window ), read_range[i] ) )
                    own[j].push_back(i);    }    }

     // Go through the alignments.

     for ( int i2 = 0; i2 < N; i2++ )
     {    int start = read_range[i2].Start( );
          int j = start/window;
          for ( int l = 0; l < own[j].isize( ); l++ )
          {    int i1 = own[j][l];
          
               // We require that the interval on the read encompassed by the first
               // alignment encompasses that of the second.
     
               if ( read_range[i1].Start( ) > read_range[i2].Start( ) ) continue;
               if ( read_range[i1].Stop( ) < read_range[i2].Stop( ) ) continue;

               int low = 0, high = x.X(i1).second.size( );
               while ( low < x.X(i1).second.isize( )
                    && x.X(i1).second[low].first < x.X(i2).second.front( ).first )
               {    low++;    }
               while ( high > 0
                    && x.X(i1).second[high-1].first > x.X(i2).second.back( ).first )
               {    high--;    }
               if ( double( high - low )
                    < min_ratio_to_kill * double( x.X(i2).second.size( ) ) )
               {    continue;    }
               to_remove[i2] = True;
               break;    }    }

     // If two alignments are inconsistent, and one is twice as big
     // as the other, delete the smaller one.

     /*
     const int min_mult = 2;
     for ( int i1 = 0; i1 < x.NAligns( ); i1++ )
     for ( int i2 = 0; i2 < x.NAligns( ); i2++ )
     {    int u1 = x.U(i1), u2 = x.U(i2);
          if ( i1 == i2 || to_remove[i1] ) continue;
          if ( x.Count(i1) < min_mult * x.Count(i2) ) continue;
          vec<int> o;
          for ( int j2 = 0; j2 < x.Count(i2); j2++ )
          {    int j1 = BinPosition1( x.X(i1).second, x.Rpos(i2,j2) );
               if ( j1 >= 0 ) o.push_back( x.Upos(i1,j1) - x.Upos(i2,j2) );    }
          UniqueSort(o);
          if ( o.empty( ) ) continue;
          Bool valid = False;
          const basevector &U1 = unibases[u1], &U2 = unibases[u2];
          for ( int k = 0; k < o.isize( ); k++ )
          {    Bool mismatch = False;
               for ( int upos1 = Max( 0, o[k] ); 
                    upos1 < Min( U1.isize( ), U2.isize( ) + o[k]  ); upos1++ )
               {    int upos2 = upos1 - o[k];
                    if ( U1[upos1] != U2[upos2] )
                    {    mismatch = True;
                         break;    }    }
               if ( !mismatch )
               {    valid = True;
                    break;    }    }
          if ( !valid ) to_remove[i2] = True;    }
     */

     // Do the deletion.

     EraseIf( x.XMutable( ), to_remove );    }

KmerAlignSet::KmerAlignSet( const basevector& x, const int K, 
     const vec< vec< pair<int,int> > >& Ulocs )
{    vec< triple<int,int,int> > aligns;
     for ( int rpos = 0; rpos <= x.isize( ) - K; rpos++ )
     {    int n = KmerId( x, K, rpos );
          for ( int l = 0; l < Ulocs[n].isize( ); l++ )
          {    int u = Ulocs[n][l].first;
               int upos = Ulocs[n][l].second;
               int offset = upos - rpos;
               aligns.push( u, offset, rpos );    }    }
     Sort(aligns);
     for ( int i = 0; i < aligns.isize( ); i++ )
     {    int j, u = aligns[i].first;
          for ( j = i + 1; j < aligns.isize( ); j++ )
               if ( aligns[j].first != u ) break;
          vec< pair<int,int> > rpos_upos;
          for ( int k = i; k < j; k++ )
               rpos_upos.push( aligns[k].third, aligns[k].second + aligns[k].third );
          AddAlign( u, rpos_upos );    
          i = j - 1;    }    }

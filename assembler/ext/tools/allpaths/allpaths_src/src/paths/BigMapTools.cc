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
#include "Equiv.h"
#include "ParallelVecUtilities.h"
#include "graph/Digraph.h"
#include "graph/DigraphTemplate.h"
#include "kmers/KmerRecord.h"
#include "paths/BigMapTools.h"
#include "paths/ProcessGap.h"
#include "paths/UnipathScaffold.h"

void NotOne( const int K, const vecbasevector& unibases, const vec<int>& to_rc,
     const digraphE<linklet>& G, const int min_kmers, const int max_link_ratio,
     const int dev_mult, const Bool RAISE_VERBOSE, vec<Bool>& not_one )
{
     // Identify unipaths that appear not to have copy number one.

     not_one.resize( unibases.size( ), False );
     for ( int u = 0; u < (int) unibases.size( ); u++ )
     {    int nkmers = unibases[u].isize( ) - K + 1;
          if ( nkmers < min_kmers ) continue;
          int max_links = 0;
          for ( int j = 0; j < G.From(u).isize( ); j++ )
          {    const linklet& l = G.EdgeObjectByIndexFrom( u, j );
               max_links = Max( max_links, l.nlinks );    }
          const int min_links_to_test = 20;
          if ( max_links >= min_links_to_test )
          {    vec<int> goods;
               for ( int j = 0; j < G.From(u).isize( ); j++ )
               {    int v = G.From(u)[j];
               const linklet& l = G.EdgeObjectByIndexFrom( u, j );
                    if ( max_links <= max_link_ratio * l.nlinks ) 
                         goods.push_back(j);    }
               for ( int i1 = 0; i1 < goods.isize( ); i1++ )
               for ( int i2 = i1 + 1; i2 < goods.isize( ); i2++ )
               {    int v1 = G.From(u)[ goods[i1] ], v2 = G.From(u)[ goods[i2] ];
                    if ( v2 == to_rc[v1] ) continue;
                    Bool link12 = False;
                    for ( int j = 0; j < G.From(v1).isize( ); j++ )
                         if ( G.From(v1)[j] == v2 ) link12 = True;
                    for ( int j = 0; j < G.To(v1).isize( ); j++ )
                         if ( G.To(v1)[j] == v2 ) link12 = True;
                    if (link12) continue;
                    const linklet& l1 = G.EdgeObjectByIndexFrom( u, goods[i1] );
                    const linklet& l2 = G.EdgeObjectByIndexFrom( u, goods[i2] );
                    int dev1 = dev_mult * l1.dev, dev2 = dev_mult * l2.dev;
                    Bool OK1 = l1.sep - dev1 + unibases[v1].isize( ) - K + 1 
                         < l2.sep + dev2;
                    Bool OK2 = l2.sep - dev2 + unibases[v2].isize( ) - K + 1 
                         < l1.sep + dev1;
                    if ( RAISE_VERBOSE && !OK1 && !OK2 )
                    {    cout << "\n";
                         PRINT(u);
                         PRINT8( v1, v2, l1.sep, l1.dev, l2.sep, l2.dev,
                              unibases[v1].size( ), unibases[v2].size( ) );    }
                    if ( !OK1 && !OK2 ) 
                    {    not_one[u] = True;    
                         not_one[ to_rc[u] ] = True;    }    }    }
          if ( RAISE_VERBOSE && not_one[u] ) 
          {    cout << u << " not of copy number one?\n";    
               cout << to_rc[u] << " not of copy number one?\n";    }    }    }

template<int K2> void ToBig( const vecbasevector& unibases1,
     const vecbasevector& unibases2, vec< vec< pair<int,int> > >& to_big )
{
     to_big.resize( unibases1.size( ) );
     cout << Date( ) << ": getting unibase2 kmers" << endl;
     vec<int64_t> starts;
     starts.push_back(0);
     for ( size_t i = 0; i < unibases2.size( ); i++ )
     {    const basevector& u = unibases2[i];
          starts.push_back( starts.back( ) + u.isize( ) - K2 + 1 );    }
     vec< triple<kmer<K2>,int,int> > kmers_plus( starts.back( ) );
     #pragma omp parallel for schedule( dynamic, 1 )
     for ( size_t i = 0; i < unibases2.size( ); i++ )
     {    const basevector& u = unibases2[i];
          for ( int j = 0; j <= u.isize( ) - K2; j++ )
          {    int64_t r = starts[i] + j;
               kmers_plus[r].first.SetToSubOf( u, j );
               kmers_plus[r].second = i, kmers_plus[r].third = j;    }    }
     ParallelSort(kmers_plus);
     cout << Date( ) << ": copying" << endl;
     vec< kmer<K2> > kmers( kmers_plus.size( ) );
     for ( size_t i = 0; i < kmers.size( ); i++ )
          kmers[i] = kmers_plus[i].first;
     cout << Date( ) << ": mapping from 1 to 2" << endl;
     #pragma omp parallel for
     for ( size_t i = 0; i < unibases1.size( ); i++ )
     {    const basevector& u1 = unibases1[i];
          vec< pair<int,int> > locs;
          kmer<K2> x;
          for ( int pos1 = 0; pos1 <= u1.isize( ) - K2; pos1++ )
          {    x.SetToSubOf( u1, pos1 );
               int64_t p = BinPosition( kmers, x );
               if ( p >= 0 )
               {    int u2 = kmers_plus[p].second, pos2 = kmers_plus[p].third;
                    int offset = pos2 - pos1;
                    locs.push( u2, offset );
                    while( pos1 < u1.isize( ) - K2 )
                    {    pos1++, pos2++;
                         if ( pos2+K2-1 >= unibases2[u2].isize( ) ) break;
                         if ( u1[pos1+K2-1] != unibases2[u2][pos2+K2-1] ) 
                         {    pos1--;
                              break;    }    
                         else locs.push( u2, offset );    }    }    }
          Sort(locs);
          vec< triple<int,int,int> > locsx;
          for ( int j = 0; j < locs.isize( ); j++ )
          {    int k = locs.NextDiff(j);
               locsx.push( k-j, locs[j].first, locs[j].second );
               j = k - 1;    }
          ReverseSort(locsx);
          if ( locsx.nonempty( ) )
          {    // cout << "\nu1 = " << i << " maps to\n";
               // for ( int j = 0; j < locsx.isize( ); j++ )
               // {    cout << "(" << locsx[j].first << ") " << locsx[j].second
               //           << "." << locsx[j].third << "\n";    }    
               
               // Note that we always take the best match.  Ignoring the other
               // matches (as we do) could be dangerous.

               to_big[i].push( locsx[0].second, locsx[0].third );    }    }    }

template void ToBig<640>( const vecbasevector& unibases1,
     const vecbasevector& unibases2, vec< vec< pair<int,int> > >& to_big );

void RemoveWeakLinks( const vecbasevector& unibases, const vec<int>& to_rc, 
     const vec<Bool>& not_one, digraphE<linklet>& G, const int max_link_ratio, 
     const int dev_mult, const Bool DELETE_VERBOSE )
{
     // If from a given CN1 unipath, we link to two unipaths, and they have to 
     // overlap, and one has far more links than the other, exclude the one having 
     // less links.  However, be more lenient in the case where there are links
     // between them.

     vec<int> to_delete;
     const int dev_mult_len = 6;
     for ( int u = 0; u < (int) unibases.size( ); u++ )
     {    if ( not_one[u] != 1 ) continue;
          for ( int i1 = 0; i1 < G.From(u).isize( ); i1++ )
          for ( int i2 = 0; i2 < G.From(u).isize( ); i2++ )
          {    const linklet& l1 = G.EdgeObjectByIndexFrom( u, i1 );
               const linklet& l2 = G.EdgeObjectByIndexFrom( u, i2 );
               if ( !( l1.nlinks >= max_link_ratio * l2.nlinks ) ) continue;
               int dev1 = dev_mult * l1.dev, dev2 = dev_mult * l2.dev;
               int v1 = G.From(u)[i1], v2 = G.From(u)[i2];
               if ( l1.sep - dev1 + unibases[v1].isize( ) < l2.sep + dev2 ) continue;
               if ( l2.sep - dev2 + unibases[v2].isize( ) < l1.sep + dev1 ) continue;
               if ( Member( G.From(v1), v2 ) )
               {    
               int dev1_len = dev_mult_len * l1.dev; 
               int dev2_len = dev_mult_len * l2.dev;
               if ( l1.sep - dev1_len + unibases[v1].isize( ) < l2.sep + dev2_len )
                    continue;
               if ( l2.sep - dev2_len + unibases[v2].isize( ) < l1.sep + dev1_len )
                    continue;
               }
               int ru = to_rc[u], rv2 = to_rc[v2];
               if (DELETE_VERBOSE)
               {    cout << "method A1, deleting link from " << u << " to " 
                         << v2 << "\n";
                    cout << "method A1, deleting link from " << rv2 << " to " 
                         << ru << "\n";    }
               to_delete.push_back( G.EdgeObjectIndexByIndexFrom(u, i2) );
               for ( int j = 0; j < G.From(rv2).isize( ); j++ )
               {    if ( G.From(rv2)[j] == ru )
                    {    to_delete.push_back( G.EdgeObjectIndexByIndexFrom( 
                              rv2, j ) );    }    }    }    }
     UniqueSort(to_delete);
     G.DeleteEdges(to_delete);

     // Look for heterozygous inversions.

     if (DELETE_VERBOSE)
     {    const double min_ratio_for_het = 0.1;
          const int min_count_for_het = 3;
          cout << "\npossible het inversions:\n";
          for ( int u1 = 0; u1 < (int) unibases.size( ); u1++ )
          {    for ( int ja = 0; ja < G.From(u1).isize( ); ja++ )
               for ( int jb = ja+1; jb < G.From(u1).isize( ); jb++ )
               {    int u2a = G.From(u1)[ja], u2b = G.From(u1)[jb];
                    const linklet& la = G.EdgeObjectByIndexFrom( u1, ja );
                    const linklet& lb = G.EdgeObjectByIndexFrom( u1, jb );
                    int na = la.nlinks, nb = lb.nlinks;
                    double ratio = double(na)/double(nb);
                    if ( Min( ratio, 1.0/ratio ) < min_ratio_for_het ) continue;
                    if ( Min( na, nb ) < min_count_for_het ) continue;
                    if ( u2a == to_rc[u2b] )
                         cout << u1 << " --> {" << u2a << "[links=" << na 
                              << ",sep=" << la.sep << "]," << u2b << "[links=" 
                              << nb << ",sep=" << lb.sep << "]}\n";    }    }
          cout << "\n";    }

     // Remove links to unipaths of predicted CN > 1.

     to_delete.clear( );
     for ( int u = 0; u < (int) unibases.size( ); u++ )
     {    for ( int j = 0; j < G.From(u).isize( ); j++ )
          {    int v = G.From(u)[j];
               if ( not_one[u] || not_one[v] )
               {    if (DELETE_VERBOSE)
                    {    cout << "method B, deleting link from " << u
                              << " to " << v << "\n";    }
                    to_delete.push_back( 
                         G.EdgeObjectIndexByIndexFrom( u, j ) );    }    }    }
     UniqueSort(to_delete);
     G.DeleteEdges(to_delete);

     // Remove links having much stronger support the other way.

     to_delete.clear( );
     for ( int u = 0; u < (int) unibases.size( ); u++ )
     {    for ( int j = 0; j < G.From(u).isize( ); j++ )
          {    int v = G.From(u)[j];
               int nlinks1 = G.EdgeObjectByIndexFrom( u, j ).nlinks;
               for ( int k = 0; k < G.From(v).isize( ); k++ )
               {    if ( G.From(v)[k] == u )
                    {    int nlinks2 = G.EdgeObjectByIndexFrom( v, k ).nlinks;
                         if ( nlinks2 >= max_link_ratio * nlinks1 )
                         {    if (DELETE_VERBOSE)
                              {    cout << "method C, deleting links from "
                                        << u << " to " << v << "\n";    }
                              to_delete.push_back(
                                   G.EdgeObjectIndexByIndexFrom( u, j ) );
                              break;    }    }    }    }    }
     UniqueSort(to_delete);
     G.DeleteEdges(to_delete);

     // Remove links that are very poorly supported relative to other links.

     const int total_ignore_ratio = 20;
     to_delete.clear( );
     for ( int u = 0; u < (int) unibases.size( ); u++ )
     {    int nlinks_max = 0;
          for ( int j = 0; j < G.From(u).isize( ); j++ )
          {    int nlinks = G.EdgeObjectByIndexFrom( u, j ).nlinks;
               nlinks_max = Max( nlinks_max, nlinks );    }
          for ( int j = 0; j < G.From(u).isize( ); j++ )
          {    int nlinks = G.EdgeObjectByIndexFrom( u, j ).nlinks;
               if ( nlinks_max >= total_ignore_ratio * nlinks )
               {    if (DELETE_VERBOSE)
                    {    cout << "method D, deleting links from " << u << " to "
                              << G.From(u)[j] << "\n";    }
                    to_delete.push_back(
                         G.EdgeObjectIndexByIndexFrom( u, j ) );    }    }    }
     UniqueSort(to_delete);
     G.DeleteEdges(to_delete);    }

void RemoveWeakLinks2( const int K2, const vecbasevector& unibases2,
     const vec<int>& to_rc2, const vec<double>& raw2, digraphE<linklet>& G2I, 
     const vec< vec< pair<int,int> > >& nexts2x, const int max_link_ratio, 
     const int dev_mult, const Bool DELETE_VERBOSE )
{
     // Let n be the maximum number of links from any unibase to unibase u.
     // Let d be a constant (taken to be 10).  Let x1,...,xn be the unibases that
     // link to u with at least n/d links.  Let y links to u with less than n/d
     // links.  Then if y does not link to any of the xi, delete the link from y
     // to u.  Repeat in the other direction.

     vec<int> to_delete0;
     const double links_div = 10.0;
     for ( int u = 0; u < (int) unibases2.size( ); u++ )
     {    vec<int> nlinks( G2I.To(u).size( ) );
          for ( int j = 0; j < G2I.To(u).isize( ); j++ )
               nlinks[j] = G2I.EdgeObjectByIndexTo( u, j ).nlinks;
          int max_links = ( nlinks.empty( ) ? 0 : Max(nlinks) );
          for ( int j = 0; j < G2I.To(u).isize( ); j++ )
          {    if ( nlinks[j] >= double(max_links)/links_div ) continue;
               int y = G2I.To(u)[j];
               Bool indirect = False;
               for ( int k = 0; k < G2I.To(u).isize( ); k++ )
               {    if ( nlinks[k] < double(max_links)/links_div ) continue;
                    int x = G2I.To(u)[k];
                    if ( Member( G2I.From(y), x ) ) indirect = True;    }
               if (indirect) continue;
               to_delete0.push_back( G2I.EdgeObjectIndexByIndexTo( u, j ) );
               int ru = to_rc2[u], ry = to_rc2[y];
               for ( int l = 0; l < G2I.From(ru).isize( ); l++ )
               {    if ( G2I.From(ru)[l] == ry )
                    {    to_delete0.push_back( G2I.EdgeObjectIndexByIndexFrom( 
                              ru, l ) );    }    }    }    }
     UniqueSort(to_delete0);
     G2I.DeleteEdges(to_delete0);

     // Some other test (I've forgotten what it is).

     int nuni2 = unibases2.size( );
     vec<Bool> not_one2( nuni2, False );
     for ( int u1 = 0; u1 < nuni2; u1++ )
     {    for ( int j = 0; j < G2I.From(u1).isize( ); j++ )
          {    int u2 = G2I.From(u1)[j];
               const linklet& l = G2I.EdgeObjectByIndexFrom( u1, j );
               int nkmers1 = unibases2[u1].isize( ) - K2 + 1;
               int nkmers2 = unibases2[u2].isize( ) - K2 + 1;
               double cn1 = raw2[u1], cn2= raw2[u2];
               double min_cn_mult = 1.5;
               double min_size_mult = 4.0;
               const int min_links = 10;
               if ( l.nlinks >= min_links && cn1 >= min_cn_mult * cn2
                    && double(nkmers1) <= min_size_mult * double(nkmers2) )
               {    
                    /*
                    if (RAISE_VERBOSE)
                    {    cout << "\n";
                         PRINT7( u1, u2, cn1, cn2, nkmers1, nkmers2, l.nlinks );
                         cout << u1 << " not of copy number one\n";    }
                    */
                    not_one2[u1] = True;    }    }    }

     // Screen links as in RemoveWeak:
     // If from a given CN1 unipath, we link to two unipaths, and they have to 
     // overlap, and one has far more links than the other, exclude the one having 
     // less links.  However, be more lenient in the case where there are links
     // between them.

     vec<int> to_delete2;
     const int dev_mult_len = 6;
     for ( int u = 0; u < (int) unibases2.size( ); u++ )
     {    if ( not_one2[u] ) continue;
          for ( int i1 = 0; i1 < G2I.From(u).isize( ); i1++ )
          for ( int i2 = 0; i2 < G2I.From(u).isize( ); i2++ )
          {    const linklet& l1 = G2I.EdgeObjectByIndexFrom( u, i1 );
               const linklet& l2 = G2I.EdgeObjectByIndexFrom( u, i2 );
               if ( !( l1.nlinks >= max_link_ratio * l2.nlinks ) ) continue;
               int dev1 = dev_mult * l1.dev, dev2 = dev_mult * l2.dev;
               int v1 = G2I.From(u)[i1], v2 = G2I.From(u)[i2];
               if ( l1.sep - dev1 + unibases2[v1].isize( ) < l2.sep + dev2 ) 
                    continue;
               if ( l2.sep - dev2 + unibases2[v2].isize( ) < l1.sep + dev1 ) 
                    continue;
               if ( Member( G2I.From(v1), v2 ) )
               {    
               int dev1_len = dev_mult_len * l1.dev; 
               int dev2_len = dev_mult_len * l2.dev;
               if ( l1.sep - dev1_len + unibases2[v1].isize( ) < l2.sep + dev2_len )
                    continue;
               if ( l2.sep - dev2_len + unibases2[v2].isize( ) < l1.sep + dev1_len )
                    continue;
               }
               int ru = to_rc2[u], rv2 = to_rc2[v2];
               /*
               if (DELETE_VERBOSE)
               {    cout << "method A1, deleting link from " << u << " to " 
                         << v2 << "\n";
                    cout << "method A1, deleting link from " << rv2 << " to " 
                         << ru << "\n";    }
               */
               to_delete2.push_back( G2I.EdgeObjectIndexByIndexFrom(u, i2) );
               for ( int j = 0; j < G2I.From(rv2).isize( ); j++ )
               {    if ( G2I.From(rv2)[j] == ru )
                    {    to_delete2.push_back( G2I.EdgeObjectIndexByIndexFrom( 
                              rv2, j ) );    }    }    }    }
     UniqueSort(to_delete2);
     G2I.DeleteEdges(to_delete2);

     // Remove ridiculously negative links.

     vec<int> to_removex;
     for ( int u1 = 0; u1 < nuni2; u1++ )
     {    int max_nlinks = 0;
          for ( int j = 0; j < G2I.From(u1).isize( ); j++ )
          {    const linklet& l = G2I.EdgeObjectByIndexFrom( u1, j );
               const int xmult = 5;
               const int xlinks = 5;
               if ( l.nlinks <= xlinks && l.sep < -K2-1 - xmult * l.dev )
               {    to_removex.push_back( 
                         G2I.EdgeObjectIndexByIndexFrom( u1, j ) );    }    }    }
     G2I.DeleteEdges(to_removex);

     // Test weak links.

     const int max_to_test_weak = 20;
     const int max_to_call_weak = 5;
     to_removex.clear( );
     #pragma omp parallel for
     for ( int u1 = 0; u1 < nuni2; u1++ )
     {    int max_nlinks = 0;
          for ( int j = 0; j < G2I.From(u1).isize( ); j++ )
          {    const linklet& l = G2I.EdgeObjectByIndexFrom( u1, j );
               max_nlinks = Max( max_nlinks, l.nlinks );    }
          if ( max_nlinks > max_to_test_weak ) continue;
          vec<Bool> supported( G2I.From(u1).size( ), False );
          for ( int j = 0; j < G2I.From(u1).isize( ); j++ )
          {    int u2 = G2I.From(u1)[j];
               const linklet& l = G2I.EdgeObjectByIndexFrom( u1, j );
               vec< vec< pair<int,int> > > walks1;
               int bad;
               vec<int> use;
               GetWalks( u1, u2, l.sep, l.dev, unibases2, K2, to_rc2, 
                    nexts2x, use, walks1, bad );
               if ( walks1.nonempty( ) ) supported[j] = True;    }
          if ( Sum(supported) == 0 ) continue;
          for ( int j = 0; j < G2I.From(u1).isize( ); j++ )
          {    if ( supported[j] ) continue;
               const linklet& l = G2I.EdgeObjectByIndexFrom( u1, j );
               if ( l.nlinks > max_to_call_weak ) continue;
               /*
               if (DELETE_VERBOSE)
               {    cout << "method X, deleting link from " << u1
                         << " to " << G2I.From(u1)[j] << "\n";    }
               */
               int u2 = G2I.From(u1)[j];
               int j2 = BinPosition( G2I.From( to_rc2[u2] ), to_rc2[u1] );
               ForceAssert( j2 >= 0 );
               #pragma omp critical
               {    to_removex.push_back( G2I.EdgeObjectIndexByIndexFrom( u1, j ) );
                    to_removex.push_back( G2I.EdgeObjectIndexByIndexFrom( 
                         to_rc2[u2], j2 ) );    }    }    }
     G2I.DeleteEdges(to_removex);    }

template void digraphE<linklet>::DeleteEdges(vec<int> const&);

vec<placementy> FindGenomicPlacementsY( const int u, basevector b, const int L,
     const vecbasevector& genome, const vec< vec< pair<int,int> > >& Glocs )
{    vec<placementy> places;
     for ( int pass = 1; pass <= 2; pass++ )
     {    if ( pass == 2 ) b.ReverseComplement( );
          int n = KmerId( b, L, 0 );
          for ( int z = 0; z < Glocs[n].isize( ); z++ )
          {    const basevector& g = genome[ Glocs[n][z].first ];
               const int gpos = Glocs[n][z].second;
               if ( !PerfectMatch( b, g, 0, gpos, b.size( ) ) ) continue;
               places.push( u, Glocs[n][z].first, 
                    gpos, gpos + b.isize( ), pass == 1, 0, 0, 0 );    }    }
     return places;    }

void RaiseLinks( const int K2, const vecbasevector& unibases1, 
     const vecbasevector& unibases2, const vec< vec< pair<int,int> > >& to_big, 
     const digraphE<linklet>& G, digraphE<linklet>& G2 )
{
     int nuni1 = unibases1.size( ), nuni2 = unibases2.size( );
     vec< triple<int,int,linklet> > links2;
     const int max_dist_mult = 15;
     const int max_self_dist_mult = 3;
     for ( int u1 = 0; u1 < nuni1; u1++ )
     {    if ( !to_big[u1].solo( ) ) continue;
          int v1 = to_big[u1][0].first, pos1 = to_big[u1][0].second;
          for ( int j = 0; j < G.From(u1).isize( ); j++ )
          {    int u2 = G.From(u1)[j];
               linklet l = G.EdgeObjectByIndexFrom( u1, j );
               if ( !to_big[u2].solo( ) ) continue;
               int v2 = to_big[u2][0].first, pos2 = to_big[u2][0].second;
               l.sep -= ( unibases2[v1].isize( ) - pos1 
                    - unibases1[u1].isize( )) + pos2;
               if ( l.sep < - K2 - 1 - max_dist_mult * l.dev ) continue;
               if ( v1 == v2 && l.sep < - K2 - 1 - max_self_dist_mult * l.dev ) 
                    continue;
               const int min_to_self_link = 20000;
               if ( v1 == v2 && unibases2[v1].isize( ) < min_to_self_link ) continue;
               links2.push( v1, v2, l );    }    }
     Sort(links2);
     vec< vec<int> > from(nuni2), to(nuni2);
     vec< vec<int> > from_edge_obj(nuni2), to_edge_obj(nuni2);
     vec<linklet> edges;
     for ( int i = 0; i < links2.isize( ); i++ )
     {    
          // In this version, for better or for worse, if there are multiple 
          // link groups, we choose the biggest one, rather than combine them.

          int j;
          for ( j = i + 1; j < links2.isize( ); j++ )
          {    if ( links2[j].first != links2[i].first ) break;
               if ( links2[j].second != links2[i].second ) break;    }
          int u1 = links2[i].first, u2 = links2[i].second;
          const linklet& l = links2[i].third;
          int sep = l.sep, dev = l.dev, nlinks = l.nlinks;
          from[u1].push_back(u2), to[u2].push_back(u1);
          from_edge_obj[u1].push_back( edges.size( ) );
          to_edge_obj[u2].push_back( edges.size( ) );
          edges.push( sep, dev, nlinks, 0 );
          i = j - 1;    }
     for ( int u = 0; u < nuni2; u++ )
     {    SortSync( from[u], from_edge_obj[u] );
          SortSync( to[u], to_edge_obj[u] );    }
     G2.Initialize( from, to, edges, to_edge_obj, from_edge_obj );    }

void Advance(
     // inputs:
     const int u, const vecbasevector& unibases2, const vec<int>& to_rc2,
     const vec<double>& raw2,
     const int K2, const digraphE<linklet> G,
     // heuristics:
     const int max_link_ratio,
     // logging:
     ostream& out,
     // outputs:
     vec<int>& unexts )
{
     // Clean links.  Suppose that when we look out from unibase u, we see two
     // unibases v1 and v2 (and possibly others).  Suppose that there are at least
     // 10 times as many links to v1 as to v2.  Suppose that there is a link from
     // v1 to v2.  Suppose that the distances imply that v1 and v2 have to overlap.
     // Then delete the link from u to v2.

     unexts.clear( );
     const int xdev_mult = 3;
     const int xdev_mult_len = 6;
     const int xmax_link_ratio = 10;
     vec<Bool> to_ignore( G.From(u).size( ), False );
     for ( int m1 = 0; m1 < G.From(u).isize( ); m1++ )
     for ( int m2 = 0; m2 < G.From(u).isize( ); m2++ )
     {    int v1 = G.From(u)[m1], v2 = G.From(u)[m2];
          const linklet& l1 = G.EdgeObjectByIndexFrom( u, m1 );
          const linklet& l2 = G.EdgeObjectByIndexFrom( u, m2 );
          if ( !( l1.nlinks >= xmax_link_ratio * l2.nlinks ) ) continue;
          int dev1 = xdev_mult * l1.dev, dev2 = xdev_mult * l2.dev;
          if ( l1.sep - dev1 + unibases2[v1].isize( ) < l2.sep + dev2 ) continue;
          if ( l2.sep - dev2 + unibases2[v2].isize( ) < l1.sep + dev1 ) continue;
          if ( Member( G.From(v1), v2 ) )
          {    int dev1_len = xdev_mult_len * l1.dev; 
               int dev2_len = xdev_mult_len * l2.dev;
               if ( l1.sep - dev1_len + unibases2[v1].isize( ) < l2.sep + dev2_len )
                    continue;
               if ( l2.sep - dev2_len + unibases2[v2].isize( ) < l1.sep + dev1_len )
                    continue;     }
          out << "(cleaning link from " << u << " to " << G.From(u)[m2] << ")\n";
          to_ignore[m2] = True;    }

     // Now do the actual cleaning.

     vec<int> linksx;
     vec<linklet> linksy;
     for ( int m = 0; m < G.From(u).isize( ); m++ )
     {    if ( to_ignore[m] ) continue;
          linksx.push_back( G.From(u)[m] );
          linksy.push_back( G.EdgeObjectByIndexFrom( u, m ) );    }

     // Announce links.

     for ( int m = 0; m < linksx.isize( ); m++ )
     {    int v = linksx[m];
          const linklet& l = linksy[m];
          out << "-- see " << v << "[" << unibases2[v].isize( ) - K2 + 1 << ",cn=" 
               << setiosflags(ios::fixed) << setprecision(2) << raw2[v] 
               << resetiosflags(ios::fixed) 
               << "], nlinks = " << l.nlinks << "\n";    }

     // Give up if no links.

     if ( linksx.empty( ) ) 
     {    out << "\ncan't see anything\n";
          out << "giving up\n";
          return;    }

     // Is there an ordering of the targets?

     restart0:
     vec< vec<int> > tried;
     restart:
     tried.push_back(linksx);
     Bool fail = False;
     for ( int m1 = 0; m1 < linksx.isize( ); m1++ )
     for ( int m2 = m1 + 1; m2 < linksx.isize( ); m2++ )
     {    int v1 = linksx[m1], v2 = linksx[m2];
          if ( Member( G.From(v2), v1 ) && !Member( G.From(v1), v2 ) )
          {    swap( linksx[m1], linksx[m2] );
               swap( linksy[m1], linksy[m2] );
               if ( Member( tried, linksx ) )
               {    fail = True;
                    goto test;    }
               tried.push_back(linksx);
               goto restart;    }    }
     test:
     if ( !fail )
     {    for ( int m = 0; m < linksx.isize( ) - 1; m++ )
          {    if ( to_rc2[ linksx[m] ] != linksx[m+1]
                    && !Member( G.From( linksx[m] ), linksx[m+1] ) )
               {    fail = True;
                    break;    }    }    }
     if ( !fail )
     {    
          // Identify and exclude links that probably go to unipaths of higher copy
          // number.  Also avoid cases where we link to both a unipath and its 
          // reverse complement, as there is likely a heterozygous inversion there, 
          // and we would prefer to link around it.

          const double max_cn_ratio = 1.5;
          double max_cn_ratio_to_use = max_cn_ratio;
          ratio_test:
          vec<Bool> exclude( linksx.size( ), False );
          for ( int m = 0; m < linksx.isize( ); m++ )
          {    int v = linksx[m];
               if ( raw2[u] > 0 && raw2[v] > 0 
                    && raw2[v]/raw2[u] > max_cn_ratio_to_use )
               {    exclude[m] = True;    }    }
          for ( int m1 = 0; m1 < linksx.isize( ); m1++ )
          for ( int m2 = m1+1; m2 < linksx.isize( ); m2++ )
          {    int v1 = linksx[m1], v2 = linksx[m2];
               if ( v2 == to_rc2[v1] )
               {    exclude[m1] = exclude[m2] = True;    }    }
          if ( Sum(exclude) == linksx.isize( ) )
          {    
               // Prevent infinite looping.  I'm not sure this is done correctly.
                    
               double M = 0;
               for ( int m = 0; m < linksx.isize( ); m++ )
               {    int v = linksx[m];
                    if ( raw2[u] > 0 && raw2[v] ) M = Max( M, raw2[v]/raw2[u] );    }
               if ( M <= max_cn_ratio_to_use ) goto delete_weakest;

               // Up the ratio.

               max_cn_ratio_to_use *= max_cn_ratio;
               goto ratio_test;    }
          int maxlinks = 0;
          for ( int m = 0; m < linksy.isize( ); m++ )
               if ( !exclude[m] ) maxlinks = Max( maxlinks, linksy[m].nlinks );
          int best = -1;
          for ( int m = 0; m < linksy.isize( ); m++ )
          {    if ( exclude[m] ) continue;
               if ( linksy[m].nlinks == maxlinks ) best = m;    }
          int unew = linksx[best];
          unexts.push_back(unew);
          out << "accepting " << unew << "\n";
          return;    }

     // Delete the least supported link and try again.

     delete_weakest:
     int min_links = 1000000000, max_links = 0;
     for ( int i = 0; i < linksy.isize( ); i++ )
     {    min_links = Min( min_links, linksy[i].nlinks );
          max_links = Max( max_links, linksy[i].nlinks );    }
     const int link_ratio_ignore = 10;
     const int link_ratio_ignore2 = 3;
     const int links_to_ignore = 5;
     if ( max_links >= link_ratio_ignore * min_links
          || ( max_links >= link_ratio_ignore2 * min_links
               && min_links <= links_to_ignore ) )
     {    vec<Bool> to_remove( linksx.size( ), False );
          for ( int i = 0; i < linksy.isize( ); i++ )
          {    if ( min_links == linksy[i].nlinks )
               {    to_remove[i] = True;
                    EraseIf( linksx, to_remove ), EraseIf( linksy, to_remove );
                    goto restart0;    }    }    }

     out << "\ncan't figure out what to delete\n";
     out << "visible:";
     for ( int m = 0; m < linksx.isize( ); m++ )
          out << " " << linksx[m];
     /*
     equiv_rel e( linksx.size( ) );
     for ( int i1 = 0; i1 < linksx.isize( ); i1++ )
     for ( int i2 = 0; i2 < linksx.isize( ); i2++ )
     {    if ( i2 != i1 )
          {    int u1 = linksx[i1], u2 = linksx[i2];
               if ( Member( G.From(u1), u2 ) ) e.Join( i1, i2 );    }    }
     vec<int> reps;
     e.OrbitReps(reps);
     for ( int j = 0; j < reps.isize( ); j++ )
     {    vec<int> o;
          e.Orbit( reps[j], o );
          // Should look at copy number!
          vec<int> nlinks( o.size( ) );
          for ( int i = 0; i < o.isize( ); i++ )
               nlinks[i] = linksy[ o[i] ].nlinks;
          ReverseSortSync( nlinks, o );
          unexts.push_back( linksx[ o[0] ] );    }
     Sort(unexts);
     out << "\nto use:";
     for ( int j = 0; j < unexts.isize( ); j++ )
          out << " " << unexts[j];
     */
     out << "\n\ngiving up\n";    }

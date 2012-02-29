///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "PairsManager.h"
#include "ParallelVecUtilities.h"
#include "Superb.h"
#include "VecUtilities.h"
#include "efasta/AmbiguityScore.h"
#include "efasta/EfastaTools.h"
#include "graph/DigraphTemplate.h"
#include "math/Functions.h"
#include "paths/ProcessGap.h"
#include "paths/LongReadTools.h"
#include "paths/Ulink.h"
#include "paths/Uniseq.h"
#include "paths/UnipathScaffold.h"

void GetWalks0( const int u1, const int u2, const int sep, const int dev,
     const vecbasevector& unibases, const int K, const vec<int>& to_rc,
     const vec< vec< pair<int,int> > >& nextsx, vec<int> use,
     vec< vec< pair<int,int> > >& walks1, int& bad )
{    
     // Heuristics.

     const int max_devs = 10;
     const double max_dev_diff = 3.0;
     const int max_iterations = 10000000;
     const int max_excess_charm = K;

     // Do two passes.  In the first pass, for now just exploratory, allow a given
     // vertex to be used at most twice.
     // FOR NOW TURNED OFF AS WE START WITH PASS 2.

     vec<double> offbys;
     for ( int pass = 1; pass <= 2; pass++ )
     {    
          double max_devs_this = max_devs;
          if ( pass == 2 )
          {    
               /*
               cout << "on screening pass, found " << walks1.size( )
                    << " walks, bad = " << bad << "\n";
               const int max_walks_to_print = 20;
               cout << "\nscreening walks from " << u1 << " to " << u2
                    << ", bad = " << bad
                    << ", nwalks = " << walks1.size( ) << "\n";
               for ( int j = 0; j < Min( max_walks_to_print, walks1.isize( ) ); j++ )
               {    cout << "[" << j << "] ";
                    const vec< pair<int,int> >& w = walks1[j];
                    vec<int> u, over;
                    for ( int m = 0; m < w.isize( ); m++ )
                    {    u.push_back( w[m].first );
                         if ( m > 0 ) over.push_back( w[m].second );    }
                    uniseq s( u, over );
                    s.Print( cout, K );
                    cout << "\n";    }
               if ( walks1.isize( ) > max_walks_to_print ) cout << "...\n";
               */
               if ( bad == 0 )
               {    use.clear( );
                    for ( int i = 0; i < walks1.isize( ); i++ )
                    {    for ( int j = 0; j < walks1[i].isize( ); j++ )
                              use.push_back( walks1[i][j].first );    }
                    UniqueSort(use);    }    }

          walks1.clear( ), offbys.clear( );
          vec< vec< pair<int,int> > > walks0;
          vec<int> charm0, charm1; // misnomer!
          vec< pair<int,int> > w;
          w.push( u1, 0 );
          walks0.push_back(w);
          charm0.push_back(0);
          int best_charm = 1000000000;
          bad = 0;
          int iterations = 0;
          while( walks0.nonempty( ) )
          {    if ( ++iterations > max_iterations )
               {    bad = 1;
                    break;    }
               vec< pair<int,int> > w = walks0.back( );
               walks0.pop_back( );
               int ch = charm0.back( );
               charm0.pop_back( );
               int v1 = w.back( ).first;
               int sepx = 0;
               for ( int j = 1; j < w.isize( ) - 1; j++ )
                    sepx += unibases[ w[j].first ].isize( ) - w[j].second;
               sepx -= w.back( ).second;
               double offby = double(sepx-sep) / double(dev);
               if ( iterations > 1 && offby > max_devs_this ) continue;
               if ( v1 == u2 && w.size( ) > 1 )
               {    if ( offby <= max_devs_this && offby >= -max_devs && 
                         ch <= best_charm + max_excess_charm )
                    {    walks1.push_back(w);
                         charm1.push_back(ch);
                         best_charm = Min( best_charm, ch );
                         offbys.push_back( Abs(offby) );    
                         max_devs_this 
                              = Min( max_devs_this, Abs(offby) + max_dev_diff );    }
                    continue;    }
               for ( int i = 0; i < nextsx[v1].isize( ); i++ )
               {    int v2 = nextsx[v1][i].first;
                    if ( use.nonempty( ) && !BinMember( use, v2 ) ) continue;
                    if ( pass == 1 )
                    {    int count = 0;
                         for ( int j = 0; j < w.isize( ); j++ )
                              if ( w[j].first == v2 ) count++;
                         if ( count > 1 ) continue;    }
                    int over = nextsx[v1][i].second;
                    int ch_new = ch + K - 1 - over;
                    if ( ch_new <= best_charm + max_excess_charm )
                    {    vec< pair<int,int> > wx(w);
                         wx.push( v2, over );
                         walks0.push_back(wx);    
                         charm0.push_back(ch_new);    }    }    }
          vec<Bool> to_delete( walks1.size( ), False );
          for ( int i = 0; i < walks1.isize( ); i++ )
               if ( charm1[i] > best_charm + max_excess_charm ) to_delete[i] = True;
          EraseIf( walks1, to_delete ), EraseIf( offbys, to_delete );    }

     // Select passing walks.

     SortSync( offbys, walks1 );
     if ( walks1.nonempty( ) )
     {    int x;
          for ( x = 1; x < walks1.isize( ); x++ )
               if ( offbys[x] > offbys[0] + max_dev_diff ) break;
          walks1.resize(x);    }    }

void GetWalks( const int u1, const int u2, const int sep, const int dev,
     const vecbasevector& unibases, const int K, const vec<int>& to_rc,
     const vec< vec< pair<int,int> > >& nextsx, const vec<int>& use,
     vec< vec< pair<int,int> > >& walks1, int& bad )
{
     GetWalks0( u1, u2, sep, dev, unibases, K, to_rc, nextsx, use, walks1, bad );
     if ( bad == 0 ) return;
     // vec< vec< pair<int,int> > > walks1_save = walks1;
     vec<int> use_rc;
     for ( int j = 0; j < use.isize( ); j++ )
          use_rc.push_back( to_rc[ use[j] ] );
     Sort(use_rc);
     GetWalks0( to_rc[u2], to_rc[u1], sep, dev, unibases, K, to_rc, nextsx, use_rc,
          walks1, bad );
     // if ( bad > 0 ) walks1 = walks1_save;
     // else
     {    for ( int j = 0; j < walks1.isize( ); j++ )
          {    walks1[j].ReverseMe( );
               for ( int l = 0; l < walks1[j].isize( ); l++ )
                    walks1[j][l].first = to_rc[ walks1[j][l].first ];
               for ( int l = walks1[j].isize( ) - 1; l >= 1; l-- )
                    walks1[j][l].second = walks1[j][l-1].second;
               walks1[j][0].second = 0;     }    }    }

// WalksToPatches: here the patch is defined to start with the first overlapping
// base and extend through the last overlapping base, as in the following pictures:
//
//    (1)          ************
//    ---------------------
//                 ------------
//                     ---------------------------------
//
//
//    (2)          *******************************
//    ---------------------
//                 -------------------------------
//                                         ---------------------------------
//

void WalksToPatches( const vec< vec< pair<int,int> > >& walks1,
     const vecbasevector& unibases, vec<basevector>& patches )
{
     patches.resize( walks1.size( ) );
     for ( int j = 0; j < walks1.isize( ); j++ )
     {    const vec< pair<int,int> >& w = walks1[j];
          for ( int pass = 1; pass <= 2; pass++ )
          {    int pos = 0;
               for ( int i = 1; i < w.isize( ); i++ )
               {    int to = w[i].second;
                    if ( i < w.isize( ) - 1 )
                         to = unibases[ w[i].first ].isize( ) - w[i+1].second;
                    if ( pass == 2 )
                    {    for ( int l = 0; l < to; l++ )
                              patches[j].Set(pos++, unibases[ w[i].first ][l]);    }
                    else pos += to;    }
               if ( pass == 1 ) patches[j].resize(pos);    }    }    }

void FilterWalksUsingJumps2( const int edge_id, snark& S, const vec<int>& to_left,
     const vec<int>& to_right, const vecbasevector& jbases,
     const PairsManager& jpairs,
     const vec< vec< triple<int,int,Bool> > >& placements_by_read,
     const vec< vec< triple<int64_t,int,Bool> > >& placements_by_unipath,
     const int verbosity )
{
     const gapster& g = S.Edge(edge_id);
     ForceAssert( g.Closed( ) );
     if ( verbosity >= 1 )
     {    cout << "\nfiltering walks from " << g.Closure(0).U( ).front( )
               << " to " << g.Closure(0).U( ).back( ) << "\n";    }

     // See if edge_id has an rc twin.

     int edge_id_rc = -1, x1 = to_left[edge_id], x2 = to_right[edge_id];
     uniseq v1rc = S.Vert(x1).Reverse( S.ToRc( ) );
     uniseq v2rc = S.Vert(x2).Reverse( S.ToRc( ) );
     for ( int e = 0; e < S.EdgeN( ); e++ )
     {    if ( S.Vert( to_left[e] ) == v2rc && S.Vert( to_right[e] ) == v1rc )
               edge_id_rc = e;    }

     vec<int> delta;
     vec< vec<int> > us( g.ClosureCount( ) );
     for ( int i = 0; i < g.ClosureCount( ); i++ )
     {    for ( int j = 0; j < g.Closure(i).N( ); j++ )
               us[i].push_back( g.Closure(i).U(j) );
          UniqueSort( us[i] );    }
     vec< vec<int> > uss(us);
     ParallelUniqueSort(uss);

     vec<int> uall;
     for ( int i = 0; i < uss.isize( ); i++ )
     {    for ( int j = 0; j < uss[i].isize( ); j++ )
               uall.push_back( uss[i][j] );    }
     UniqueSort(uall);
     #pragma omp parallel for
     for ( int i = 0; i < uall.isize( ); i++ )
     {    int u = uall[i];
          Bool missing = False;
          for ( int j = 0; j < uss.isize( ); j++ )
          {    if ( !BinMember( uss[j], u ) )
               {    missing = True;
                    break;    }    }
          if (missing)
          {   
               #pragma omp critical
               {    delta.push_back(u);    }    }    }
               
     UniqueSort(delta);
     vec<int64_t> pids;
     for ( int i = 0; i < delta.isize( ); i++ )
     {    int u = delta[i];
          for ( int j = 0; j < placements_by_unipath[u].isize( ); j++ )
          {    int64_t id = placements_by_unipath[u][j].first;
               pids.push_back( jpairs.getPairID(id) );    }    }
     UniqueSort(pids);

     // Compute u_pos_apos.

     vec< vec< triple<int,int,int> > > u_pos_apos( S.EdgeN( ) );
     for ( int gi = 0; gi < S.EdgeN( ); gi++ )
     {    const gapster& g = S.Edge(gi);
          if ( g.Open( ) ) continue;

          // Find all the unibases in g, along with their positions and
          // antipositions.

          const vec<uniseq>& closures = g.Closures( );
          for ( int cj = 0; cj < closures.isize( ); cj++ )
          {    const uniseq& r = closures[cj];
               int len = r.Len( ), pos = 0;
               for ( int j = 0; j < r.N( ); j++ )
               {    int u = r.U(j);
                    if ( j > 0 )
                         pos += S.Unibase( r.U(j-1) ).isize( ) - r.Over(j-1);
                    int antipos = len - pos;
                    u_pos_apos[gi].push( u, pos, antipos );    }    }
          UniqueSort( u_pos_apos[gi] );    }

     vec< vec<int> > HITS( g.ClosureCount( ) );
     vec<String> reports( pids.size( ) );

     const int dev_mult = 3;
     vec<int> min_len, max_len;
     for ( int e = 0; e < S.EdgeN( ); e++ )
     {    min_len.push_back( S.Edge(e).MinLen(dev_mult) );
          max_len.push_back( S.Edge(e).MaxLen(dev_mult) );    }

     #pragma omp parallel for
     for ( int i = 0; i < pids.isize( ); i++ )
     {    int64_t pid = pids[i];
          vec< vec<placement_on> > pos_on;
          vec< triple<int,int,Bool> > placements;
          GetPairPlacements( S, to_right, min_len, max_len, pid, jbases, jpairs, 
               placements_by_read, u_pos_apos, pos_on, placements );
          Bool found_alt_placement = False, found_good_placement = False;
          Bool found_double_placement = False;
          for ( int j = 0; j < placements.isize( ); j++ )
          {    int j1 = placements[j].first, j2 = placements[j].second;
               int oid1 = pos_on[0][j1].id, oid2 = pos_on[1][j2].id;
               if ( edge_id + S.VertN( ) != oid1 && edge_id + S.VertN( ) != oid2
                    && ( edge_id_rc < 0 || edge_id_rc + S.VertN( ) != oid1 )
                    && ( edge_id_rc < 0 || edge_id_rc + S.VertN( ) != oid2 ) )
               {    found_alt_placement = True;
                    break;    }
               /*
               if ( edge_id + S.VertN( ) == oid1 && edge_id + S.VertN( ) != oid2 )
                    found_good_placement = True;
               if ( edge_id + S.VertN( ) != oid1 && edge_id + S.VertN( ) == oid2 )
                    found_good_placement = True;
               */

               // Call it good if at least one end lands on the edge:

               Bool on1 = ( edge_id + S.VertN( ) == oid1 );
               Bool on2 = ( edge_id + S.VertN( ) == oid2 );
               if ( on1 || on2 ) found_good_placement = True;
               if ( on1 && on2 ) found_double_placement = True;

               /*
               if ( edge_id + S.VertN( ) == oid1 || edge_id + S.VertN( ) == oid2 )
                    found_good_placement = True;
               if ( edge_id + S.VertN( ) == oid1 && edge_id + S.VertN( ) == oid2 )
                    found_double_placement = True;    
               */
                    }
          if ( verbosity <= 1 &&
               ( found_alt_placement || !found_good_placement ) )
          {     continue;   }
          if ( verbosity >= 1 )
          {    ostringstream out;
               int64_t id1 = jpairs.ID1(pid), id2 = jpairs.ID2(pid);
               for ( int pass = 0; pass < 2; pass++ )
               {    int64_t id = ( pass == 0 ? id1 : id2 );
                    out << "\nplacements of " << ( pass == 0 ? "first" : "second" ) 
                         << " read " << ( pass == 0 ? id1 : id2 ) << " of pair " 
                         << pid << "\n";
                    for ( int j = 0; j < placements_by_read[id].isize( ); j++ )
                    {    int u = placements_by_read[id][j].first;
                         int pos = placements_by_read[id][j].second;
                         Bool fw = placements_by_read[id][j].third;
                         out << "[" << j+1 << "] " << u << "." << pos << "-" 
                              << (fw ? "fw" : "rc" ) << "\n";    }    }
               for ( int pass = 0; pass < 2; pass++ )
               {    out << "\npos_on for " << ( pass == 0 ? "first" : "second" ) 
                         << " read " << ( pass == 0 ? id1 : id2 ) << " of pair " 
                         << pid << "\n";
                    for ( int j = 0; j < pos_on[pass].isize( ); j++ )
                    {    const placement_on& p = pos_on[pass][j];
                         out << "[" << j+1 << "] ";
                         pos_on[pass][j].Print( 
                              out, S, to_left, to_right );    
                         out << "\n";    }    }
               out << "\nplacements of pair " << pid << "\n";
               for ( int j = 0; j < placements.isize( ); j++ )
               {    int j1 = placements[j].first, j2 = placements[j].second;
                    int oid1 = pos_on[0][j1].id, oid2 = pos_on[1][j2].id;
                    out << "[" << j+1 << "] left on ";
                    pos_on[0][j1].Print( out, S, to_left, to_right );
                    out << ", right on ";
                    pos_on[1][j2].Print( out, S, to_left, to_right );
                    out << "\n";    }    
               reports[i] = out.str( );    }
          if ( found_alt_placement || !found_good_placement ) continue;

          // Second version rather primitive, not really right:

          vec<Bool> supported( g.ClosureCount( ), False );

          for ( int j = 0; j < g.ClosureCount( ); j++ )
          {    int64_t id1 = jpairs.ID1(pid), id2 = jpairs.ID2(pid);
               for ( int pass = 1; pass <= 2; pass++ )
               {    if ( pass == 2 ) swap( id1, id2 );
                    Bool found1 = False, found2 = False;
                    for ( int l1 = 0; 
                         l1 < placements_by_read[id1].isize( ); l1++ )
                    {    if ( !placements_by_read[id1][l1].third ) continue;
                         for ( int m = 0; m < g.Closure(j).U( ).isize( ); m++ )
                         {    if ( placements_by_read[id1][l1].first
                                   == g.Closure(j).U(m) )
                              {    found1 = True;
                                   break;    }    }
                         if (found1) break;    }
                    for ( int l2 = 0; 
                         l2 < placements_by_read[id2].isize( ); l2++ )
                    {    if ( placements_by_read[id2][l2].third ) continue;
                         for ( int m = 0; m < g.Closure(j).U( ).isize( ); m++ )
                         {    if ( placements_by_read[id2][l2].first
                                   == g.Closure(j).U(m) )
                              {    found2 = True;
                                   break;    }    }
                         if (found2) break;    }
                    if (found_double_placement)
                    {    if ( !found1 || !found2 ) continue;    }
                    else
                    {    if ( !found1 && !found2 ) continue;    }
                    supported[j] = True;
                    break;    }    }
     
          /*
          for ( int j = 0; j < delta.isize( ); j++ )
          {    int u = delta[j];
               Bool have_pid = False;
               for ( int j = 0; j < placements_by_unipath[u].isize( ); j++ )
               {    int64_t id = placements_by_unipath[u][j].first;
                    int64_t pidx = jpairs.getPairID(id);
                    if ( pidx == pid ) 
                    {    have_pid = True;
                         break;    }    }
               if (have_pid)
               {    for ( int j = 0; j < g.ClosureCount( ); j++ )
                         if ( BinMember( us[j], u ) ) supported[j] = True;    }    }
          */

          #pragma omp critical
          {    for ( int k = 0; k < g.ClosureCount( ); k++ )
                    if ( supported[k] ) HITS[k].push_back(i);    }    }
     for ( int i = 0; i < HITS.isize( ); i++ )
          Sort( HITS[i] );
     if ( verbosity >= 1 )
     {    for ( int i = 0; i < reports.isize( ); i++ )
               cout << reports[i];    }

     // Make decisions.  I'm not sure that the parallel for loop is deterministic.

     if ( verbosity >= 1 )
     {    cout << "\n";
          for ( int i = 0; i < g.ClosureCount( ); i++ )
          {    for ( int j = 0; j < HITS[i].isize( ); j++ )
               {    int k = HITS[i][j];
                    int64_t pid = pids[k];
                    cout << "walk " << i << " supported by jump pair " << pid 
                         << ", id1 = " << jpairs.ID1(pid) << ", id2 = " 
                         << jpairs.ID2(pid) << "\n";    }    }    }
     vec<Bool> to_delete( g.ClosureCount( ), False );
     const int min_mult = 4;
     const int min_count = 2;
     vec< vec<Bool> > h( g.ClosureCount( ), vec<Bool>( pids.size( ), False ) );
     for ( int i = 0; i < g.ClosureCount( ); i++ )
     {    for ( int j = 0; j < HITS[i].isize( ); j++ )
               h[i][ HITS[i][j] ] = True;    }
     vec<int> ids( g.ClosureCount( ), vec<int>::IDENTITY );
     ParallelSortSync( h, ids );
     vec< vec<Bool> > hred;
     vec< vec<int> > ids_hred;
     vec<int> nhred;
     for ( int i = 0; i < h.isize( ); i++ )
     {    int j = h.NextDiff(i);
          hred.push_back( h[i] );
          nhred.push_back( Sum( h[i] ) );
          vec<int> x;
          for ( int k = i; k < j; k++ )
               x.push_back( ids[k] );
          ids_hred.push_back(x);
          i = j - 1;    }
     vec<Bool> hred_to_delete( hred.size( ), False );
     #pragma omp parallel for
     for ( int i1 = 0; i1 < hred.isize( ); i1++ )
     {    if ( hred_to_delete[i1] ) continue;
          for ( int i2 = 0; i2 < hred.isize( ); i2++ )
          {    if ( hred_to_delete[i2] ) continue;
               if ( !( nhred[i1] > nhred[i2] ) ) continue;
               int better1 = 0, better2 = 0;
               for ( int j = 0; j < pids.isize( ); j++ )
               {    if ( hred[i1][j] && !hred[i2][j] ) better1++;
                    if ( hred[i2][j] && !hred[i1][j] ) better2++;    }
               if ( better1 >= min_count && better1 >= min_mult * better2 ) 
               {    if ( verbosity >= 1 )
                    {    cout << "walks";
                         for ( int j = 0; j < ids_hred[i2].isize( ); j++ )
                              cout << " " << ids_hred[i2][j];
                         cout << " are beaten by walks";
                         for ( int j = 0; j < ids_hred[i1].isize( ); j++ )
                              cout << " " << ids_hred[i1][j];
                         cout << "\n";    }
                    hred_to_delete[i2] = True;
                    #pragma omp critical
                    {    for ( int j = 0; j < ids_hred[i2].isize( ); j++ )
                              to_delete[ ids_hred[i2][j] ] = True;   }   }   }   }

     // Also remove closures having identical bases.

     vec<basevector> cbases;
     for ( int i = 0; i < g.ClosureCount( ); i++ )
          cbases.push_back( g.Closure(i).Bases( ) );
     #pragma omp parallel for
     for ( int i1 = 0; i1 < g.ClosureCount( ); i1++ )
     {    if ( to_delete[i1] ) continue;
          for ( int i2 = i1 + 1; i2 < g.ClosureCount( ); i2++ )
          {    if ( to_delete[i2] ) continue;
               if ( cbases[i1] == cbases[i2] ) 
               {    
                    #pragma omp critical
                    {    to_delete[i2] = True;    }    }    }    }

     // Make edit.

     S.Gmutable( ).EdgeObjectMutable(edge_id).RemoveSomeClosures(to_delete);    }

void FilterWalksUsingJumps( const int u1, const int u2,
     const vec<basevector>& jbases_sorted, const vec<int64_t>& jbases_sorted_id,
     const PairsManager& jpairs, const vec< triple<int64_t,int,int> >& jaligns,
     const vecbasevector& unibases, vec<basevector>& patches,
     vec< vec< pair<int,int> > >& walks1 )
{
     int npatches = patches.size( );
     if ( jbases_sorted.empty( ) || npatches < 2 ) return;
     const int verbosity = 0;
     int L = jbases_sorted[0].size( );
     vec<int> hits( npatches, 0 ), id( npatches, vec<int>::IDENTITY );
     vec< vec<int> > HITS(npatches);
     #pragma omp parallel for
     for ( int i = 0; i < npatches; i++ )
     {    for ( int j = 0; j <= patches[i].isize( ) - L; j++ )
          {    basevector b( patches[i], j, L );
               for ( int pass = 1; pass <= 2; pass++ )
               {    if ( pass == 2 ) b.ReverseComplement( );
                    int64_t low = LowerBound( jbases_sorted, b );
                    int64_t high = UpperBound( jbases_sorted, b );
                    for ( int64_t x = low; x < high; x++ )
                    {    int64_t id = jbases_sorted_id[x];
                         int64_t idp = jpairs.getPartnerID(id);
                         int xp = BinPosition1( jaligns, idp );
                         if ( xp < 0 ) continue;
                         int mx = jaligns[xp].second, p = jaligns[xp].third;
                         const int max_dist = 10000; // ****************************
                         if ( mx == u1 && unibases[u1].isize( )
                              - p <= max_dist && pass == 2 )
                         {    HITS[i].push_back(id);    }
                         if ( mx == -u2-1 && p <= max_dist && pass == 1 )
                              HITS[i].push_back(id);    }    }    }    }
     for ( int i = 0; i < npatches; i++ )
          UniqueSort( HITS[i] );
     Bool verbose = False;
     if (verbose)
     {    for ( int i = 0; i < npatches; i++ )
          {    for ( int j = 0; j < HITS[i].isize( ); j++ )
               {    cout << "walk " << i << " supported by jump read "
                         << HITS[i][j] << "\n";    }    }    }
     vec<Bool> to_delete( npatches, False );
     for ( int i1 = 0; i1 < npatches; i1++ )
     {    if ( to_delete[i1] ) continue;
          for ( int i2 = 0; i2 < npatches; i2++ )
          {    if ( !( HITS[i1].size( ) > HITS[i2].size( ) ) ) continue;
               int better1 = 0, better2 = 0;
               for ( int j = 0; j < HITS[i1].isize( ); j++ )
               {    if ( !BinMember( HITS[i2], HITS[i1][j] ) ) better1++;    }
               for ( int j = 0; j < HITS[i2].isize( ); j++ )
               {    if ( !BinMember( HITS[i1], HITS[i2][j] ) ) better2++;    }
               const int min_mult = 4;
               const int min_count = 2;
               if ( better1 >= min_count && better1 >= min_mult * better2 ) 
                    to_delete[i2] = True;    }    }
     EraseIf( walks1, to_delete );
     EraseIf( patches, to_delete );    }

void FilterWalksByIllegalCN1( const superb& s,
     vec<basevector>& patches, vec< vec< pair<int,int> > >& walks1 )
{    
     vec<int> s_used;
     for ( int j = 0; j < s.Ntigs( ); j++ )
          s_used.push_back( s.Tig(j) );
     UniqueSort(s_used);
     vec<Bool> to_delete( patches.size( ), False );
     for ( int j = 0; j < patches.isize( ); j++ )
     {    for ( int l = 1; l < walks1[j].isize( ) - 1; l++ )
          {    if ( BinMember( s_used, walks1[j][l].first ) )
                    to_delete[j] = True;    }    }
     EraseIf( walks1, to_delete );
     EraseIf( patches, to_delete );    }

Bool Follow( const placementy& p1, const placementy& p2, const int K )
{
     if ( p2.g != p1.g || p2.fw != p1.fw ) return False;
     if ( p2.fw )
     {    if ( p1.Pos - p2.pos != K - 1 ) return False;    }
     else if ( p2.Pos - p1.pos != K - 1 ) return False;
     return True;    }

void Print( ostream& out, const placementy& p, const int len )
{    out << " " << p.g << "." << p.pos << "-" << p.Pos
          << " " << ( p.fw ? "fw" : "rc" );    }

void ProcessGap(

     // definition of the gap

     const superb& s, const int u1, const int u2, const int sep, const int dev,

     // global structures

     const vecbasevector& unibases, const vec<int>& to_rc, 
     const vec< vec< pair<int,int> > >& nextsx, 
     const vec< vec< pair<int,int> > >& nextsy, const int K,
     const vec<int>& use,

     // jump stuff to filter

     const vec<basevector>& jbases_sorted, const vec<int64_t>& jbases_sorted_id,
     const PairsManager& jpairs, const vec< triple<int64_t,int,int> >& jaligns,
     const vec< vec< triple<int,int,Bool> > >& placements_by_read,

     // for evaluation and logging

     ostream& out, const int verbosity, const Bool validate, const int LG, 
     const vecbasevector& genome, const vec< vec< pair<int,int> > >& Glocs,
     Bool& have_lastp, placementy& lastp, int& gapcount, const digraphE<linklet>& G,

     // output

     Bool& found_patch, efasta& epatch, int& lastover,
     vec<int>& extras    // unibases1 used in every patch

)

{    
     // Get walks, then the sequences associated to the patches.  Filter walks that 
     // would violate copy number 1 nonredundancy.  Filter walks using jump 
     // alignments.

     vec< vec< pair<int,int> > > walks1;
     int bad;
     int walk_attempts = 1;
     GetWalks( u1, u2, sep, dev, unibases, K, to_rc, nextsy, use, walks1, bad );
     vec<basevector> patches;
     if ( bad == 0 ) 
     {    WalksToPatches( walks1, unibases, patches );
          FilterWalksByIllegalCN1( s, patches, walks1 );
          FilterWalksUsingJumps( u1, u2, jbases_sorted, jbases_sorted_id, jpairs, 
               jaligns, unibases, patches, walks1 );    }

     // Try again if no walks, this time using nextsx.

     if ( bad == 0 && walks1.empty( ) )
     {    walk_attempts = 2;
          GetWalks( u1, u2, sep, dev, unibases, K, to_rc, nextsx, use, walks1, bad );
          if ( bad == 0 )
          {    WalksToPatches( walks1, unibases, patches );
               FilterWalksByIllegalCN1( s, patches, walks1 );
               FilterWalksUsingJumps( u1, u2, jbases_sorted, jbases_sorted_id, 
                    jpairs, jaligns, unibases, patches, walks1 );    }    }

     // Process walks.

     String single_line = "-------------------------------------------------------"
          "-----------------------------\n";
     if ( bad == 0 && walks1.solo( ) )
     {    const vec< pair<int,int> >& w = walks1[0];
          for ( int j = 1; j < w.isize( ); j++ )
               extras.push_back( w[j].first );
          if ( verbosity >= 1 )
          {    for ( int j = 1; j < w.isize( ) - 1; j++ )
               {    int v = w[j].first;
                    if ( j > 1 ) out << "+";
                    out << v;    }    }
          int nkmers = patches[0].isize( ) - (K-1);
          ForceAssertGe( w.isize( ), 2 );
          found_patch = True;
          if ( w.size( ) > 2 )
          {    if ( verbosity >= 1 ) out << "[" << nkmers << "]";
               if ( validate && verbosity >= 1 )
               {    out << " (";
                    vec<placementy> goods = FindGenomicPlacementsY( 
                         0, patches[0], LG, genome, Glocs );
                    Bool have_good = False;
                    if (have_lastp)
                    {    for (int l = 0; l < goods.isize( ); l++)
                         {    if ( Follow( lastp, goods[l], w[1].second + 1 ) )
                              {    have_good = True;
                                   goods[0] = goods[l];
                                   goods.resize(1);
                                   break;    }    }    }
                    Bool print_line = ( have_lastp && !have_good );
                    for ( int l = 0; l < goods.isize( ); l++ )
                         Print( out, goods[l], patches[0].size( ) );
                    if ( goods.solo( ) )
                    {    have_lastp = True;
                         lastp = goods[0];
                         lastover = w.back( ).second;    }
                    else have_lastp = False;
                    out << " )\n";    
                    if (print_line) out << single_line;    }
               if ( !validate && verbosity >= 1 ) out << "\n";
               found_patch = True;
               epatch = patches[0].ToString( );    }    }
     else 
     {    if (validate) have_lastp = False;
          if ( verbosity >= 1 )
          {    out << "\n[gap " << ++gapcount << "] ";
               int nlinks = -1;
               for ( int j = 0; j < G.From(u1).isize( ); j++ )
               {    if ( G.From(u1)[j] == u2 )
                    {    nlinks = G.EdgeObjectByIndexFrom( u1, j ).nlinks;    }    }
               PRINT3_TO( out, sep, dev, nlinks );
               PRINT3_TO( out, walk_attempts, walks1.size( ), bad );    }
          if ( bad != 0 || walks1.empty( ) )
          {    found_patch = False;
               return;    }
          if ( verbosity >= 1 ) out << "walks:\n";
          for ( int z = 0; z < walks1.isize( ); z++ )
          {    const vec< pair<int,int> >& w = walks1[z];
               if ( verbosity >= 1 && z <= 5 )
               {    out << "<" << z << "> ";
                    for ( int j = 1; j < w.isize( ) - 1; j++ )
                    {    int v = w[j].first;
                         if ( j > 1 ) out << "+";
                         out << v;    }
                    out << "[" << patches[z].isize( ) - (K-1) << "]";    }
               if ( validate && verbosity >= 1 && z <= 5 )
               {    out << " (";
                    vec<placementy> goods = FindGenomicPlacementsY( 
                         z, patches[z], LG, genome, Glocs );
                    for ( int l = 0; l < goods.isize( ); l++ )
                         Print( out, goods[l], patches[z].size( ) );
                    out << " )";    }
               if ( verbosity >= 1 && z <= 5 ) out << "\n";
               if ( verbosity >= 1 && z == 5 && z < walks1.isize( ) - 1 )
                    out << "...\n";    }

          // Form efasta.  First get shared sequence in patches.

          int npatches = patches.size( );
          int left_share, right_share;
          GetShares( patches, left_share, right_share );
          basevector Left_share( patches[0], 0, left_share );
          basevector Right_share( patches[0], 
               patches[0].isize( ) - right_share, right_share );
     
          // Find the variant part "patches_diff" of the patches.  After this,
          // the joined contigs may be represented as
          // tigse[m1] 
          // + Left_share + {patches_diff} + Right_share 
          // + tigse[m2] (after deleting over2+K bases on left).
          
          vec<basevector> patches_diff(npatches);
          for ( int i = 0; i < npatches; i++ )
          {    patches_diff[i].SetToSubOf( patches[i], left_share, 
                    patches[i].isize( ) - left_share - right_share );    }
          vec<int> varsize(npatches);
          for ( int i = 0; i < npatches; i++ )
               varsize[i] = patches_diff[i].size( );
     
          // Generate report.
     
          Sort(varsize);
          String varsizes = "{";
          for ( int i = 0; i < npatches; i++ )
          {    if ( i > 0 ) varsizes += ",";
               varsizes += ToString( varsize[i] );    
               int j;
               for ( j = i + 1; j < npatches; j++ )
                    if ( varsize[j] != varsize[i] ) break;
               if ( j - i > 1 ) varsizes += "^" + ToString(j-i);
               i = j - 1;    }
          varsizes += "}";
          String share = "{" + ToString(left_share) + ","
               + ToString(right_share) + "}";
          if ( verbosity >= 1 ) PRINT2_TO( out, share, varsizes );

          // Print patches.

          for ( int i = 0; i < npatches; i++ )
          {    if ( verbosity >= 2 )
               {    patches_diff[i].Print( 
                         out, "patches_diff[" + ToString(i) + "]" );    }    }

          // Test for unacceptable patch.  This test really doesn't make sense.

          const int max_var_sum = 500;
          const int max_score = 150;
          const int64_t max_cost = 100 * 1000 * 1000;
          int64_t cost = AmbiguityScoreCost(patches_diff);
          if ( cost > max_cost )
          {    found_patch = False;
               if ( verbosity >= 1 ) out << "efasta patch rejected, cost too high\n";
               return;    }
          int score = AmbiguityScore(patches_diff);
          if ( verbosity >= 1 ) out << "ambiguity score = " << score << "\n";
          if ( score > max_score ) 
          {    found_patch = False;
               if ( verbosity >= 1 ) 
                    out << "efasta patch rejected, score too high\n";
               return;    }
          for ( int i = 0; i < walks1.isize( ); i++ )
          {    for ( int j = 1; j < walks1[i].isize( ) - 1; j++ )
                    extras.push_back( walks1[i][j].first );    }

          // Build and filter efasta patches.  The current filter is very stringent.

          epatch = Left_share.ToString( ) + ( npatches > 1 ? "{" : "" );
          for ( int i = 0; i < npatches; i++ )
               epatch += ( i > 0 ? "," : "" ) + patches_diff[i].ToString( );
          epatch += ( npatches > 1 ? "}" : "" ) + Right_share.ToString( );
          found_patch = True;
          if ( verbosity >= 1 ) out << "efasta patch generated\n\n";    }    }

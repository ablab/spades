///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#undef NDEBUG

#include <omp.h>

#include "Basevector.h"
#include "CoreTools.h"
#include "graph/Digraph.h"
#include "graph/DigraphTemplate.h"
#include "paths/CorrectLongReads1.h"
#include "paths/KmerAlignSet.h"
#include "paths/LongReadTools.h"
#include "paths/Uniseq.h"

void Phase1( 

     // inputs:

     const int NUM_THREADS, const vecbasevector& unibases, const vec<int>& to_rc,
     const vec< vec<int> >& nexts, const int K, const int L,
     const vec< vec< pair<int,int> > >& Ulocs, const vecbasevector& longreads,
     const vec<int>& reads_to_use, const heuristics& heur, 

     // algorithms:

     const Bool USE_SHORTEST, const Bool CLUSTER_ALIGNS_NEW_CLEAN,

     // control:

     const Bool PATCHES_ONLY, const int MIN_PATCH1, const int MIN_PATCH2,
     const Bool STANDARD_ALIGNS, const Bool FILTER,
     const Bool CLEAN_GRAPH, const Bool CYCLIC_BAIL,
     const Bool KILL_INFERIOR_CLUSTERS_NEW, const int MIN_SPREAD1,

     // logging:

     const vec<int>& reads_to_trace, const Bool VERBOSE1, const Bool SKIP_SILENT, 
     const Bool DOT1, const Bool QLT1, const String& data_dir, 
     const String& run_dir, const Bool ABBREVIATE_ALIGNMENTS,
     const Bool PRINT_RAW_ALIGNS,

     // outputs:

     vec<uniseq>& UNISEQ, vec< vec<int> >& UNISEQ_ID, vec<Bool>& COMPUTED, 
     vecbasevector& all, vec<GapPatcher>& patchers, vec< digraphVE<int,int> >& Hall )
{
     // Set up.

     int nreads = longreads.size( );
     vec<Bool> done( nreads, False );
     int done_ptr = 0;
     vec<String> reports( longreads.size( ) );
     if (QLT1) omp_set_num_threads(1);

     // Phase 1.  Now go through the reads.

     #pragma omp parallel for schedule(dynamic, 1)
     for ( int idx = 0; idx < reads_to_use.isize( ); idx++ )
     {    int id = reads_to_use[idx];
          const basevector& x = longreads[id];

          // Find the L-mer matches of the read to the unibases.

          vec< triple<int,int,int> > aligns, aligns2;
          KmerAlignSet alx( x, L, Ulocs );
          for ( int i = 0; i < alx.NAligns( ); i++ )
          {    for ( int j = 0; j < alx.Count(i); j++ )
                    aligns.push( alx.U(i), alx.Offset(i,j), alx.Rpos(i,j) );    }

          if (PRINT_RAW_ALIGNS)
          {    ostringstream hout;
               hout << "\nraw alignments of read " << id << "\n";
               for ( int j = 0; j < aligns.isize( ); j++ )
               {    int u = aligns[j].first, offset = aligns[j].second;
                    int rpos = aligns[j].third;
                    hout << u << " " << rpos << " " << offset << "\n";    }
               #pragma omp critical
               {    cout << hout.str( );    }    }

          // Filter the matches.

          if ( !FILTER ) aligns2 = aligns;
          else
          {    for ( int i = 0; i < aligns.isize( ); i++ )
               {    int offset = aligns[i].second, rpos = aligns[i].third;
                    int upos = rpos + offset;
                    int u = aligns[i].first;
                    int ulen = unibases[u].size( );
                    int overlap = Max( Min( rpos, upos ),
                         Min( longreads[id].isize( ) - rpos, ulen - upos ) );
                    Bool confirmed = True;
                    if ( overlap >= heur.min_overlap_to_see_other )
                    {    confirmed = False;
                         for ( int j = i - 1; j >= 0; j-- )
                         {    if ( aligns[j].first != u ) break;
                              int diff = Abs( offset - aligns[j].second );
                              if ( diff > heur.max_offset_diff_for_other ) break;
                              if ( Abs( rpos - aligns[j].third ) 
                                   < heur.min_rdist_for_other )
                              {    continue;    }
                              confirmed = True;
                              break;    }
                         for ( int j = i + 1; j < aligns.isize( ); j++ )
                         {    if ( aligns[j].first != u ) break;
                              int diff = Abs( offset - aligns[j].second );
                              if ( diff > heur.max_offset_diff_for_other ) break;
                              if ( Abs( rpos - aligns[j].third ) 
                                   < heur.min_rdist_for_other )
                              {    continue;    }
                              confirmed = True;
                              break;    }
                         if ( !confirmed ) continue;    }
                    if (confirmed) aligns2.push_back( aligns[i] );    }    }

          // Cluster alignments.
          
          KmerAlignSet aligns_in, aligns_out;
          for ( int i = 0; i < aligns2.isize( ); i++ )
          {    int u = aligns2[i].first, j;
               for ( j = i + 1; j < aligns2.isize( ); j++ )
                    if ( aligns2[j].first != u ) break;
               vec< pair<int,int> > rpos_upos;
               for ( int k = i; k < j; k++ )
               {    rpos_upos.push( aligns2[k].third, 
                         aligns2[k].second + aligns2[k].third );    }
               aligns_in.AddAlign( u, rpos_upos );
               i = j - 1;    }
          ClusterAlignsNew( aligns_in, aligns_out, CLUSTER_ALIGNS_NEW_CLEAN,
               MIN_SPREAD1 );
          vec< pair< int, vec< pair<int,int> > > > ALIGNS = aligns_out.X( );

          // For efficiency, test to see if the read is unambiguously aligned
          // internal to a unipath.  In that case we don't do anything.

          const int winf = 4;
          const int min_dist = 200;
          int M = 0, Mid = -1;
          for ( int j = 0; j < aligns_out.NAligns( ); j++ )
          {    if ( aligns_out.Count(j) >= M )
               {    M = aligns_out.Count(j);
                    Mid = j;    }    }
          Bool far = True;
          for ( int j = 0; j < aligns_out.NAligns( ); j++ )
               if ( j != Mid && aligns_out.Count(j) >= M/winf ) far = False;
          ostringstream hout;
          if ( far && Mid >= 0 )
          {    int u = aligns_out.U(Mid);
               int p1 = aligns_out.Upos( Mid, 0 );
               int p2 = aligns_out.Upos( Mid, M-1 ) + L;
               if ( p1 >= min_dist && p2 <= unibases[u].isize( ) - min_dist )
               {    hout << "\n0 Nothing to do.\n";
                    if ( !SKIP_SILENT ) reports[id] = hout.str( );
                    done[id] = True;
                    continue;    }    }

          // Kill inferior clusters.

          KmerAlignSet aly(ALIGNS);
          if (KILL_INFERIOR_CLUSTERS_NEW) 
          {    const double min_ratio_to_kill = 2.0;
               KillInferiorClustersNew( aly, unibases, min_ratio_to_kill );    }
          else KillInferiorClusters( aly );
          ALIGNS = aly.X( );

          // Print alignments.

          hout << "\n=========================================================="
               << "==========================\n";
          hout << "\nALIGNMENTS OF READ " << id << ", LENGTH = " 
               << longreads[id].size( ) << endl;
          int count = 0;

          // Make alignment graph.  This graph has one vertex for each alignment
          // group, and one edge for each adjacency.  The edge describes the
          // relative offset.

          const int verbosity = 0;
          digraphE<int> G;
          vec< pair< int, vec< pair<int,int> > > > ALIGNSX(ALIGNS);
          for ( int j = 0; j < ALIGNS.isize( ); j++ )
          {    for ( int l = 0; l < ALIGNS[j].second.isize( ); l++ )
               {    int rpos = ALIGNS[j].second[l].first;
                    int upos = ALIGNS[j].second[l].second;
                    ALIGNSX[j].second[l].first = upos;
                    ALIGNSX[j].second[l].second = rpos;    }    }
          MakeAlignmentGraph( ALIGNSX, unibases, verbosity, hout, G );
          if (STANDARD_ALIGNS)
          {    vec<int> to_delete;
               for ( int v = 0; v < G.N( ); v++ )
               {    for ( int j = 0; j < G.From(v).isize( ); j++ )
                    {    int e = G.EdgeObjectIndexByIndexFrom( v, j );
                         int u = ALIGNS[v].first;
                         int overlap = unibases[u].isize( ) - G.EdgeObject(e);
                         if ( overlap != K - 1 ) to_delete.push_back(e);    }    }
               G.DeleteEdges(to_delete);    }

          // Remove transitive edges.

          vec<int> edges_to_remove;
          for ( int x = 0; x < G.N( ); x++ )
          for ( int j1 = 0; j1 < G.From(x).isize( ); j1++ )
          {    int y = G.From(x)[j1];
               int e1 = G.EdgeObjectIndexByIndexFrom( x, j1 );
               if ( Member( edges_to_remove, e1 ) ) continue;
               int off1 = G.EdgeObject(e1);
               for ( int j2 = 0; j2 < G.From(y).isize( ); j2++ )
               {    int z = G.From(y)[j2];
                    int e2 = G.EdgeObjectIndexByIndexFrom( y, j2 );
                    if ( Member( edges_to_remove, e2 ) ) continue;
                    int off2 = G.EdgeObject(e2);
                    for ( int j3 = 0; j3 < G.From(x).isize( ); j3++ )
                    {    if ( G.From(x)[j3] == z )
                         {    int e3 = G.EdgeObjectIndexByIndexFrom( x, j3 );
                              int off3 = G.EdgeObject(e3);
                              if ( off3 == off1 + off2 )
                                   edges_to_remove.push_back(e3);    }    }    }    }
          UniqueSort(edges_to_remove);
          hout << "\nremoving " << edges_to_remove.size( ) 
               << " edges by transitivity" << endl;
          G.DeleteEdges(edges_to_remove);

          // Remove triple transitive edges.

          edges_to_remove.clear( );
          for ( int x = 0; x < G.N( ); x++ )
          for ( int j1 = 0; j1 < G.From(x).isize( ); j1++ )
          {    int y = G.From(x)[j1];
               int e1 = G.EdgeObjectIndexByIndexFrom( x, j1 );
               if ( Member( edges_to_remove, e1 ) ) continue;
               int off1 = G.EdgeObject(e1);
               for ( int j2 = 0; j2 < G.From(y).isize( ); j2++ )
               {    int z = G.From(y)[j2];
                    int e2 = G.EdgeObjectIndexByIndexFrom( y, j2 );
                    if ( Member( edges_to_remove, e2 ) ) continue;
                    int off2 = G.EdgeObject(e2);
                    for ( int j3 = 0; j3 < G.From(z).isize( ); j3++ )
                    {    int w = G.From(z)[j3];
                         int e3 = G.EdgeObjectIndexByIndexFrom( z, j3 );
                         if ( Member( edges_to_remove, e3 ) ) continue;
                         int off3 = G.EdgeObject(e3);
                         for ( int j4 = 0; j4 < G.From(x).isize( ); j4++ )
                         {    if ( G.From(x)[j4] == w )
                              {    int e4 = G.EdgeObjectIndexByIndexFrom( x, j4 );
                                   int off4 = G.EdgeObject(e4);
                                   if ( off4 == off1 + off2 + off3 )
                                        edges_to_remove.push_back(e4);    
                                        }    }    }    }    }
          UniqueSort(edges_to_remove);
          hout << "removing " << edges_to_remove.size( ) 
               << " edges by triple transitivity\n" << endl;
          G.DeleteEdges(edges_to_remove);

          // Clean graph.

          if (CLEAN_GRAPH)
          {    const int min_size = 10;
               vec<int> keep;
               vec<Bool> to_delete( G.N( ), False );
               for ( int v = 0; v < G.N( ); v++ )
               {    if ( G.From(v).empty( ) && G.To(v).empty( )
                         && ALIGNS[v].second.isize( ) < min_size )
                    {    to_delete[v] = True;    }
                    else keep.push_back(v);    }
               G = digraphE<int>( G, keep );
               EraseIf( ALIGNS, to_delete );
               EraseIf( ALIGNSX, to_delete );    }

          // Display alignment graph.

          for ( int i = 0; i < ALIGNS.isize( ); i++ )
          {    int u = ALIGNS[i].first;
               const vec< pair<int,int> >& x = ALIGNS[i].second;
               hout << "[" << count++ << "] u = " << u << "[l=" 
                    << unibases[u].size( ) << "]";
               if (ABBREVIATE_ALIGNMENTS)
               {    hout << " (" << x.front( ).first << "," << x.front( ).second 
                         << "," << x.front( ).second - x.front( ).first
                         << ") ...[" << x.size( ) << "]..."
                         << " (" << x.back( ).first << "," << x.back( ).second 
                         << "," << x.back( ).second - x.back( ).first << ")";    }
               else
               {    for ( int j = 0; j < x.isize( ); j++ )
                    {    hout << " (" << x[j].first << "," << x[j].second 
                              << "," << x[j].second - x[j].first << ")";    }    }
               hout << "\n";    }
          hout << "\n";
          vec<double> overlaps( G.EdgeObjectCount( ) );
          for ( int v = 0; v < G.N( ); v++ )
          {    for ( int j = 0; j < G.From(v).isize( ); j++ )
               {    int w = G.From(v)[j];
                    int e = G.EdgeObjectIndexByIndexFrom( v, j );
                    int u = ALIGNS[v].first;
                    overlaps[e] = unibases[u].isize( ) - G.EdgeObject(e);    }    }
          for ( int v = 0; v < G.N( ); v++ )
          {    for ( int j = 0; j < G.From(v).isize( ); j++ )
               {    int w = G.From(v)[j];
                    int e = G.EdgeObjectIndexByIndexFrom( v, j );
                    PRINT3_TO( hout, v, w, overlaps[e] );    }    }
          if (DOT1)
          {    Ofstream( out, "xxx.dot" );
               digraph(G).PrettyDOT( out, False, True );    
               // G.PrettyDOT( out, overlaps, False, True, True );
                    }
          if ( ALIGNS.size( ) <= 1 ) 
          {    hout << "\n1 Nothing to do.\n";
               if ( !SKIP_SILENT ) reports[id] = hout.str( );
               done[id] = True;
               continue;    }

          // Build alternate graph.

          digraphE<int> Galt(G);
          vec<int> Galt_U;
          for ( int j = 0; j < G.N( ); j++ )
               Galt_U.push_back( ALIGNS[j].first );
          for ( int x = 0; x < G.N( ); x++ )
          {    for ( int j = 0; j < G.From(x).isize( ); j++ )
               {    int y = G.From(x)[j], pos1;
                    int& over = Galt.EdgeObjectByIndexFromMutable( x, j );
                    for ( pos1 = 0; pos1 < ALIGNS[y].second.isize( ); pos1++ )
                    {    if ( ALIGNS[y].second[pos1].first 
                              == ALIGNS[x].second.back( ).first )
                         {    break;    }    }
                    if ( pos1 == ALIGNS[y].second.isize( ) ) over = 0;
                    else
                    {    int rpos1 = ALIGNS[y].second[pos1].first;
                         int upos1 = ALIGNS[y].second[pos1].second;
                         int rpos2 = ALIGNS[y].second.back( ).first;
                         int upos2 = ALIGNS[y].second.back( ).second;
                         if ( rpos1 == rpos2 ) over = 0;
                         else
                         {    int u2 = ALIGNS[y].first;
                              basevector z1( longreads[id], rpos1, rpos2 - rpos1 );
                              basevector z2( unibases[u2], upos1, upos2 - upos1 );
                              align a;
                              if ( z1.size( ) == 0 || z2.size( ) == 0 )
                                   over = 0;
                              else
                              {
                              over = SmithWatFreeSym( 
                                   z1, z2, a, True, True, 1, 1 );    }    }    }    }
                              }
          vec<int> verts( ALIGNS.size( ) );
          for ( int i = 0; i < ALIGNS.isize( ); i++ )
               verts[i] = ALIGNS[i].first;
          digraphVE<int,int> H( Galt, verts );
          Hall[id] = H;

          // Bail if graph is cyclic.  Very sloppy.

          if (CYCLIC_BAIL)
          {    if ( !G.Acyclic( ) )
               {    hout << "\nGraph is cyclic." << endl;
                    PRINT_TO( hout, id );
                    G.Clear( );    }    }

          // Find all paths.

          hout << "\n" << Date( ) << ": creating edgepaths" << endl;
          vec<int> sources, sinks;
          G.Sources(sources), G.Sinks(sinks);
          vec<int> to_left, to_right;
          G.ToLeft(to_left), G.ToRight(to_right);
          vec< vec<int> > epaths;
          vec<int> MATCHES;
          vec<basevector> BPATHS;
          vec< vec<int> > successors( sources.size( ) );
          for ( int i = 0; i < sources.isize( ); i++ )
          {    G.GetSuccessors1( sources[i], successors[i] );
               Sort( successors[i] );    }
          for ( int i1 = 0; i1 < sources.isize( ); i1++ )
          for ( int i2 = 0; i2 < sinks.isize( ); i2++ )
          {    if ( !BinMember( successors[i1], sinks[i2] ) ) continue;
               vec< vec<int> > epathsi;

               // Find shortest path.  This is an alternative to
               // finding all paths.  

               vec<int> path;
               Galt.ShortestPath( sources[i1], sinks[i2], path );
               vec<int> epath;
               for ( int j = 0; j < path.isize( ) - 1; j++ )
               {    int x = path[j], y = path[j+1];
                    int min_edge = 1000000000;
                    for ( int l = 0; l < Galt.From(x).isize( ); l++ )
                    {    if ( Galt.From(x)[l] == y )
                         {    min_edge = Min( min_edge,
                                   Galt.EdgeObjectByIndexFrom( x, l ) );    }    }
                    for ( int l = 0; l < G.From(x).isize( ); l++ )
                    {    if ( G.From(x)[l] == y )
                         {    if ( Galt.EdgeObjectByIndexFrom( x, l ) == min_edge )
                              {    epath.push_back( G.EdgeObjectIndexByIndexFrom( 
                                        x, l ) );
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
               int matches = M.size( );
               vec< vec<int> > epathsi_alt;
               if ( epath.size( ) > 0 )
               {    epathsi_alt.push_back(epath);
                    hout << "\nshortest path: ";
                    int len = 
                         unibases[ ALIGNS[ to_left[ epath[0] ] ].first ].size( );
                    int errs = 0;
                    for ( int r = 0; r < epath.isize( ); r++ )
                    {    int e = epath[r];
                         errs += Galt.EdgeObject(e);
                         int v = to_left[e], w = to_right[e];
                         int u = ALIGNS[v].first;
                         int o = unibases[u].isize( ) - G.EdgeObject(e);
                         if ( r == 0 ) hout << v << "<" << u << ">";;
                         hout << " --[" << o << "]--> " << w
                              << "<" << ALIGNS[w].first << ">";
                         len += unibases[ ALIGNS[w].first ].isize( ) - o;    }
                    hout << ", length = " << len
                         << ", read range = " 
                         << ALIGNS[ to_left[ epath.front( ) ] ].second.front( ).first
                         << "-"
                         << ALIGNS[ to_right[ epath.back( ) ] ].second.back( ).first
                              + L
                         << ", matches = " << matches 
                         << ", errs = " << errs << endl;    }

               Bool OK = True;
               if ( !USE_SHORTEST )
               {    OK = G.EdgePaths( sources[i1], sinks[i2], epathsi, -1,
                         heur.max_paths, heur.max_iterations );    }

               Bool found_epath = False;

               if ( !OK ) hout << "edge path computation exploded" << endl;
               if (USE_SHORTEST)
               {    epathsi.clear( ); 
                    if ( epath.size( ) > 0 ) epathsi.push_back(epath);    }

               {    vec<int> vpaths;
                    vec<basevector> bpaths;

                    // Remove nonsensical paths.

                    vec<Bool> to_remove( epathsi.size( ), False );
                    for ( int j = 0; j < epathsi.isize( ); j++ )
                    {    if ( !ConvertPathIntoBasesValid( L, G, epathsi[j],
                              ALIGNSX, unibases ) )
                         {    to_remove[j] = True;    }    }
                    EraseIf( epathsi, to_remove );

                    int verbosity = 0;
                    if ( epathsi.empty( ) ) continue;
                    hout << "have " << epathsi.size( ) << " epaths" << endl;
                    ConvertPathsIntoBases( L, longreads[id], G, epathsi, vpaths, 
                         ALIGNSX,
                         unibases, bpaths, verbosity, hout, run_dir, data_dir,
                         "Translate", True );    
                    vec<int> matches( bpaths.size( ) );
                    for ( int j = 0; j < bpaths.isize( ); j++ )
                    {    vec<int> M;
                         for ( int r = 0; r < epathsi[j].isize( ); r++ )
                         {    int e = epathsi[j][r];
                              int v = to_left[e], w = to_right[e];
                              for ( int m = 0; m < ALIGNS[v].second.isize( ); m++ )
                                   M.push_back( ALIGNS[v].second[m].first );
                              for ( int m = 0; m < ALIGNS[w].second.isize( ); m++ )
                                   M.push_back( ALIGNS[w].second[m].first );     }
                         UniqueSort(M);
                         matches[j] = M.size( );    }
                    int max_matches = Max(matches);
                    int mfudge = 0;
                    for ( int j = 0; j < bpaths.isize( ); j++ )
                    {    if ( matches[j] < max_matches - mfudge ) continue;
                         hout << "<" << j << "> ";
                         for ( int r = 0; r < epathsi[j].isize( ); r++ )
                         {    int e = epathsi[j][r];
                              int v = to_left[e], w = to_right[e];
                              int u = ALIGNS[v].first;
                              int o = unibases[u].isize( ) - G.EdgeObject(e);
                              if ( r == 0 ) hout << v << "<" << u << ">";;
                              hout << " --[" << o << "]--> " << w
                                   << "<" << ALIGNS[w].first << ">";    }
                         hout << ", matches = " << matches[j] << endl;    

                         if ( epath == epathsi[j] ) found_epath = True;
                         epaths.push_back( epathsi[j] );
                         MATCHES.push_back( matches[j] );
                         BPATHS.push_back( bpaths[j] );    }    }    
               if ( !found_epath ) hout << "DIFFERENT FROM SHORTEST PATH!\n";    }
          hout << "found " << epaths.size( ) << " epaths" << endl;

          // Call an alignment strong if for the range of read bases it matches,
          // there is no other alignment having more matches.  Also check epaths.

          vec<Bool> strong( ALIGNS.size( ), True );
          for ( int i1 = 0; i1 < ALIGNS.isize( ); i1++ )
          for ( int i2 = 0; i2 < ALIGNS.isize( ); i2++ )
          {    if ( i2 == i1 ) continue;
               const vec< pair<int,int> >& x1 = ALIGNS[i1].second;
               const vec< pair<int,int> >& x2 = ALIGNS[i2].second;
               int count1 = 0, count2 = x2.size( );
               int start2 = x2.front( ).first, stop2 = x2.back( ).first;
               for ( int j = 0; j < x1.isize( ); j++ )
                    if ( start2 <= x1[j].first && x1[j].first <= stop2 ) count1++;
               if ( count1 > count2 ) strong[i2] = False;    }
          for ( int j = 0; j < epaths.isize( ); j++ )
          {    const vec<int>& e = epaths[j];
               vec<int> v;
               v.push_back( to_left[ e[0] ] );
               for ( int l = 0; l < e.isize( ); l++ )
                    v.push_back( to_right[ e[l] ] );
               vec<int> r;
               for ( int l = 0; l < v.isize( ); l++ )
               {    for ( int m = 0; m < ALIGNS[ v[l] ].second.isize( ); m++ )
                         r.push_back( ALIGNS[ v[l] ].second[m].first );    }
               UniqueSort(r);
               for ( int i2 = 0; i2 < ALIGNS.isize( ); i2++ )
               {    const vec< pair<int,int> >& x2 = ALIGNS[i2].second;
                    int count1 = 0, count2 = x2.size( );
                    int start2 = x2.front( ).first, stop2 = x2.back( ).first;
                    for ( int z = 0; z < r.isize( ); z++ )
                         if ( start2 <= r[z] && r[z] <= stop2 ) count1++;
                    if ( count1 > count2 ) strong[i2] = False;    }    }

          // Look for good alignments that go up to the end of a unibase.

          Bool found_patch = False;
          const int max_dist_from_end = 100; // should probably depend on K
          const int min_span = 200;
          const int min_hits = 20;
          const int max_iterations = 20000;
          const double err_fudge = 60;
          const double max_rel_err = 0.3;
          vec<int> plefts, prights;
          for ( int i = 0; i < ALIGNS.isize( ); i++ )
          {    if ( !PATCHES_ONLY ) continue;
               if ( !strong[i] ) continue;
               int u = ALIGNS[i].first;
               const vec< pair<int,int> >& x = ALIGNS[i].second;
               if ( x.isize( ) < min_hits ) continue;
               if ( x.back( ).second + L - x.front( ).second < min_span ) continue;
               if ( x.front( ).second <= max_dist_from_end ) prights.push_back(i);
               if ( unibases[u].isize() - x.back().second - L <= max_dist_from_end )
                    plefts.push_back(i);    }
          for ( int j1 = 0; j1 < plefts.isize( ); j1++ )
          {    for ( int j2 = 0; j2 < prights.isize( ); j2++ )
               {    int i1 = plefts[j1], i2 = prights[j2];
                    int u1 = ALIGNS[i1].first, u2 = ALIGNS[i2].first;
                    if ( u2 == u1 || u2 == to_rc[u1] ) continue;
                    vec< pair<int,int> > x1 = ALIGNS[i1].second;
                    vec< pair<int,int> > x2 = ALIGNS[i2].second;
          
                    // Don't allow big jumps in offset near ends.

                    /*
                    {    const int max_jump = 10;
                         const int end_prox = 50;
                         int trunc1 = 0, trunc2 = 0;
                         int n1 = x1.size( ), N1 = unibases[u1].size( );
                         for ( int l = 1; l < x1.isize( ); l++ )
                         {    if ( N1 - x1[n1-l].second - L > end_prox ) break;
                              int offset1 = x1[n1-l-1].first - x1[n1-l-1].second;
                              int offset2 = x1[n1-l].first - x1[n1-l].second;
                              if ( Abs( offset1 - offset2 ) > max_jump )
                                   trunc1 = l;    }
                         x1.resize( x1.isize( ) - trunc1 );
                         for ( int l = 1; l < x2.isize( ); l++ )
                         {    if ( x2[l].second > end_prox ) break;
                              int offset1 = x2[l].first - x2[l].second;
                              int offset2 = x2[l-1].first - x2[l-1].second;
                              if ( Abs( offset1 - offset2 ) > max_jump )
                                   trunc2 = l;    }
                         vec< pair<int,int> > x2_new;
                         x2_new.SetToSubOf( x2, trunc2, x2.isize( ) - trunc2 );
                         x2 = x2_new;    }
                    */

                    int read_gap = x2.front( ).first - x1.back( ).first - L;
                    int left_overhang 
                         = unibases[u1].isize( ) - x1.back( ).second - L;
                    int right_overhang = x2.front( ).second;
                    int expected_gap = read_gap - left_overhang - right_overhang;

                    // Special handling for overlaps.

                    if ( expected_gap <= -L )
                    {
                         // First find all overlaps between u1 and u2.

                         vec<int> o;
                         for ( int m = L; m < K; m++ )
                         {    if ( unibases[u1].Overlap( unibases[u2], m ) )
                                   o.push_back(m);    }
                         if ( o.nonempty( ) )
                         {

                         // Find 'trusted matches' by moving back from the ends.

                         pair<int,int> left, right;
                         const int advance = 100;
                         for ( int l = x1.isize( ) - 1; l >= 0; l-- )
                         {    if ( l == 0 || 
                                   unibases[u1].isize( ) - x1[l].second >= advance )
                              {    left = x1[l];
                                   break;    }    }
                         for ( int l = 0; l < x2.isize( ); l++ )
                         {    if ( l == x2.isize( ) - 1
                                   || x2[l].second + L >= advance )
                              {    right = x2[l];
                                   break;    }    }
                         if ( right.first + L - left.first <= 0 ) continue;

                         // For each overlap, Smith-Waterman to the joined unibases.

                         vec<align> sw( o.size( ) );
                         vec<int> errs( o.size( ) );
                         vec< vec<ho_interval> > perf1( o.size( ) );
                         vec< vec<ho_interval> > perf2( o.size( ) );
                         for ( int l = 0; l < o.isize( ); l++ )
                         {    basevector u = basevector( unibases[u1], left.second,
                                   unibases[u1].isize( ) - left.second - o[l] );
                              u = Cat( u, basevector( unibases[u2], 0,
                                   right.second + L ) );
                              basevector r( 
                                   x, left.first, right.first + L - left.first );
                              errs[l] = SmithWatFreeSym( 
                                   r, u, sw[l], True, True, 1, 1 );
                              sw[l].PerfectIntervals1( r, u, perf1[l] );
                              sw[l].PerfectIntervals2( r, u, perf2[l] );    }
     
                         // Find the best overlap.

                         int me = Min(errs), bo;
                         for ( bo = 0; bo < o.isize( ); bo++ )
                              if ( errs[bo] == me ) break;

                         // Use the best overlap to replace x1 and x2.

                         vec< pair<int,int> > x1_old(x1), x2_old(x2); // XXXXXXXXXXX
                         vec< pair<int,int> > x1_new, x2_new;
                         for ( int m = 0; m < perf1[bo].isize( ); m++ )
                         {    const ho_interval& h1 = perf1[bo][m];
                              const ho_interval& h2 = perf2[bo][m];
                              for ( int l = 0; l <= h1.Length( ) - L; l++ )
                              {    int p1 = h1.Start( ) + l + left.first; 
                                   int p2 = h2.Start( ) + l + left.second;
                                   if ( p2 + L <= unibases[u1].isize( ) )
                                        x1_new.push( p1, p2 );
                                   if ( p2 >= unibases[u1].isize( ) - o[bo] )
                                   {    x2_new.push( p1, p2 - ( unibases[u1].isize( )
                                             - o[bo] ) );    }    }    }
                         if ( x1_new.nonempty( ) && x2_new.nonempty( ) )
                         {    x1 = x1_new;
                              x2 = x2_new;    }    }    }

                    Bool bridged = False;
                    for ( int pass = 1; pass <= 3; pass++ )
                    {    if (bridged) break;
                         int source, target;
                         if ( pass == 1 )
                         {    source = u1, target = u2;    }
                         else if ( pass == 2 )
                         {    source = u1, target = to_rc[u2];    }
                         else
                         {    source = to_rc[u2], target = u1;    }
                         double max_gap = ( 1.0 + max_rel_err ) 
                              * double(expected_gap) + err_fudge;
                         double min_gap;
                         if ( pass == 1 )
                         {    min_gap = ( 1.0 - max_rel_err ) 
                                   * double(expected_gap) - err_fudge;    }
                         else min_gap = - K - 1;
                         vec< pair<int,int> > partials;
                         partials.push( source, - K + 1 );
                         int iterations = 0;
                         while( partials.nonempty( ) )
                         {    if ( ++iterations > max_iterations ) 
                              {    break;    }
                              int v1 = partials.back( ).first;
                              int gap1 = partials.back( ).second;
                              partials.pop_back( );
                              if ( gap1 > max_gap ) continue;
                              if ( v1 == target && gap1 >= min_gap )
                              {    bridged = True;
                                   break;    }
                              if( v1 == target ) 
                                   gap1 += unibases[target].isize( ) - K + 1;
                              for ( int j = 0; j < nexts[v1].isize( ); j++ )
                              {    int v2 = nexts[v1][j];
                                   int gap2 = gap1;
                                   if ( v2 != target ) 
                                        gap2 += unibases[v2].isize( ) - K + 1;
                                   partials.push( v2, gap2 );    }    }    }
                    int n1 = unibases[u1].size( ), n2 = unibases[u2].size( );
                    if ( Max(n1,n2) >= MIN_PATCH1 && Min(n1,n2) >= MIN_PATCH2 )
                    {
                         #pragma omp critical
                         {    found_patch = True;
                              int best_len = 0;

                              // tpos1 is the beginning of the last L-mer on u1
                              // tpos2 is the end of the first L-mer on u2
                              // rpos1 and rpos2 are the corresponding 
                              // positions on r

                              int tpos1 = x1.back( ).second;
                              int tpos2 = x2.front( ).second + L;
                              int rpos1 = x1.back( ).first;
                              int rpos2 = x2.front( ).first + L;

                              // For now, if r would be of negative length,
                              // we move the flanking L-mers away.  This is not
                              // necessarily the right thing to do.

                              int k1 = x1.isize( ) - 1;
                              int k2 = 0;
                              Bool fail = False;
                              while( rpos2 - rpos1 < 0 )
                              {    k1--;
                                   k2++;
                                   if ( k1 < 0 || k2 >= x2.isize( ) )
                                   {    fail = True;
                                        break;    }
                                   tpos1 = x1[k1].second;
                                   tpos2 = x2[k2].second + L;
                                   rpos1 = x1[k1].first;
                                   rpos2 = x2[k2].first + L;
                                   left_overhang 
                                        = unibases[u1].isize( ) - x1[k1].second - L;
                                   right_overhang = x2[k2].second;
                                   read_gap = x2[k2].first - x1[k1].first - L;
                                   expected_gap = read_gap 
                                        - left_overhang - right_overhang;    }
                              if (fail) hout << "can't avoid negative r\n";
                              else
                              {

                              // Now define r.

                              basevector r( longreads[id], rpos1, rpos2 - rpos1 );

                              // Save the patch, and then its reverse complement.

                              patchers.push( u1, u2, r, 2*id, tpos1, tpos2,
                                   best_len, expected_gap );
                              hout << "possible gap " << u1 << " --(" 
                                   << expected_gap << ")--> " << u2 << ", from [" 
                                   << i1 << "] --> [" << i2 << "], fw\n";
                              r.ReverseComplement( );
                              int tpos1x = unibases[u2].isize( ) - tpos2;
                              int tpos2x = unibases[u1].isize( ) - tpos1;
                              tpos1 = tpos1x, tpos2 = tpos2x;
                              patchers.push( to_rc[u2], to_rc[u1], r, 2*id + 1,
                                   tpos1, tpos2, best_len, expected_gap );    
                              hout << "possible gap " << to_rc[u2] << " --(" 
                                   << expected_gap << ")--> " << to_rc[u1]
                                   << ", from [" << i1 << "] --> [" << i2 
                                   << "], rc\n";    }    }    }    }    }
          if (PATCHES_ONLY)
          {    reports[id] = hout.str( );
               continue;    }

          // Find vpaths.

          vec<int> vpaths;
          vec<Bool> used( G.N( ), False );
          for ( int j = 0; j < epaths.isize( ); j++ )
          {    for ( int l = 0; l < epaths[j].isize( ); l++ )
               {    int e = epaths[j][l];
                    used[ to_left[e] ] = True, used[ to_right[e] ] = True;    }    }
          for ( int v = 0; v < G.N( ); v++ )
               if ( !used[v] ) vpaths.push_back(v);
          hout << "found " << vpaths.size( ) << " vpaths" << endl;

          // Print best paths.

          ReverseSortSync( MATCHES, epaths, BPATHS );
          hout << "\nBEST MATCHES FOR READ " << id << ":\n\n";
          vec<basevector> B;
          vec<int> E;
          for ( int j = 0; j < epaths.isize( ); j++ )
          {    if ( MATCHES[j] < MATCHES[0] ) break;
               B.push_back( BPATHS[j] );
               E.push_back(j);    }
          UniqueSortSync( B, E );
          for ( int j = 0; j < epaths.isize( ); j++ )
          {    if ( MATCHES[j] < MATCHES[0] ) break;
               hout << "<" << j << "> ";
               for ( int r = 0; r < epaths[j].isize( ); r++ )
               {    int e = epaths[j][r];
                    int v = to_left[e], w = to_right[e];
                    int u = ALIGNS[v].first;
                    int o = unibases[u].isize( ) - G.EdgeObject(e);
                    if ( r == 0 ) hout << v << "<" << u << ">";;
                    hout << " --[" << o << "]--> " << w
                         << "<" << ALIGNS[w].first << ">";    }
               hout << ", matches = " << MATCHES[j] 
                    << ", bases = " << BPATHS[j].size( ) << endl;    
               if (QLT1)
               {    {    Ofstream( out, "slobber.fasta" );
                         BPATHS[j].Print( out, j );    }
                    hout << AllOfOutput1( "QueryLookupTable K=12 MM=12 MC=0.15 "
                         "SEQS=slobber.fasta L=" + data_dir + "/genome.lookup "
                         "KB=1 VISUAL=True NH=True QUIET=True "
                         "QUERY_NAMING=from_record" );    }    }

          // Check the vpaths.

          Bool reject = False;
          for ( int j = 0; j < vpaths.isize( ); j++ )
          {    int v = vpaths[j];
               int matches = ALIGNS[v].second.size( );
               int u = ALIGNS[v].first;
               if ( MATCHES.empty( ) || matches >= MATCHES[0] )
               {    hout << "\nalignment of unipath " << u << " has "
                         << matches << " matches\n";    
                    reject = True;    }    }

          // Save sequences.

          if ( !reject )
          {    
               #pragma omp critical
               {    // all.append(B);    

                    if ( B.size( ) > 0 ) 
                    {    all.push_back_reserve( B[0] );
                         COMPUTED[id] = True;
                         hout << "\nSAVING ";
                         int j = E[0];
                         vec<int> U, overlap;
                         int nbases = 0;
                         vec<int> Uid;
                         for ( int r = 0; r < epaths[j].isize( ); r++ )
                         {    int e = epaths[j][r];
                              int v = to_left[e], w = to_right[e];
                              int u = ALIGNS[v].first;
                              U.push_back(u);
                              Uid.push_back(v);
                              if ( r == epaths[j].isize( ) - 1 )
                              {    U.push_back( ALIGNS[w].first );
                                   Uid.push_back(w);    }
                              int o = unibases[u].isize( ) - G.EdgeObject(e);
                              nbases += unibases[u].isize( ) - o;
                              if ( r == epaths[j].isize( ) - 1 )
                                   nbases += unibases[ ALIGNS[w].first ].size( );
                              overlap.push_back(o);
                              if ( r == 0 ) hout << u;
                              hout << " --[" << o << "]--> " << ALIGNS[w].first;    }
                         hout << ", bases = " << nbases << "\n";    
                         UNISEQ[id] = uniseq( U, overlap );    
                         UNISEQ_ID[id] = Uid;    }    }    }
          else
          {    hout << "\n3 Nothing to do.\n";
               if ( !SKIP_SILENT || found_patch ) reports[id] = hout.str( );
               done[id] = True;
               continue;    }

          // Print.

          reports[id] = hout.str( );
          done[id] = True;
          if (VERBOSE1)
          {    
               #pragma omp critical
               {    if ( done[done_ptr] )
                    {    while(1)
                         {    cout << reports[done_ptr];
                              done_ptr++;
                              if ( done_ptr == done.isize( ) || !done[done_ptr] ) 
                                   break;    }    }    }    }    }

     // Finish printing reports.

     for ( int id = done_ptr; id < done.isize( ); id++ )
          if ( VERBOSE1 || BinMember( reads_to_trace, id ) ) cout << reports[id];
     cout << "\n";
     DPRINT( all.size( ) );
     omp_set_num_threads(NUM_THREADS);    }

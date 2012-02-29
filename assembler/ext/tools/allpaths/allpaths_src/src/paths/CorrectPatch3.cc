///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// CorrectPatch3 uses corrected fragment read pair data to correct long read 
// patches.
//
// Known issues:
//
// 1. Because we place read pairs by finding the best placement, it follows that
// if the genome has a sufficiently long perfect tandem duplication, and one copy is 
// sequenced perfectly, whereas the second is sequenced imperfectly, the second copy
// will not be corrected.

#include "Basevector.h"
#include "CoreTools.h"
#include "MergeReads2.h"
#include "PackAlign.h"
#include "PairsManager.h"
#include "PrintAlignment.h"
#include "Set.h"
#include "VecUtilities.h"
#include "graph/DigraphTemplate.h"
#include "kmers/KmerRecord.h"
#include "math/Functions.h"
#include "math/HoInterval.h"
#include "pairwise_aligners/SmithWatBandedA.h"
#include "paths/AssemblyEdit.h"
#include "paths/CorrectPatch3.h"
#include "paths/KmerBaseBroker.h"
#include "paths/LongReadTools.h"
#include "paths/GetNexts.h"
#include "paths/ReadsToPathsCoreX.h"
#include "paths/Unipath.h"

// Trim.  Chew back sources and sinks that are not the original source/sink.
// Also delete spurious components.

int Trim( const vec< vec< triple<int,int,Bool> > >& sux, 
     digraphV< triple<int,int,Bool> >& M )
{    int total_trims = 0;
     while(1)
     {    int trims = 0;
          for ( int v = 0; v < M.N( ); v++ )
          {    
               // We use a more general notion of source and sink that takes
               // account of self loops.

               if ( ( M.To(v).empty( ) 
                    || M.To(v).solo( ) && BinSubset( M.To(v), M.From(v) ) )
                    && !( M.Vert(v) == sux[0].front( ) && M.Vert(v).third ) )
               {    M.DeleteVertex(v);
                    trims++;
                    break;    }
               if ( ( M.From(v).empty( ) 
                    || M.From(v).solo( ) && BinSubset( M.From(v), M.To(v) ) )
                    && !( M.Vert(v) == sux[0].back( ) && M.Vert(v).third ) )
               {    M.DeleteVertex(v);
                    trims++;
                    break;    }    }
          total_trims += trims;
          if ( trims == 0 ) break;    }
     int start = -1, stop = -1;
     for ( int j = 0; j < M.N( ); j++ )
     {    if ( M.Vert(j) == sux[0].front( ) ) 
          {    start = j;
               break;    }    }
     for ( int j = M.N( ) - 1; j >= 0; j-- )
     {    if ( M.Vert(j) == sux[0].back( ) ) 
          {    stop = j;
               break;    }    }
     vec< vec<int> > comp;
     vec<int> to_delete;
     M.Components(comp);
     for ( int i = 0; i < comp.isize( ); i++ )
     {    if ( !BinMember( comp[i], start ) || !BinMember( comp[i], stop ) )
               to_delete.append( comp[i] );    }
     Sort(to_delete);
     M.DeleteVertices(to_delete);
     return total_trims + to_delete.isize( );    }

// Compute coverage.  Note that the way we compute coverage awards points to
// pairs even if foreign material has been removed from them.  One can imagine
// a smarter algorithm that would kill these reads entirely.

void GetPairCoverage( const vec< vec< triple<int,int,Bool> > >& sux, 
     const digraphV< triple<int,int,Bool> >& M, vec<double>& covx )
{    covx.resize_and_set( M.N( ), 0 );
     for ( int id = 0; id < sux.isize( ); id++ )
     {    for ( int j = 0; j < sux[id].isize( ); j++ )
          {    for ( int l = 0; l < M.N( ); l++ )
                    if ( M.Vert(l) == sux[id][j] ) covx[l]++;    }    }    }

// Find the placements of one read.  A placement is presented as a sequence of
// integers (x1,...,xn) equal to the length of the read (measured by its number of
// unipaths).  If xi >= 0, it represents a vertex id on the graph.  If xi < 0,
// then -xi-1 is a unipath id.  We require that the read is placed starting at its
// beginning, so we always have x1 >= 0.

void FindPlacements( const int id, const vec< vec<int> >& sux, const vec<int>& dx, 
     digraphV< triple<int,int,Bool> >& M, vec< vec<int> >& placements )
{
     placements.clear( );
     const vec<int>& x = sux[id];
     if ( x.empty( ) ) return;
     const int max_dist_delta = 5;
     for ( int v = 0; v < M.N( ); v++ )
     {    if ( M.Vert(v).first != x[0] ) continue;
          vec<int> place, distx;
          place.push_back(v), distx.push_back( M.Vert(v).second );
          for ( int z = 1; z < x.isize( ); z++ )
          {    int w = x[z];
          
               // There are three possibilities for what we can do with w.
               // 1. Check for direct placement.
     
               int vc = place.back( );

               if ( vc >= 0 )
               {    for ( int i = 0; i < M.From(vc).isize( ); i++ )
                    {    int vn = M.From(vc)[i];
                         if ( M.Vert(vn).first == w )
                         {    place.push_back(vn); 
                              distx.push_back( M.Vert(vn).second );
                              break;    }    }
                    if ( place.isize( ) == z + 1 ) continue;    }

               // Nope.  Look off into space, adding a new edge.

               int d = distx.back( ) + dx[ x[z-1] ];
               for ( int t = 0; t < M.N( ); t++ )
               {    if ( M.Vert(t).first != w ) continue;

                    // 2. We accept a given graph vertex if the distance difference
                    // s at most a fixed constant.  Note that this is a messy point.

                    if ( Abs( d - M.Vert(t).second ) > max_dist_delta ) continue;
                    place.push_back(t), distx.push_back( M.Vert(t).second );
                    break;    }
               if ( place.isize( ) == z + 1 ) continue;

               // 3. Now we have to stick the vertex out in space.

               place.push_back(-w-1), distx.push_back(d);    }

          // Record possible placement.

          placements.push_back(place);    }    }

// InsertRead.  Given a read (defined by x), insert it into the graph at v.

void InsertRead( const vec<int>& x, const int v, const vec<int>& dx,
     digraphV< triple<int,int,Bool> >& M )
{
     const int max_dist_delta = 5; // note multiple definitions
     vec<int> place, distx;
     place.push_back(v), distx.push_back( M.Vert(v).second );
     for ( int z = 1; z < x.isize( ); z++ )
     {    int w = x[z], vc = place.back( );
          for ( int i = 0; i < M.From(vc).isize( ); i++ )
          {    int vn = M.From(vc)[i];
               if ( M.Vert(vn).first == w )
               {    place.push_back(vn); 
                    distx.push_back( M.Vert(vn).second );
                    break;    }    }
          if ( place.isize( ) == z + 1 ) continue;
          int d = distx.back( ) + dx[ x[z-1] ];
          for ( int t = 0; t < M.N( ); t++ )
          {    if ( M.Vert(t).first != w ) continue;
               if ( Abs( d - M.Vert(t).second ) > max_dist_delta ) continue;

               // Add an edge to the graph.

               M.AddEdge( place.back( ), t );
               place.push_back(t); 
               distx.push_back( M.Vert(t).second );
               break;    }
          if ( place.isize( ) == z + 1 ) continue;

          // Add a vertex and an edge to the graph.

          M.AddVertex( make_triple( w, d, False ) );
          M.AddEdge( vc, M.N( ) - 1 );
          place.push_back( M.N( ) - 1 ); 
          distx.push_back(d);    }    }

// Place one read, requiring that its first position is in alignment.  Return True
// if the read is placed at least once.

Bool PlaceRead( const int id, const vec< vec<int> >& sux, const vec<int>& dx, 
     digraphV< triple<int,int,Bool> >& M )
{
     // First find all the placements.

     vec< vec<int> > placements;
     FindPlacements( id, sux, dx, M, placements );

     // Now score the placements, as measured by the number of vertices and edges 
     // added (smaller is better).  

     vec<int> where;
     vec< pair<int,int> > score;
     for ( int i = 0; i < placements.isize( ); i++ )
     {    const vec<int>& p = placements[i];
          where.push_back( p[0] );
          int edges_added = 0, verts_added = 0;
          for ( int j = 0; j < p.isize( ); j++ )
          {    if ( p[j] >= 0 && j > 0 && p[j-1] >= 0 )
               {    if ( !Member( M.From( p[j-1] ), p[j] ) ) edges_added++;    }
               if ( p[j] < 0 ) 
               {    edges_added++;
                    verts_added++;    }    }
          score.push( verts_added, edges_added );    }

     // Select the best placements and install them.

     if ( where.empty( ) ) return False;
     const vec<int>& x = sux[id];
     SortSync( score, where );
     int j;
     for ( j = 1; j < where.isize( ); j++ )
          if ( score[j] != score[0] ) break;
     for ( int k = 0; k < j; k++ )
     {    
          // Go back through the same process as above to place the read, but this
          // time edit the graph.  It is conceivable that the earlier placements of
          // the same read will affect the outcome.

          InsertRead( x, where[k], dx, M );    }

     return placements.nonempty( );    }

// Iteratively add reads, requiring that each can be placed starting at its 
// beginning.

void PlaceReads( const vec< vec<int> >& sux, const vec<int>& dx, 
     digraphV< triple<int,int,Bool> >& M )
{    int nseq = sux.size( );
     vec<Bool> placed( nseq, False );
     placed[0] = True;
     while(1)
     {    int adds = 0;
          for ( int id = 1; id < nseq; id++ )
          {    if ( placed[id] ) continue;
               if ( PlaceRead( id, sux, dx, M ) ) 
               {    placed[id] = True;    
                    adds++;    }    }    
          if ( adds == 0 ) break;    }    }

// RemoveBubbles.  First look for bubbles in which one branch has low coverage and
// the other branch has high coverage, and delete the low coverage branch.  Then
// look for 'groups' in the graph (defined below), try to find a 'reference patch'
// across each, and if we can find it, delete the rest of the group.

void RemoveBubbles( const vec< vec< triple<int,int,Bool> > >& sux, 
     digraphV< triple<int,int,Bool> >& M,
     const Bool verbose, ostream& out )
{
     // Iteratively remove bubbles.

     vec<double> covx;
     const double min_cov = 5.0;
     const double cov_mult = 5.0;
     while(1)
     {
          // Get coverage.

          GetPairCoverage( sux, M, covx );

          // Look for bubbles for which one branch has low minimum coverage and the 
          // other branch has high minimum coverage.  Delete the weak branch.

          vec<int> verts_to_delete;
          for ( int v = 0; v < M.N( ); v++ )
          {    if ( M.From(v).size( ) != 2 ) continue;
               vec< vec<int> > b(2);
               vec<Bool> clean( 2, True );
               Bool bad = False;
               for ( int l = 0; l < 2; l++ )
               {    b[l].push_back( v, M.From(v)[l] );
                    while(1)
                    {    if ( M.From( b[l].back( ) ).empty( ) ) break;

                         // Prevent infinite looping.  Not sure how this can happen.

                         int x = M.From( b[l].back( ) )[0];
                         if ( Member( b[l], x ) ) break;

                         b[l].push_back(x);
                         if ( !M.To( b[l].back( ) ).solo( ) ) break;    }
                    if ( M.To( b[l].back( ) ).size( ) != 2 ) bad = True;
                    for ( int j = 1; j < b[l].isize( ) - 1; j++ )
                    {    if ( M.From( b[l][j] ).size( ) 
                              + M.To( b[l][j] ).size( ) != 2 )
                         {    clean[l] = False;    }    }    }
               if ( b[0].back( ) != b[1].back( ) ) bad = True;
               if (bad) continue;
               if (verbose)
               {    out << "found bubble from " << b[0].front( ) << " = "
                         << M.Vert( b[0].front( ) ).first << " to " << b[0].back( ) 
                         << " = " << M.Vert( b[0].back( ) ).first << endl;    }
               vec<double> mc(2);
               int infinity = 1000000000;
               mc[0] = infinity, mc[1] = infinity;
               for ( int l = 0; l < 2; l++ )
               {    for ( int j = 0; j < b[l].isize( ); j++ )
                         mc[l] = Min( mc[l], covx[ b[l][j] ] );    }
               if ( mc[1] <= 2.0 ) 
               {    swap( mc[0], mc[1] );
                    swap( b[0], b[1] );
                    swap( clean[0], clean[1] );    }
               if ( mc[0] > 2.0 ) continue;
               if ( mc[1] < min_cov ) continue;
               if ( mc[0] && mc[1] < cov_mult * mc[0] ) continue;
               if ( !clean[0] ) continue;
               for ( int j = 1; j < b[0].isize( ) - 1; j++ )
                    verts_to_delete.push_back( b[0][j] );    }
          if (verbose)
               out << "deleting " << verts_to_delete.size( ) << " vertices" << endl;
          UniqueSort(verts_to_delete);
          M.DeleteVertices(verts_to_delete);

          // Get coverage.

          GetPairCoverage( sux, M, covx );

          // Try a different bubble removing method.  First we define the notion
          // of a 'group' in a graph.  It is all vertices between two vertices v
          // and w, subject to the following conditions on v and w:
          // - |from(v)| > 1
          // - consider all x in succ(from(v)) such that
          //   for all b in betw(v,x) - {v,x}, from(b) and to(b) are in betw(v,x)
          // - w is the unique source of these x.

          if (verbose) out << "\n" << Date( ) << ": finding groups\n";
          vec< pair<int,int> > groups;
          for ( int v = 0; v < M.N( ); v++ )
          {    if ( M.From(v).size( ) < 2 ) continue;
               vec<int> succ_v;
               M.GetSuccessors1( v, succ_v );
               vec<int> succ_from_v;
               vec< vec<int> > s;
               for ( int j = 0; j < M.From(v).isize( ); j++ )
               {    vec<int> n;
                    M.GetSuccessors1( M.From(v)[j], n );
                    s.push_back(n);    }
               Intersection( s, succ_from_v );
               vec<int> good_x;
               for ( int j = 0; j < succ_from_v.isize( ); j++ )
               {    int x = succ_from_v[j];
                    vec<int> pred_x;
                    M.GetPredecessors1( x, pred_x );
                    vec<int> betw = Intersection( succ_v, pred_x );
                    Bool x_good = True;
                    for ( int l = 0; l < betw.isize( ); l++ )
                    {    int b = betw[l];
                         if ( b == v || b == x ) continue;
                         if ( !BinSubset( M.From(b), betw )
                              || !BinSubset( M.To(b), betw ) )
                         {    x_good = False;
                              break;    }    }
                    if (x_good) good_x.push_back(x);    }
               vec<int> sources0, sources;
               M.SubgraphSources( good_x, sources0 );
               for ( int j = 0; j < sources0.isize( ); j++ )
               {    int s = sources0[j];
                    vec<int> ss, ssp;
                    for ( int k = 0; k < sources0.isize( ); k++ )
                         if ( k != j ) ss.push_back( sources0[k] );
                    M.GetSuccessors( ss, ssp );
                    if ( !BinMember( ssp, s ) ) sources.push_back(s);    }
               if ( verbose && sources.solo( ) ) PRINT2_TO( out, v, sources[0] );
               if ( sources.solo( ) ) groups.push( v, sources[0] );    }

          // Study groups.

          if (verbose) out << Date( ) << ": studying groups" << endl;
          vec<int> to_kill;
          for ( int gi = 0; gi < groups.isize( ); gi++ )
          {    int v = groups[gi].first, w = groups[gi].second;
               vec<int> from_v, to_w, betw;
               M.GetSuccessors1( v, from_v );
               M.GetPredecessors1( w, to_w );
               betw = Intersection( from_v, to_w );

               // Find all paths through the group.  Dangerous.

               if (verbose)
               {    out << Date( ) << ": getting all paths for group " << gi 
                         << " of " << groups.size( ) << ", from " << v << " to " 
                         << w << endl;    }
               vec< vec<int> > paths;
               M.AllPaths( v, w, paths );
               if ( paths.size( ) == 1 ) continue; // not sure how this can happen!

               // See if there is a unique path that best follows the patch.

               if (verbose) out << Date( ) << ": looking for refpath" << endl;
               vec<int> refpath;
               vec< vec<int> > refhits( paths.size( ) );
               vec<int> nrefhits( paths.size( ), 0 ); 
               vec<int> pid( paths.size( ), vec<int>::IDENTITY );
               for ( int j = 0; j < paths.isize( ); j++ )
               {    const vec<int>& p = paths[j];
                    for ( int l = 1; l < p.isize( ) - 1; l++ )
                    {    if ( M.Vert( p[l] ).third ) 
                         {    refhits[j].push_back( p[l] );
                              nrefhits[j]++;    }    }
                    UniqueSort( refhits[j] );    }
               ReverseSortSync( nrefhits, refhits, pid );
               if ( nrefhits[0] > nrefhits[1] ) 
               {    Bool OK = True;
                    for ( int l = 1; l < paths.isize( ); l++ )
                    {    if ( nrefhits[l] == 0 ) break;
                         if ( !BinSubset( refhits[l], refhits[0] ) ) 
                         {    OK = False;
                              break;    }    }
                    if (OK) refpath = paths[ pid[0] ];    }
               if ( verbose && refpath.size( ) > 0 )
               {    out << "refpath =";
                    for ( int j = 0; j < refpath.isize( ); j++ )
                         out << " " << refpath[j];
                    out << "\n";    }

               // Find the best coverage from v to w.

               double best = 0.0;
               const double best_div = 3.0;
               for ( int j = 0; j < paths.isize( ); j++ )
               {    double cov = 1000000000.0;
                    for ( int l = 0; l < paths[j].isize( ); l++ )
                         cov = Min( cov, covx[ paths[j][l] ] );
                    best = Max( best, cov );    }
               if (verbose) PRINT3_TO( out, v, w, best );

               // If the refpath exists, and has good coverage, and beats all our
               // paths, delete them.

               const double min_cov_for_refpath = 1.0;
               Bool refpath_wins = False;
               if ( refpath.nonempty( ) )
               {    double min_cov = covx[ refpath[0] ];
                    for ( int j = 0; j < refpath.isize( ); j++ )
                         min_cov = Min( min_cov, covx[ refpath[j] ] );
                    if ( min_cov >= Max( min_cov_for_refpath, best/best_div ) )
                    {    refpath_wins = True;
                         if (verbose) out << "refpath wins" << endl;
                         for ( int j = 0; j < betw.isize( ); j++ )
                         {    int x = betw[j];
                              if ( !Member( refpath, x ) )
                                  to_kill.push_back(x);    }    }    }

               // Find vertices in the group that are highly undercovered.

               if ( !refpath_wins )
               {    for ( int j = 0; j < betw.isize( ); j++ )
                    {    int x = betw[j];
                         if ( best >= min_cov && best >= cov_mult * covx[x] )
                         {    if (verbose) 
                                   out << "should kill vertex " << x << endl;
                              to_kill.push_back(x);    }    }    }    }

          // Kill vertices.

          UniqueSort(to_kill);
          M.DeleteVertices(to_kill);
          int trims = Trim( sux, M );

          // Are we done?

          if ( verts_to_delete.empty( ) && to_kill.empty( ) && trims == 0 ) 
               break;    }    }

// PickBranches.  At a branch vertex, try to pick the branch that is right.

void PickBranches( const vec< vec< triple<int,int,Bool> > >& sux, 
     digraphV< triple<int,int,Bool> >& M )
{    const double min_branch_ratio = 5.0;
     vec<double> covx;
     GetPairCoverage( sux, M, covx );
     while(1)
     {    vec<int> verts_to_go;
          for ( int v = 0; v < M.N( ); v++ )
          {    if ( M.From(v).size( ) > 1 )
               {    vec<int> ref;
                    for ( int j = 0; j < M.From(v).isize( ); j++ )
                         if ( M.Vert( M.From(v)[j] ).third ) ref.push_back(j);
                    if ( ref.solo( ) )
                    {    int r = ref[0];
                         vec<int> bads;
                         for ( int j = 0; j < M.From(v).isize( ); j++ )
                         {    if ( j == r ) continue;
                              if ( covx[ M.From(v)[r] ] 
                                   >= min_branch_ratio * covx[ M.From(v)[j] ] )
                              {    bads.push_back(j);    }    }
                         for ( int z = 0; z < bads.isize( ); z++ )
                         {    verts_to_go.push_back( 
                                   M.From(v)[ bads[z] ] );    }    }    }
               if ( M.To(v).size( ) > 1 )
               {    vec<int> ref;
                    for ( int j = 0; j < M.To(v).isize( ); j++ )
                         if ( M.Vert( M.To(v)[j] ).third ) ref.push_back(j);
                    if ( ref.solo( ) )
                    {    int r = ref[0];
                         vec<int> bads;
                         for ( int j = 0; j < M.To(v).isize( ); j++ )
                         {    if ( j == r ) continue;
                              if ( covx[ M.To(v)[r] ] 
                                   >= min_branch_ratio * covx[ M.To(v)[j] ] )
                              {    bads.push_back(j);    }    }
                         for ( int z = 0; z < bads.isize( ); z++ )
                         {    verts_to_go.push_back( 
                                   M.To(v)[ bads[z] ] );    }    }    }    }
          if ( verts_to_go.empty( ) ) break;
          UniqueSort(verts_to_go);
          M.DeleteVertices(verts_to_go);
          Trim( sux, M );
          GetPairCoverage( sux, M, covx );    }    }

Bool CorrectPatch3( const basevector& LEFT, const basevector& RIGHT,
     const vecbasevector& fbases, const vecqualvector& fquals, 
     const PairsManager& fpairs, const vec< kmer<20> >& fheads, 
     const vec<int64_t>& fids, assembly_edit& e, ostringstream& out,
     const Bool verbose, const vecbasevector& genome2, const int LG,
     const vec< vec< pair<int,int> > > & Glocs )
{
     // Heuristics.

     const int F = 20;
     const int flank = 200;
     const int max_iterations = 10000;
     const int min_overlap_add = 5;
     const int bandwidth = 10;

     // First we define an extended version of the patch.

     int start1 = e.Start1( ), stop2 = e.Stop2( );
     ForceAssertEq( e.Nreps( ), 1 );
     basevector& patch = e.Rep(0); 
     basevector left, right;
     left.SetToSubOf( LEFT, 0, start1 );
     right.SetToSubOf( RIGHT, stop2, (int) RIGHT.size( ) - stop2 );
     int left_ext = Min( flank, left.isize( ) );
     int right_ext = Min( flank, right.isize( ) );
     basevector epatch = Cat( basevector( left, left.isize( ) - left_ext, left_ext ),
          patch, basevector( right, 0, right_ext ) );
     if (verbose)
     {    out << "\n";
          patch.Print( out, "patch" );
          epatch.Print( out, "epatch" );    }

     // Now find all the error-corrected fragment reads whose first 20 bases matches
     // the extended patch perfectly.

     vec< pair<int64_t,int> > fw_hits, rc_hits;
     vec<ho_interval> cov;
     kmer<F> f;
     for ( int pass = 1; pass <= 2; pass++ )
     {    for ( int j = 0; j <= epatch.isize( ) - F; j++ )
          {    f.SetToSubOf( epatch, j );
               if ( pass == 2 ) f.ReverseComplement( );
               int64_t low = LowerBound( fheads, f ), high = UpperBound( fheads, f );
               for ( int64_t i = low; i < high; i++ )
                    ( pass == 1 ? fw_hits : rc_hits ).push( fids[i], j );    }    }

     // Find aligned pairs.  Along with the original patch, we put all the reads
     // into a pot called seq.  We call them all 'reads'.

     vecbasevector seq;
     vec<int64_t> seq_id;
     seq.push_back(epatch);
     seq_id.push_back(-1);
     set<int64_t> used;
     for ( int i1 = 0; i1 < fw_hits.isize( ); i1++ )
     for ( int i2 = 0; i2 < rc_hits.isize( ); i2++ )
     {    int64_t id1 = fw_hits[i1].first, id2 = rc_hits[i2].first;
          if ( Member( used, id1 ) ) continue;
          if ( fpairs.getPartnerID(id1) != id2 ) continue;
          int start = fw_hits[i1].second, stop = rc_hits[i2].second + F;
          if ( stop < start ) continue;
          used.insert(id1);
          if (verbose) out << "using read " << id1 << "\n";
          if (verbose) out << "using read " << id2 << "\n";
          seq.push_back_reserve( fbases[id1] );
          basevector b = fbases[id2];
          b.ReverseComplement( );
          seq.push_back_reserve(b);
          seq_id.push_back( id1, id2 );    }

     // Find alignments.  We keep only the best alignment of a pair.  In case of a
     // tie, we keep the 'first' one, in some ill-defined sense.  Note that we don't
     // check for excess separation of the pair, and we don't use the quality 
     // scores.

     vec< pair<int64_t,int64_t> > id1_id2;
     vec< pair<align,align> > aligns;
     vec<int> errs;
     for ( int i1 = 0; i1 < fw_hits.isize( ); i1++ )
     for ( int i2 = 0; i2 < rc_hits.isize( ); i2++ )
     {    int64_t id1 = fw_hits[i1].first, id2 = rc_hits[i2].first;
          if ( fpairs.getPartnerID(id1) != id2 ) continue;
          int start = fw_hits[i1].second, stop = rc_hits[i2].second + F;
          if ( stop < start ) continue;
          id1_id2.push( id1, id2 );
          align a1, a2;
          int errors1, errors2;
          SmithWatBandedA( 
               fbases[id1], epatch, -start, bandwidth, a1, errors1, 0, 1, 1 );
          CenterMobileGaps( a1, fbases[id1], epatch, false, out );
          basevector b = fbases[id2];
          b.ReverseComplement( );
          SmithWatBandedA( 
               b, epatch, -( stop - b.isize( ) ), bandwidth, a2, errors2, 0, 1, 1 );
          CenterMobileGaps( a2, b, epatch, false, out );
          aligns.push( a1, a2 );
          errs.push_back( errors1 + errors2 );    }
     SortSync( id1_id2, errs, aligns );

     // Get best alignments.

     vec<align> seq_aligns( seq.size( ) );
     avector<int> gaps(1), lengths(1);
     gaps(0) = 0, lengths(0) = epatch.size( );
     seq_aligns[0] = align( 0, 0, gaps, lengths );
     for ( int i = 0; i < id1_id2.isize( ); i++ )
     {    int j = id1_id2.NextDiff(i);
          vec<int> this_errs, ids( j - i, vec<int>::IDENTITY );
          for ( int k = i; k < j; k++ )
               this_errs.push_back( errs[k] );
          SortSync( this_errs, ids );
          int64_t id1 = id1_id2[i].first, id2 = id1_id2[i].second;
          int p1 = Position(seq_id, id1), p2 = Position(seq_id, id2);
          seq_aligns[p1] = aligns[ i + ids[0] ].first;
          seq_aligns[p2] = aligns[ i + ids[0] ].second;
          i = j - 1;    }

     // Print alignments. (NO LONGER FUNCTIONAL.)

     /*
     for ( int z = 0; z < aligns.isize( ); z++ )
     {    int64_t id1 = id1_id2[z].first, id2 = id1_id2[z].second;
          const align &a1 = aligns[z].first, &a2 = aligns[z].second;
          if (verbose)
          {    out << "\n";
               PRINT_TO( out, id1 );
               out << a1.pos2( ) << "-" << a1.Pos2( ) << "\n";
               PrintVisualAlignment(
                    True, out, fbases[id1], epatch, a1, fquals[id1] );    }
          basevector b = fbases[id2];
          b.ReverseComplement( );
          PRINT_TO( out, id2 );
          out << a2.pos2( ) << "-" << a2.Pos2( ) << "\n";
          if (verbose) PrintVisualAlignment( True, out, b, epatch, a2 );    }
     */

     // Form the unipath graph.

     vecKmerPath paths, pathsrc, unipaths;
     vec<tagged_rpint> pathsdb, unipathsdb;
     ReadsToPathsCoreY( seq, F, paths, pathsrc, pathsdb );
     Unipath( paths, pathsrc, pathsdb, unipaths, unipathsdb );
     KmerBaseBroker kbb( F, paths, pathsrc, pathsdb, seq );
     vecbasevector unibases;
     for ( size_t i = 0; i < unipaths.size( ); i++ )
          unibases.push_back_reserve( kbb.Seq( unipaths[i] ) );

     // Unroll each read along the unipath graph.

     vec< vec< pair<int,int> > > unipos( seq.size( ) );
     for ( int i = 0; i < (int) seq.size( ); i++ )
     {    
          // First find the unipaths that comprise the read.

          vec<int> u;
          const KmerPath& p = paths[i];
          vec< pair<int,int> > uo;
          for ( int j = 0; j < p.NSegments( ); j++ )
          {    const KmerPathInterval& I = p.Segment(j);
               vec<longlong> locs;
               Contains( unipathsdb, I, locs );
               for ( int u = 0; u < locs.isize( ); u++ )
               {    const tagged_rpint& t = unipathsdb[ locs[u] ];
                    int uid = t.PathId( );
                    longlong offset = t.Start( ) - I.Start( );
                    for ( int r = 0; r < j; r++ )
                         offset += p.Segment(r).Length( );
                    for ( int r = 0; r < t.PathPos( ); r++ )
                         offset -= unipaths[uid].Segment(r).Length( );
                    uo.push_back( make_pair( int(offset), uid ) );    }    }
          UniqueSort(uo);
          for ( int j = 0; j < uo.isize( ); j++ )
               u.push_back( uo[j].second );

          // Find the read start on the first unipath.

          String u_s = unibases[ u[0] ].ToString( );
          String r_s = basevector( seq[i], 0, F ).ToString( );
          int read_start = u_s.Position(r_s);

          // Now assign coordinates to these unipaths.

          vec<int> coord;
          const align& a = seq_aligns[i];
          // PRINT4_TO( out, a.pos1( ), a.Pos1( ), a.pos2( ), a.Pos2( ) );

          coord.push_back( a.pos2( ) - read_start );


          vec<int> to_p2;
          int p1 = a.pos1( ), p2 = a.pos2( );
          to_p2.push_back(p2);
          for ( int j = 0; j < a.Nblocks( ); j++ ) 
          {    if ( a.Gaps(j) > 0 ) p2 += a.Gaps(j);
               if ( a.Gaps(j) < 0 ) 
               {    p1 -= a.Gaps(j);
                    for ( int l = 0; l < -a.Gaps(j); l++ )
                         to_p2.push_back(p2);    }
               for ( int x = 0; x < a.Lengths(j); x++ ) 
               {    ++p1; ++p2;
                    to_p2.push_back(p2);    }    }

          // Let start[l] denote the implied start position of the lth unibase on 
          // the read.

          int start = -read_start;
          for ( int j = 1; j < u.isize( ); j++ )
          {    start += unibases[ u[j-1] ].isize( ) - (F-1);
               if ( start < 0 ) coord.push_back( to_p2[0] + start );
               else if ( start >= to_p2.isize( ) ) // should be rare
               {    coord.push_back(               // also maybe should just punt
                         to_p2.back( ) + start - to_p2.isize( ) + 1 );    }
               else coord.push_back( to_p2[start] );    }

          if (verbose)
          {    out << "read " << i << " = " << seq_id[i] << " -->";
               for ( int j = 0; j < u.isize( ); j++ )
                    out << " " << u[j] << "." << coord[j];
               out << "\n\n";    }

          for ( int j = 0; j < u.isize( ); j++ )
               unipos[i].push( u[j], coord[j] );    }

     // Form the unipos elements into a graph.  Note that the way we assign 
     // positions could probably be improved.

     const int max_dist_delta = 5;
     vec< pair<int,int> > verts_init;
     vec< triple<int,int,Bool> > verts;
     vec<int> verts_map;
     int nseq = seq.size( );
     for ( int id = 0; id < nseq; id++ )
     {    for ( int j = 0; j < unipos[id].isize( ); j++ )
               verts_init.push( unipos[id][j].first, unipos[id][j].second );    }
     UniqueSort(verts_init);
     for ( int i = 0; i < verts_init.isize( ); i++ )
     {    int j;
          for ( j = i + 1; j < verts_init.isize( ); j++ )
          {    if ( verts_init[j].first != verts_init[i].first
                    || verts_init[j].second - verts_init[j-1].second 
                    > max_dist_delta )
               {    break;    }    }
          for ( int k = i; k < j; k++ )
               verts_map.push_back( verts.size( ) );
          Bool orig = False;
          for ( int k = i; k < j; k++ )
               if ( Member( unipos[0], verts_init[k] ) ) orig = True;
          int m = i + (j-i)/2;
          verts.push( verts_init[m].first, verts_init[m].second, orig );
          i = j - 1;    }
     vec< vec<int> > sux0(nseq);
     vec< vec<int> > from( verts.size( ) ), to( verts.size( ) );
     for ( int id = 0; id < unipos.isize( ); id++ )
     {    for ( int j = 0; j < unipos[id].isize( ); j++ )
          {    int v = verts_map[ BinPosition( verts_init, unipos[id][j] ) ];
               sux0[id].push_back(v);    }    }
     for ( int id = 0; id < unipos.isize( ); id++ )
     {    for ( int j = 0; j < unipos[id].isize( ) - 1; j++ )
          {    int v1 = sux0[id][j], v2 = sux0[id][j+1];
               from[v1].push_back(v2), to[v2].push_back(v1);    }    }
     for ( int v = 0; v < verts.isize( ); v++ )
     {    UniqueSort(from[v]), UniqueSort(to[v]);    }
     digraphV< triple<int,int,Bool> > M( from, to, verts );
     vec< vec< triple<int,int,Bool> > > sux2(nseq);
     for ( int id = 0; id < nseq; id++ )
     {    for ( int j = 0; j < sux0[id].isize( ); j++ )
               sux2[id].push_back( M.Vert( sux0[id][j] ) );    }
     if (verbose)
     {    out << "\nvertices:\n";
          for ( int v = 0; v < verts.isize( ); v++ )
               out << verts[v].first << "." << verts[v].second << "\n";
          out << "\n";
          out << "\nsux2:\n";
          for ( int id = 0; id < nseq; id++ )
          {    out << "id = " << id << ":";
               for ( int j = 0; j < sux2[id].isize( ); j++ )
               {    out << " " << sux2[id][j].first << "."
                         << sux2[id][j].second;    }
               out << "\n";    }    }

     // Trim graph.

     Trim( sux2, M );
     
     // Get coverage.

     vec<double> covx2;
     GetPairCoverage( sux2, M, covx2 );

     // Trim to remove stupid stuff.  We look for a sequence of four original
     // vertices, which we assume all have at least a certain minimum coverage.
     // Then we take the middle two, and get their positions relative to the
     // origin, which we call low and high.  Then we kill any non-original vertex
     // that is positioned between low and high.

     vec<int> stupid2;
     const double min_covs2 = 10.0;
     for ( int v = 0; v < M.N( ); v++ )
     {    if ( !M.Vert(v).third ) continue;
          if ( covx2[v] < min_covs2 ) continue;
          vec<int> mini;
          int x = v;
          mini.push_back(x);
          for ( int m = 0; m < 3; m++ )
          {    for ( int j = 0; j < M.From(x).isize( ); j++ )
               {    int y = M.From(x)[j];
                    if ( M.Vert(y).third && covx2[y] >= min_covs2 )
                    {    x = y;
                         mini.push_back(x);    }    }    }
          if ( mini.size( ) == 4 )
          {    int x = mini[1], y = mini[2];
               int low = M.Vert(x).second, high = M.Vert(y).second;
               for ( int w = 0; w < M.N( ); w++ )
               {    int d = M.Vert(w).second;
                    if ( low <= d && d <= high && !M.Vert(w).third )
                         stupid2.push_back(w);    }    }    }
     UniqueSort(stupid2);
     M.DeleteVertices(stupid2);
     Trim( sux2, M );
     GetPairCoverage( sux2, M, covx2 );

     // Remove bubbles.

     RemoveBubbles( sux2, M, verbose, out );
     GetPairCoverage( sux2, M, covx2 );

     // Pick branches and again remove bubbles.

     PickBranches( sux2, M );
     RemoveBubbles( sux2, M, verbose, out );
     GetPairCoverage( sux2, M, covx2 );

     // Print graph.

     if (verbose)
     {    vec<String> vertex_labels;
          for ( int i = 0; i < M.N( ); i++ )
          {    ostringstream out;
               out << setiosflags(ios::fixed) << setprecision(1) << covx2[i]
                    << resetiosflags(ios::fixed);
               vertex_labels.push_back( ToString(i) + "="
                    + ToString( M.Vert(i).first ) + "." 
                    + ToString( M.Vert(i).second ) + "[" + out.str( ) + "]" 
                    + ( M.Vert(i).third ? "###" : "" ) );    }
          M.DOT_vl( out, vertex_labels, "", vec<String>( ) );
          flush(out);    }

     // For convenience, write down the lengths of the unipaths in kmers.

     vec<int> dx( unibases.size( ) );
     for ( int i = 0; i < (int) unibases.size( ); i++ )
          dx[i] = unibases[i].isize( ) - (F-1);

     // Compare to truth data.

     int components = M.NComponents( ), branch_points = 0;
     for ( int v = 0; v < M.N( ); v++ )
          if ( M.From(v).size( ) > 1 || M.To(v).size( ) > 1 ) branch_points++;
     int start = -1, stop = -1;
     for ( int j = 0; j < M.N( ); j++ )
     {    if ( M.Vert(j) == sux2[0].front( ) ) 
          {    start = j;
               break;    }    }
     for ( int j = M.N( ) - 1; j >= 0; j-- )
     {    if ( M.Vert(j) == sux2[0].back( ) ) 
          {    stop = j;
               break;    }    }
     int ustart = M.Vert(start).first;
     if (verbose)
     {    out << "\nGraph has " << components << " components" 
               << " and " << branch_points << " branch points.\n";
          PRINT3_TO( out, start, stop, ustart );    }
     if ( verbose && genome2.size( ) > 0 )
     {    vec< vec<placementx> > locs( unibases.size( ) );
          for ( int i = 0; i < (int) unibases.size( ); i++ )
          {    locs[i] = FindGenomicPlacements( 
                    unibases[i], LG, genome2, Glocs );    }
          vec< vec<int> > places;
          for ( int i = 0; i < locs[ustart].isize( ); i++ )
          {    placementx p = locs[ustart][i];
               vec<int> place;
               place.push_back(ustart);
               int lj = ustart;
               restart:
               Bool extended = False;
               for ( int j = 0; j < (int) unibases.size( ); j++ )
               {    for ( int k = 0; k < locs[j].isize( ); k++ )
                    {    placementx pn = locs[j][k];
                         if ( pn.g != p.g || pn.fw != p.fw ) continue;
                         if ( pn.fw && p.pos + dx[lj] == pn.pos )
                         {    extended = True;
                              place.push_back(j);
                              p = pn;
                              lj = j;
                              goto next;    }
                         if ( !pn.fw && pn.pos + dx[j] == p.pos )
                         {    extended = True;
                              place.push_back(j);
                              p = pn;
                              lj = j;
                              goto next;    }    }    }
               next: 
               if ( !extended ) places.push_back(place);
               else goto restart;    }
          UniqueSort(places);
          out << "\nThere are " << places.size( ) << " genomic placements:\n";
          for ( int j = 0; j < places.isize( ); j++ )
          {    const vec<int>& p = places[j];

               // Display the placement, mapping to the graph.  Determine the 
               // position f at which the placement fails to map to the assembly 
               // graph.
     
               out << "\n[" << j+1 << "]";
               int v = start, f;
               for ( f = 1; f < p.isize( ); f++ )
               {    int l, w = -1;
                    out << ( f % 8 == 0 ? "\n" : " " ) << p[f-1] << "(" << v << ")";
                    for ( l = 0; l < M.From(v).isize( ); l++ )
                    {    w = M.From(v)[l];
                         if ( M.Vert(w).first == p[f] ) break;    }
                    if ( l == M.From(v).isize( ) ) break;
                    v = w;    }
               for ( int l = f; l < p.isize( ); l++ )
                    out << ( l % 8 == 0 ? "\n" : " " ) << p[l];
               out << "\n";    }

          out << endl;    }

     // Make edits.  We give up if the graph is not a line.

     if ( components != 1 || branch_points != 0 ) return False;
     int vv = start;
     patch = unibases[ M.Vert(vv).first ];
     while( M.From(vv).nonempty( ) )
     {    vv = M.From(vv)[0];
          patch.resize( patch.isize( ) - (F-1) );
          patch = Cat( patch, unibases[ M.Vert(vv).first ] );    }
     if (verbose) patch.Print( out, "patch0" );
     if (verbose) 
          PRINT2_TO( out, M.Vert(start).first, M.Vert(stop).first );

     int u1 = -1, u2 = -1, start1e = -1, stop2e = -1;
     String s1 = basevector( epatch, 0, F ).ToString( );
     String s2 = basevector( epatch, epatch.isize( ) - F, F ).ToString( );
     for ( int u = 0; u < (int) unibases.size( ); u++ )
     {    String s = unibases[u].ToString( );
          if ( s.Contains(s1) )
          {    u1 = u;
               start1e = s.Position(s1);    }
          if ( s.Contains(s2) )
          {    u2 = u;
               stop2e = s.Position(s2) + F;    }    }

     int patch_start = left_ext + start1e;
     int patch_stop = patch.isize( ) - right_ext
          - ( unibases[ M.Vert(stop).first ].isize( ) - stop2e );
     if ( patch_start > patch_stop || patch_stop > patch.isize( ) )
     {    if (verbose)
          {    out << "\npatch boundaries don't make sense\n";
               PRINT5_TO( out, patch.size( ), left_ext, right_ext, start1e, stop2e );
               PRINT_TO( out, unibases[ M.Vert(stop).first ].size( ) );    }
          return False;    }
     patch = basevector( patch, left_ext + start1e, 
          patch.isize( ) - left_ext - right_ext - start1e
               - ( unibases[ M.Vert(stop).first ].isize( ) - stop2e ) );
     if (verbose) patch.Print( out, "patch" );
     return True;    }

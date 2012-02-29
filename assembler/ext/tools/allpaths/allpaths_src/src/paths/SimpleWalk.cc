/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include <queue>

#include "CoreTools.h"
#include "math/Functions.h"
#include "graph/Digraph.h"
#include "graph/DigraphTemplate.h"
#include "paths/PairedPair.h"
#include "paths/SimpleWalk.h"

class walkleaf {

     public:

     walkleaf( ) { }
     walkleaf( int e, int v, int plen, Bool open ) 
          : e(e), v(v), plen(plen), open(open) { }

     int e;     // edge id
     int v;     // vertex id at source of edge
     int plen;  // length of path through e from root, in kmers
     Bool open; // if open for further extension

     friend Bool operator<( const walkleaf& l1, const walkleaf& l2 )
     {    return l1.plen < l2.plen;    }

};

void SimpleWalkRight( const pp_pair& p, const vec<pp_read>& reads, 
     const vec<int>& L, const double dmult, vec<pp_closure>& closures,
     const int max_opens, Bool create_closures_if_fail, Bool& fail, int verbosity,
     Bool depth_first, int max_nodes )
{    
     // Sanity checks, etc.

     ForceAssertGt( p.LeftSize( ), 0 );
     double clock = 0.0;
     if ( verbosity >= 1 ) clock = WallClockTime( );
     fail = False;
     closures.clear( );
     int max_opens_seen = 1;
     if ( verbosity >= 2 ) cout << "Entering SimpleWalkRight, p = " << p << "\n";

     // Set up index.

     vec< vec< pair<int,int> > > rind( L.size( ) );
     for ( int i = 0; i < reads.isize( ); i++ )
     {    const pp_read& r = reads[i];
          for ( int j = 0; j < r.isize( ) - 1; j++ )
               rind[ r[j] ].push_back( make_pair( i, j ) );    }

     // Determine longest length of path (walkleaf.plen) that we allow.

     int maxgap = int( floor( p.Gap( ) + dmult * p.Dev( ) ) );
     int leftsize = 0;
     for ( int i = 0; i < p.LeftSize( ); i++ )
          leftsize += L[ p.Left(i) ];
     int maxplen = maxgap + leftsize;

     // We keep a tree T as a digraphE.  This is presumably inefficient.

     digraphE<walkleaf> T;
     T.AddVertices(1);

     // In depth first case, keep track of open leaves via a priority queue.

     priority_queue<walkleaf> open_leaves;

     // Load p.Left onto the tree.

     int g = 0;
     for ( int i = 0; i < p.LeftSize( ); i++ )
     {    g += L[ p.Left(i) ];
          T.AddVertices(1);
          Bool open = ( i == p.LeftSize( ) - 1 && g <= maxplen );
          T.AddEdge( i, i+1, walkleaf( p.Left(i), i, g, open ) );    }
     if (depth_first) open_leaves.push( T.EdgeObject( p.LeftSize( ) - 1 ) );

     // Iteratively grow the tree.

     int lf = 0, opens = 1;
     while(1)
     {
          // Find an active leaf x.

          if ( max_nodes > 0 && T.N( ) > max_nodes )
          {    fail = True;
               goto finale;    }
          walkleaf x;
          int l = 0;
          if (depth_first)
          {    if ( open_leaves.empty( ) ) break;
               x = open_leaves.top( );
               open_leaves.pop( );
               int v = x.v;
               for ( int j = 0; j < T.From(v).isize( ); j++ )
               {    if ( T.EdgeObjectByIndexFrom( v, j ).e == x.e )
                    {    l = T.From(v)[j];
                         T.EdgeObjectByIndexFromMutable( v, j ).open = False;
                         break;    }    }    }
          else
          {    for ( lf++; lf < T.N( ); lf++ )
               {    const walkleaf& x = T.EdgeObjectByIndexTo( lf, 0 );
                    if ( x.open ) break;    }
               if ( lf == T.N( ) ) break;
               l = lf;
               T.EdgeObjectByIndexToMutable( l, 0 ).open = False;
               x = T.EdgeObjectByIndexToMutable( l, 0 );    }
          --opens;

          /*
          // Find extensions.

          static vec<int> nexts;
          nexts.clear( );
          for ( int u = 0; u < rind[x.e].isize( ); u++ )
          {    int i = rind[x.e][u].first, j = rind[x.e][u].second;
               const pp_read& r = reads[i];
               {    int m = x.v;
                    Bool failx = False;
                    if ( m != 0 )
                    {    walkleaf y = T.EdgeObjectByIndexTo( m, 0 );
                         for ( int k = j - 1; k >= 0; k-- )
                         {    if ( y.e != r[k] )
                              {    failx = True;
                                   break;    }
                              m = y.v;
                              if ( m == 0 ) break;
                              y = T.EdgeObjectByIndexTo( m, 0 );    }    }
                    if ( !failx ) nexts.push_back( r[j+1] );    }    }
          UniqueSort(nexts);

          // Add nodes to graph.

          for ( int i = 0; i < nexts.isize( ); i++ )
          {    int n = nexts[i];
               int v = T.N( );
               if ( verbosity >= 1 && v % 100000 == 0 )
               {    cout << "v = " << v << ", opens = " << opens
                         << ", time used = " << TimeSince(clock) << endl;    }
               T.AddVertices(1);
               int plen = x.plen + L[n];
               T.AddEdge( l, v, walkleaf( n, l, plen, plen <= maxplen ) );    
               if ( plen <= maxplen ) 
               {    ++opens;    
                    if (depth_first)
                         open_leaves.push( walkleaf( n, l, plen, True ) );    }
               max_opens_seen = Max( max_opens_seen, opens );
               if ( max_opens > 0 && opens > max_opens )
               {    fail = True;
                    goto finale;    }    }    }
          */

/* ...............................................................................*/

          // Find reads that extend the current node.

          static vec<int> nexts;
          nexts.clear( );
          for ( int u = 0; u < rind[x.e].isize( ); u++ )
          {    int i = rind[x.e][u].first, j = rind[x.e][u].second;
               const pp_read& r = reads[i];
               {    int m = x.v;
                    Bool failx = False;
                    if ( m != 0 )
                    {    walkleaf y = T.EdgeObjectByIndexTo( m, 0 );
                         for ( int k = j - 1; k >= 0; k-- )
                         {    if ( y.e != r[k] )
                              {    failx = True;
                                   break;    }
                              m = y.v;
                              if ( m == 0 ) break;
                              y = T.EdgeObjectByIndexTo( m, 0 );    }    }
                    if ( !failx ) nexts.push_back(u);    }    }

          // Find the extensions implied by the aligning reads.

          vec< vec<int> > exts;
          for ( int z = 0; z < nexts.isize( ); z++ )
          {    int u = nexts[z];
               int i = rind[x.e][u].first, j = rind[x.e][u].second;
               const pp_read& r = reads[i];
               static vec<int> e;
               e.clear( );
               int plen = x.plen;
               for ( int q = j + 1; q < r.isize( ); q++ )
               {    e.push_back( r[q] );
                    plen += L[ r[q] ];
                    if ( plen > maxplen ) break;    }
               exts.push_back(e);    }
          UniqueSort(exts);

          // Remove extensions that are extensions of other extensions.

          vec<Bool> remove( exts.size( ), False );
          for ( int i1 = 0; i1 < exts.isize( ); i1++ )
          {    if ( remove[i1] ) continue;
               for ( int i2 = i1 + 1; i2 < exts.isize( ); i2++ )
               {    if ( remove[i2] ) continue;
                    if ( exts[i2].size( ) > exts[i1].size( ) )
                    {    int j;
                         for ( j = 0; j < exts[i1].isize( ); j++ )
                              if ( exts[i1][j] != exts[i2][j] ) break;
                         if ( j == exts[i1].isize( ) ) 
                              remove[i2] = True;    }    }    }
          EraseIf( exts, remove );

          // Add nodes to graph.

          for ( int i = 0; i < exts.isize( ); i++ )
          {    int vlast = l;
               const vec<int>& e = exts[i];
               int plen = x.plen;
               for ( int j = 0; j < e.isize( ); j++ )
               {    plen += L[ e[j] ];
                    int v = T.N( );
                    if ( verbosity >= 1 && v % 100000 == 0 )
                    {    cout << "v = " << v << ", opens = " << opens
                              << ", time used = " << TimeSince(clock) << endl;    }
                    T.AddVertices(1);
                    Bool open = ( j == e.isize( ) - 1 && plen <= maxplen );
                    T.AddEdge( vlast, v, walkleaf( e[j], vlast, plen, open ) );    
                    if (open) 
                    {    ++opens;    
                         if (depth_first)
                              open_leaves.push( 
                                   walkleaf( e[j], vlast, plen, True ) );    }
                    vlast = v;    }
               max_opens_seen = Max( max_opens_seen, opens );
               if ( max_opens > 0 && opens > max_opens )
               {    fail = True;
                    goto finale;    }    }    }

/*................................................................................*/

     finale:
     if ( fail && !create_closures_if_fail ) 
     {    cout << "\nv = " << T.N( ) << ", opens = " << max_opens_seen << "\n";
          return;    }

     // Get the paths.

     vec<pp_read> paths;
     vec<int> plens;
     for ( int l = 0; l < T.N( ); l++ )
     {    if ( T.From(l).empty( ) )
          {    static pp_read q;
               q.clear( );
               walkleaf y = T.EdgeObjectByIndexTo( l, 0 );
               int plen = y.plen;
               int m = l;
               while(1)
               {    q.push_back(y.e);
                    m = y.v;
                    if ( m == 0 ) break;
                    y = T.EdgeObjectByIndexTo( m, 0 );    }
               q.ReverseMe( );
               paths.push_back(q);
               plens.push_back(plen);    }    }
     if ( verbosity == 2 ) cout << paths.size( ) << " paths found\n";
     if ( verbosity >= 4 )
     {    cout << paths.size( ) << " paths found:\n";
          for ( int i = 0; i < paths.isize( ); i++ )
               cout << "[" << i << "] " << paths[i] << "\n";    }

     // Find and validate the closures.

     for ( int i = 0; i < paths.isize( ); i++ )
     {    int plen = plens[i];
          pp_pair q( paths[i], p.Right( ), p.Gap( ) - plen + leftsize, p.Dev( ) );
          static vec<pp_closure> qclosures;
          GetClosures( q, L, dmult, qclosures );
          for ( int j = 0; j < qclosures.isize( ); j++ )
          {    const pp_closure& c = qclosures[j];
               static vec<Bool> covered;
               covered.resize_and_set( c.size( ), False );
               for ( int k = 0; k < reads.isize( ); k++ )
               {    static vec<int> offsets;
                    GetOverlaps( c, reads[k], offsets );
                    for ( int u = 0; u < offsets.isize( ); u++ )
                    {    int o = offsets[u];
                         int start = Max( 0, o );
                         int stop = Min( c.isize( ), o + reads[k].isize( ) );
                         if ( stop == c.isize( ) 
                              && ( start < stop - 1 || start == 0 ) )
                         {    ++stop;    }
                         for ( int x = start; x < stop - 1; x++ )
                              covered[x] = True;    }    }
               if ( Sum(covered) == c.isize( ) ) closures.push_back(c);    }    }
     Sort(closures);
     if ( verbosity >= 1 ) 
     {    cout << "v = " << T.N( ) << ", opens = " << max_opens_seen
               << ", found " << closures.size( ) << " closures\n";    }    
     if ( verbosity >= 3 )
     {    for ( int i = 0; i < closures.isize( ); i++ )
               cout << "[" << i << "] " << closures[i] << "\n";   
          cout << "\n";    }
     if ( verbosity >= 1 )
          cout << "SimpleWalk time used = " << TimeSince(clock) << "\n";    }

template void digraphE<walkleaf>::AddEdge(int, int, walkleaf const&);

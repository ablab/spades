///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// This file contains some template functions from Digraph.h.  They are here in
// a separate file so that these functions do not have to be inlined, thereby 
// allowing for reduction of compilation time and executable size (in principle).
// Generally, the less places this file is included, the better.

#ifndef DIGRAPH_TEMPLATE_H
#define DIGRAPH_TEMPLATE_H

#include "CoreTools.h"
#include "Equiv.h"
#include "FeudalMimic.h"
#include "Set.h"
#include "VecUtilities.h"
#include "graph/Digraph.h"
#include <cstddef>

template<class E> void digraphE<E>::Initialize( 
     const digraphE& g, const vec<int>& v )
{    from_.resize( v.size( ) ), to_.resize( v.size( ) );
     from_edge_obj_.resize( v.size( ) ), to_edge_obj_.resize( v.size( ) );
     int edgecount = 0;
     vec<int> vsorted(v), vindex( v.size( ), vec<int>::IDENTITY );
     SortSync( vsorted, vindex );
     for ( int i = 0; i < v.isize( ); i++ )
     {    int x = v[i];
          for ( int j = 0; j < g.From(x).isize( ); j++ )
          {    int y = g.From(x)[j];
               int p2 = BinPosition( vsorted, y );
               if ( p2 < 0 ) continue;
               int i2 = vindex[p2];
               from_[i].push_back(i2);    
               to_[i2].push_back(i);
               from_edge_obj_[i].push_back(edgecount);
               to_edge_obj_[i2].push_back(edgecount);
               ++edgecount;    }    }
     edges_.reserve(edgecount);
     for ( int i = 0; i < v.isize( ); i++ )
     {    int x = v[i];
          for ( int j = 0; j < g.From(x).isize( ); j++ )
          {    int y = g.From(x)[j];
               int p2 = BinPosition( vsorted, y );
               if ( p2 < 0 ) continue;
               int i2 = vindex[p2];
               edges_.push_back( g.EdgeObjectByIndexFrom( x, j ) );    }    }
     for ( int i = 0; i < v.isize( ); i++ )
     {    SortSync( from_[i], from_edge_obj_[i] );
          SortSync( to_[i], to_edge_obj_[i] );    }    }

template<class E> digraphE<E>::digraphE( const digraphE& g, const vec<int>& v )
{    Initialize( g, v );    }

template<class E> digraphE<E>::digraphE( 
     const digraphE& g, const vec< vec<int> >& C )
{    int nedges = 0;
     for ( int i = 0; i < C.isize( ); i++ )
          nedges += C[i].size( );
     edges_.reserve(nedges);
     vec<int> to_left, to_right;
     g.ToLeft(to_left), g.ToRight(to_right);
     for ( int i = 0; i < C.isize( ); i++ )
     {    for ( int j = 0; j < C[i].isize( ); j++ )
               edges_.push_back( g.EdgeObject( C[i][j] ) );    }
     for ( int pass = 1; pass <= 2; pass++ )
     {    int nverts = 0, nedges = 0;
          for ( int i = 0; i < C.isize( ); i++ )
          {    vec<int> verts;
               for ( int j = 0; j < C[i].isize( ); j++ )
                    verts.push_back( to_left[ C[i][j] ], to_right[ C[i][j] ] );
               UniqueSort(verts);
               if ( pass == 2 )
               {    for ( int j = 0; j < C[i].isize( ); j++ )
                    {    int v = BinPosition( verts, to_left[ C[i][j] ] );
                         int w = BinPosition( verts, to_right[ C[i][j] ] );
                         from_[ nverts + v ].push_back( nverts + w );
                         to_[ nverts + w ].push_back( nverts + v );
                         from_edge_obj_[ nverts + v ].push_back(nedges + j);
                         to_edge_obj_[ nverts + w ].push_back(nedges + j);    }    }
               nverts += verts.size( );
               nedges += C[i].size( );    }
          if ( pass == 1 )
          {    from_.resize(nverts), to_.resize(nverts);
               from_edge_obj_.resize(nverts), to_edge_obj_.resize(nverts);    }    }
     for ( int v = 0; v < N( ); v++ )
     {    SortSync( from_[v], from_edge_obj_[v] );
          SortSync( to_[v], to_edge_obj_[v] );    }    }

template<class E> digraphE<E>::digraphE( const digraphE& g, int n )
{    equiv_rel e;
     g.ComponentRelation(e);
     vec<int> reps, o;
     e.OrbitRepsAlt(reps);
     ForceAssertLt( n, reps.isize( ) );
     e.Orbit( reps[n], o );
     Initialize( g, o );    }

template<class E> 
void BinaryWrite( int fd, const digraphE<E>& g )
{    
    BinaryWrite(fd, (const digraph&) g);
    int n = g.N();
    for ( int v = 0; v < n; v++ )
    {
        BinaryWrite( fd, g.from_edge_obj_[v] );
        BinaryWrite( fd, g.to_edge_obj_[v] );
    }
    int e = g.edges_.size();
    WriteBytes(fd, &e, sizeof(int));
    for ( int i = 0; i < e; i++ )
    {
        BinaryWrite( fd, g.EdgeObject(i) );
    }
}

template<class E> void BinaryRead( int fd, digraphE<E>& g )
{    BinaryRead( fd, (digraph&) g );
     int n = g.N( );
     g.from_edge_obj_.resize(n), g.to_edge_obj_.resize(n);
     for ( int v = 0; v < n; v++ )
     {    BinaryRead( fd, g.from_edge_obj_[v] );
          BinaryRead( fd, g.to_edge_obj_[v] );    }
     int e;
     ReadBytes( fd, &e, sizeof(int) );
     g.edges_.resize(e);
     for ( int i = 0; i < e; i++ )
          BinaryRead( fd, g.EdgeObjectMutable(i) );    }

template<class E> void digraphE<E>::EdgeEquivConstructor( 
     const vec<E>& edges, const equiv_rel& e )
{    edges_ = edges;
     int ne = edges.size( );
     vec<int> reps;
     e.OrbitReps(reps);
     int nv = 2 * reps.isize( );
     to_edge_obj_.resize(nv);
     from_edge_obj_.resize(nv);
     to_.resize(nv);
     from_.resize(nv);
     for ( int i = 0; i < reps.isize( ); i++ )
     {    vec<int> o;
          e.Orbit( reps[i], o );
          for ( int j = 0; j < o.isize( ); j++ )
          {    from_[ 2*i ].push_back( 2*i + 1 );
               from_edge_obj_[ 2*i ].push_back( o[j] );
               to_[ 2*i + 1 ].push_back( 2*i );
               to_edge_obj_[ 2*i + 1 ].push_back( o[j] );    }    }    }

template<class E> digraphE<E>::digraphE( const vec<E>& edges, const equiv_rel& e )
{    EdgeEquivConstructor( edges, e );    }

template<class E> digraphE<E>::digraphE( const vec<E>& edges, 
     const ConstructorBehavior constructor_type )
{    edges_ = edges;
     int ne = edges.size( );
     int nv = ( constructor_type == EDGES_SEPARATE ? ne * 2 : ne + 1 );
     to_edge_obj_.resize(nv);
     from_edge_obj_.resize(nv);
     to_.resize(nv);
     from_.resize(nv);
     if ( constructor_type == EDGES_SEPARATE )
     {    for ( int i = 0; i < ne; i++ )
          {    from_[ 2*i ].push_back( 2*i + 1 );
               from_edge_obj_[ 2*i ].push_back(i);
               to_[ 2*i + 1 ].push_back( 2*i );
               to_edge_obj_[ 2*i + 1 ].push_back(i);    }    }
     else if ( constructor_type == EDGES_IN_LINE )
     {    for ( int i = 0; i <= ne; i++ )
          {    if ( i < ne )
               {    from_[i].push_back(i+1);
                    from_edge_obj_[i].push_back(i);    }
               if ( i > 0 ) 
               {    to_[i].push_back(i-1);
                    to_edge_obj_[i].push_back(i-1);    }    }    }
     else ForceAssert( 0 == 1 );    }

template<class E> void digraphE<E>::Used( vec<Bool>& used ) const
{    used.resize_and_set( edges_.size( ), False );
     for ( int i = 0; i < N( ); i++ )
     {    for ( int j = 0; j < to_edge_obj_[i].isize( ); j++ )
               used[ to_edge_obj_[i][j] ] = True;    }    }

template<class E> void digraphE<E>::JoinEdges( int x, const E& e )
{    ForceAssert( from_[x].size( ) == 1 && to_[x].size( ) == 1 );
     int v = to_[x][0], w = from_[x][0];
     ForceAssert( x != v || x != w );
     from_[x].clear( ), from_edge_obj_[x].clear( );
     to_[x].clear( ), to_edge_obj_[x].clear( );
     for ( int i = 0; i < from_[v].isize( ); i++ )
     {    if ( from_[v][i] == x )
          {    from_[v].erase( from_[v].begin( ) + i );
               from_edge_obj_[v].erase( from_edge_obj_[v].begin( ) + i );
               break;    }    }
     for ( int i = 0; i < to_[w].isize( ); i++ )
     {    if ( to_[w][i] == x )
          {    to_[w].erase( to_[w].begin( ) + i );
               to_edge_obj_[w].erase( to_edge_obj_[w].begin( ) + i );
               break;    }    }
     AddEdge( v, w, e );    }

template<class E> void digraphE<E>::RemoveUnneededVertices( )
{    for ( int i = 0; i < N( ); i++ )
     {    if ( From(i).size( ) == 1 && To(i).size( ) == 1 && From(i)[0] != i )
          {    E p = EdgeObjectByIndexTo( i, 0 );
               p.append( EdgeObjectByIndexFrom( i, 0 ) );
               JoinEdges( i, p );    }    }
     RemoveEdgelessVertices( );    }

// Input is a set of vertices v.  Each v must be located at the opening of
// a bubble, with exactly two edges that lead to the same successor w:
//        _-_
//  --> v     w -->
//        -_-
template<class E> void digraphE<E>::PopBubbles( const vec<int> & bubble_vs )
{
  vec<int> bubble_edges;
  bubble_edges.reserve( bubble_vs.size() );
  
  for ( int i = 0; i < bubble_vs.isize(); i++ ) {
    int v = bubble_vs[i];
    ForceAssertEq( from_[v].size(), 2u );
    ForceAssertEq( from_[v][0], from_[v][1] );
    
    // Choose one of the edges that make up this bubble, and delete it.
    // Arbitrarily, we choose the higher-indexed path.
    bubble_edges.push_back( from_edge_obj_[v][1] );
  }
  
  DeleteEdges( bubble_edges );
  // Combine edges.  For bubbles v->w in which v had only one predecessor
  // and/or w had only one successor, this will combine the remaining edge in
  // the bubble with the edge leading to/from the bubble.
  RemoveUnneededVertices( );
  // Clear out edges that have been removed from the graph.
  RemoveDeadEdgeObjects( );
}

template<class E> void digraphE<E>::RemoveEdgelessVertices( 
     const vec<int>& to_remove )
{    vec<Bool> remove( N( ), False );
     for ( int i = 0; i < to_remove.isize( ); i++ )
          remove[ to_remove[i] ] = True;
     vec<int> new_vertex_id( N( ), -1 );
     int id = 0;
     for ( int i = 0; i < N( ); i++ )
     {    if ( remove[i] )
          {    ForceAssert( from_[i].empty( ) );
               ForceAssert( to_[i].empty( ) );    }
          else
          {    new_vertex_id[i] = id;
               ++id;    }    }
     EraseIf( from_, remove ), EraseIf( from_edge_obj_, remove );
     EraseIf( to_, remove ), EraseIf( to_edge_obj_, remove );
     for ( int i = 0; i < N( ); i++ )
     {    for ( int j = 0; j < from_[i].isize( ); j++ )
               from_[i][j] = new_vertex_id[ from_[i][j] ];
          for ( int j = 0; j < to_[i].isize( ); j++ )
               to_[i][j] = new_vertex_id[ to_[i][j] ];    }    }

template<class E> void digraphE<E>::RemoveEdgelessVertices( )
{    vec<int> to_remove;
     for ( int i = 0; i < N( ); i++ )
          if ( from_[i].empty( ) && to_[i].empty( ) ) to_remove.push_back(i);
     RemoveEdgelessVertices(to_remove);    }

template<class E> void digraphE<E>::Reverse( )
{    for ( int i = 0; i < N( ); i++ )
     {    swap( from_[i], to_[i] );
          swap( from_edge_obj_[i], to_edge_obj_[i] );    }    }

template<class E> void digraphE<E>::ReverseComponent( int x )
{    equiv_rel e( N( ) );
     for ( int v = 0; v < N( ); v++ )
     {    for ( int i = 0; i < from_[v].isize( ); i++ )
          {    int w = from_[v][i];
               e.Join( v, w );    }    }
     vec<int> o;
     e.Orbit( x, o );
     for ( int j = 0; j < o.isize( ); j++ )
     {    int i = o[j];
          swap( from_[i], to_[i] );
          swap( from_edge_obj_[i], to_edge_obj_[i] );    }    }

template<class E> void digraphE<E>::ReorderVertices( const vec<int>& new_order )
{    ForceAssertEq( new_order.isize( ), N( ) );
     vec<int> order_new( N( ) );
     for ( int i = 0; i < N( ); i++ )
          order_new[ new_order[i] ] = i;
     PermuteVec( from_, order_new );
     PermuteVec( from_edge_obj_, order_new );
     PermuteVec( to_, order_new );
     PermuteVec( to_edge_obj_, order_new );
     for ( int v = 0; v < N( ); v++ )
     {    for ( int j = 0; j < from_[v].isize( ); j++ )
               from_[v][j] = order_new[ from_[v][j] ];
          for ( int j = 0; j < to_[v].isize( ); j++ )
               to_[v][j] = order_new[ to_[v][j] ];
          SortSync( from_[v], from_edge_obj_[v] );
          SortSync( to_[v], to_edge_obj_[v] );    }    }

template<class E> void digraphE<E>::ReorderComponents( const vec<int>& new_order )
{    equiv_rel e( N( ) );
     for ( int v = 0; v < N( ); v++ )
     {    for ( int i = 0; i < from_[v].isize( ); i++ )
          {    int w = from_[v][i];
               e.Join( v, w );    }    }
     vec<int> reps;
     for ( int v = 0; v < N( ); v++ )
          if ( e.Representative(v) ) reps.push_back(v);
     ForceAssertEq( new_order.size( ), reps.size( ) );
     vec<int> new_vertex_order;
     for ( int i = 0; i < reps.isize( ); i++ )
     {    int v = reps[ new_order[i] ];
          vec<int> o;
          e.Orbit( v, o );
          new_vertex_order.append(o);    }
     ReorderVertices(new_vertex_order);    }

template<class E> void
digraphE<E>::ComponentEdges( vec< vec<edge_t> >& edges ) const
{
  vec<vec<int> > vertices;
  Components( vertices );
  int n = vertices.isize( );
  
  edges.resize( 0 );
  edges.resize( n );
  for ( int i = 0; i < n; i++ ) {
    for ( int j = 0; j < vertices[i].isize( ); j++ )
      edges[i].append( FromEdgeObj( vertices[i][j] ) );
    UniqueSort( edges[i] );
  }
}

template<class E> void digraphE<E>::Append( const digraphE<E>& D )
{    int nedges = edges_.size( );
     edges_.append( D.edges_ );
     int nvertices = from_.size( );
     from_.append( D.from_ );
     to_.append( D.to_ );
     from_edge_obj_.append( D.from_edge_obj_ );
     to_edge_obj_.append( D.to_edge_obj_ );
     for ( int i = nvertices; i < N( ); i++ )
     {    for ( int j = 0; j < from_[i].isize( ); j++ )
          {    from_[i][j] += nvertices;
               from_edge_obj_[i][j] += nedges;    }
          for ( int j = 0; j < to_[i].isize( ); j++ )
          {    to_[i][j] += nvertices;
               to_edge_obj_[i][j] += nedges;    }    }    }

template<class E>
void digraphE<E>::SplitEdge( int v, int j, const E& e1, const E& e2 )
{    int n = N( );
     int ne = edges_.size( );
     edges_.push_back( e1, e2 );
     int w = from_[v][j];
     int we = from_edge_obj_[v][j];
     int i = InputFromOutputTo( v, j );
     from_[v].erase( from_[v].begin( ) + j );
     from_edge_obj_[v].erase( from_edge_obj_[v].begin( ) + j );
     to_[w].erase( to_[w].begin( ) + i );
     to_edge_obj_[w].erase( to_edge_obj_[w].begin( ) + i );
     from_[v].push_back(n), from_edge_obj_[v].push_back(ne);
     vec<int> nfrom, nto;
     vec<int> nfrom_edge_obj, nto_edge_obj;
     nfrom.push_back(w), nfrom_edge_obj.push_back(ne+1);
     nto.push_back(v), nto_edge_obj.push_back(ne);
     from_.push_back(nfrom), to_.push_back(nto);
     from_edge_obj_.push_back(nfrom_edge_obj);
     to_edge_obj_.push_back(nto_edge_obj);
     for ( int u = 0; u < to_[w].isize( ); u++ )
     {    if ( to_[w][u] == v && we == to_edge_obj_[w][u] )
          {    to_.erase( to_.begin( ) + u );
               to_edge_obj_.erase( to_edge_obj_.begin( ) + u );
               break;    }    }
     to_[w].push_back(n), to_edge_obj_[w].push_back(ne+1);    }

template<class E> void digraphE<E>::Glue( const EmbeddedSubPath<E>& a,
     const EmbeddedSubPath<E>& b, const vec<int>& EE, const vec<int>& FF, 
     const digraphE<E>& c )
{    
     // Sanity check.

     ForceAssertGe( a.NVertices( ), 2 ); ForceAssertGe( b.NVertices( ), 2 );
     ForceAssert( !HasSharedEdge(a, b) );
     ForceAssertEq( EE.isize( ), a.NVertices( ) ); 
     ForceAssertEq( FF.isize( ), b.NVertices( ) );
     vec<int> Esort = EE, Fsort = FF;
     Sort(Esort), Sort(Fsort);
     ForceAssert( Esort.UniqueOrdered( ) );
     ForceAssert( Fsort.UniqueOrdered( ) );
     ForceAssertEq( EE.front( ), 0 ); ForceAssertEq( EE.back( ), c.N( ) - 1 );
     ForceAssertEq( FF.front( ), 0 ); ForceAssertEq( FF.back( ), c.N( ) - 1 );

     // Delete edges appearing in a and b.

     for ( int i = 0; i < a.NVertices( ) - 1; i++ )
     {    int v = a.Vertex(i), w = a.Vertex(i+1);
          int e = a.EdgeObjectIndexAbs(i);
          int ef = EdgeObjectIndexToFromIndex( v, e );
          int et = InputFromOutputTo( v, ef );
          from_[v].erase( from_[v].begin( ) + ef );
          from_edge_obj_[v].erase( from_edge_obj_[v].begin( ) + ef );
          to_[w].erase( to_[w].begin( ) + et );
          to_edge_obj_[w].erase( to_edge_obj_[w].begin( ) + et );    }
     for ( int i = 0; i < b.NVertices( ) - 1; i++ )
     {    int v = b.Vertex(i), w = b.Vertex(i+1);
          int e = b.EdgeObjectIndexAbs(i);
          int ef = EdgeObjectIndexToFromIndex( v, e );
          int et = InputFromOutputTo( v, ef );
          from_[v].erase( from_[v].begin( ) + ef );
          from_edge_obj_[v].erase( from_edge_obj_[v].begin( ) + ef );
          to_[w].erase( to_[w].begin( ) + et );
          to_edge_obj_[w].erase( to_edge_obj_[w].begin( ) + et );    }

     // Attach c.

     int nvertices = N( );
     Append(c);
     for ( int i = 0; i < a.NVertices( ); i++ )
          TransferEdges( a.Vertex(i), EE[i] + nvertices );
     for ( int i = 0; i < b.NVertices( ); i++ )
          TransferEdges( b.Vertex(i), FF[i] + nvertices );    

     // If data implies that some vertices in c should be identified, do so.

     vec< vec<int> > sources( c.N( ) );
     for ( int i = 0; i < a.NVertices( ); i++ )
          sources[ EE[i] ].push_back( a.Vertex(i) );
     for ( int i = 0; i < b.NVertices( ); i++ )
          sources[ FF[i] ].push_back( b.Vertex(i) );
     for ( int i = 0; i < c.N( ); i++ )
          Sort( sources[i] );
     for ( int i1 = 0; i1 < c.N( ); i1++ )
     {    for ( int i2 = i1 + 1; i2 < c.N( ); i2++ )
          {    if ( Meet( sources[i1], sources[i2] ) )
                    TransferEdges( i1 + nvertices, i2 + nvertices );    }    }    }

template<class E> void digraphE<E>::TransferEdges( int v, int w, 
     const Bool enter_only )
{    ForceAssert( v != w );

     // Change edges v --> v to edges w --> w.

     if ( !enter_only )
     {
     vec<Bool> remove_from_v;
     remove_from_v.resize_and_set( from_[v].size( ), False );
     for ( int i = 0; i < from_[v].isize( ); i++ )
     {    if ( from_[v][i] == v )
          {    from_[w].push_back(w);
               from_edge_obj_[w].push_back( from_edge_obj_[v][i] );
               to_[w].push_back(w);
               to_edge_obj_[w].push_back( from_edge_obj_[v][i] );
               remove_from_v[i] = True;
               int j = InputFromOutputTo( v, i );
               to_[v].erase( to_[v].begin( ) + j );
               to_edge_obj_[v].erase( to_edge_obj_[v].begin( ) + j );    }    }
     EraseIf( from_[v], remove_from_v );
     EraseIf( from_edge_obj_[v], remove_from_v );
     SortSync( from_[w], from_edge_obj_[w] );
     SortSync( to_[w], to_edge_obj_[w] );
     }

     // Change edges u --> v to edges u --> w.
     
     for ( int i = 0; i < to_[v].isize( ); i++ )
     {    int u = to_[v][i];
          int j = InputToOutputFrom( v, i );
          from_[u][j] = w;
          SortSync( from_[u], from_edge_obj_[u] );    }

     // Change edges v --> x to edges w --> x.

     // if ( !enter_only )
     {    for ( int i = 0; i < from_[v].isize( ); i++ )
          {    int x = from_[v][i];
               int j = InputFromOutputTo( v, i );
               if ( !enter_only ) to_[x][j] = w;
               else to_[x][j] = v;
               SortSync( to_[x], to_edge_obj_[x] );    }    }

     // Do the rest.

     if ( !enter_only )
     {    from_[w].append( from_[v] );
          from_edge_obj_[w].append( from_edge_obj_[v] );    }
     SortSync( from_[w], from_edge_obj_[w] );
     to_[w].append( to_[v] );
     to_edge_obj_[w].append( to_edge_obj_[v] );
     SortSync( to_[w], to_edge_obj_[w] );
     to_[v].clear( ), to_edge_obj_[v].clear( );    
     if ( !enter_only ) { from_[v].clear( ), from_edge_obj_[v].clear( ); }    }

#define VALID_FAILS(REASON)                \
{    cout << "\nDigraph is invalid.\n";    \
     cout << REASON << "\n";               \
     cout << "Abort.\n";                   \
     TracebackThisProcess( );    }

template<class E> void digraphE<E>::RemoveDuplicateEdges( )
{    for ( int v = 0; v < N( ); v++ )
     {    vec<Bool> remove;
          remove.resize_and_set( from_[v].size( ), False );
          for ( int j = 0; j < from_[v].isize( ); j++ )
          {    int k;
               for ( k = j + 1; k < from_[v].isize( ); k++ )
                    if ( from_[v][k] != from_[v][j] ) break;
               for ( int u1 = j; u1 < k; u1++ )
               {    if ( remove[u1] ) continue;
                    for ( int u2 = u1 + 1; u2 < k; u2++ )
                    {    if ( edges_[ from_edge_obj_[v][u1] ]
                                   == edges_[ from_edge_obj_[v][u2] ] )
                         {    remove[u2] = True;    }    }    }
               j = k - 1;    }
          for ( int i = 0; i < remove.isize( ); i++ )
          {    if ( remove[i] )
               {    int w = from_[v][i];
	            int j = InputFromOutputTo( v, i );
                    to_[w].erase( to_[w].begin( ) + j );
                    to_edge_obj_[w].erase( to_edge_obj_[w].begin( ) + j );    }    }
          EraseIf( from_[v], remove );
          EraseIf( from_edge_obj_[v], remove );    }    }

template<class E> void digraphE<E>::DeleteEdgesAtVertex( int v )
{    for ( int i = 0; i < from_[v].isize( ); i++ )
     {    int w = from_[v][i];
          int j = InputFromOutputTo( v, i );
          if ( v == w ) continue;
          to_[w].erase( to_[w].begin( ) + j );
          to_edge_obj_[w].erase( to_edge_obj_[w].begin( ) + j );    }
     for ( int i = 0; i < to_[v].isize( ); i++ )
     {    int w = to_[v][i];
          int j = InputToOutputFrom( v, i );
          if ( v == w ) continue;
          from_[w].erase( from_[w].begin( ) + j );
          from_edge_obj_[w].erase( from_edge_obj_[w].begin( ) + j );    }
     from_[v].clear( ), from_edge_obj_[v].clear( );
     to_[v].clear( ), to_edge_obj_[v].clear( );    }
                    
template<class E> void digraphE<E>::RemoveDeadEdgeObjects( )
{    vec<Bool> used;
     Used(used);
     int count = 0;
     vec<int> to_new_id( edges_.size( ) );
     for ( int i = 0; i < edges_.isize( ); i++ )
     {    if ( used[i] )
          {    if ( count != i ) edges_[count] = edges_[i];
               to_new_id[i] = count;
               ++count;    }    }
     edges_.resize(count);
     for ( int v = 0; v < N( ); v++ )
     {    for ( int i = 0; i < from_[v].isize( ); i++ )
               from_edge_obj_[v][i] = to_new_id[ from_edge_obj_[v][i] ];
          for ( int i = 0; i < to_[v].isize( ); i++ )
               to_edge_obj_[v][i] = to_new_id[ to_edge_obj_[v][i] ];    }    }

template<class E> void digraphE<E>::TestValid( ) const
{    digraph(*this).TestValid( );
     if ( from_edge_obj_.size( ) != to_edge_obj_.size( ) )
          VALID_FAILS( "sizes of from_edge_obj_ and to_edge_obj_ are different" );
     if ( from_.size( ) != from_edge_obj_.size( ) )
          VALID_FAILS( "sizes of from_ and from_edge_obj_ are different" );
     for ( int v = 0; v < N( ); v++ )
     {    if ( from_[v].size( ) != from_edge_obj_[v].size( ) )
          {    VALID_FAILS( "sizes of from_[" << v << "] and "
                    << "from_edge_obj_[" << v << "] are different" );    }    }
     for ( int v = 0; v < N( ); v++ )
     {    if ( to_[v].size( ) != to_edge_obj_[v].size( ) )
          {    VALID_FAILS( "sizes of to_[" << v << "] and "
                    << "to_edge_obj_[" << v << "] are different" );    }    }
     for ( int v = 0; v < N( ); v++ )
     {    for ( int j = 0; j < from_[v].isize( ); j++ )
          {    int w = from_[v][j];
               int ei = from_edge_obj_[v][j];
               Bool found = False;
               for ( int r = 0; r < to_[w].isize( ); r++ )
                    if ( to_[w][r] == v && to_edge_obj_[w][r] == ei ) found = True;
               if ( !found )
               {    VALID_FAILS( "There is an edge from " << v << " to " << w
                         << " in from_[" << v 
                         << "], but not in to_[" << w << "]." );    }    }    }
     for ( int v = 0; v < N( ); v++ )
     {    for ( int j = 0; j < to_[v].isize( ); j++ )
          {    int w = to_[v][j];
               int ei = to_edge_obj_[v][j];
               Bool found = False;
               for ( int r = 0; r < from_[w].isize( ); r++ )
               {    if ( from_[w][r] == v && from_edge_obj_[w][r] == ei ) 
                         found = True;    }
               if ( !found )
               {    VALID_FAILS( "There is an edge from " << v << " to " << w
                         << " in to_[" << v << "], but not in from_[" << w 
                         << "]." );    }    }    }    }

template<class E>
void digraphE<E>::Initialize( const vec< vec<int> >& from, const vec< vec<int> >& to,
          const vec<E>& edges, const vec< vec<int> >& to_edge_obj,
          const vec< vec<int> >& from_edge_obj )
{    digraph::Initialize( from, to );
     edges_ = edges; 
     to_edge_obj_ = to_edge_obj;
     from_edge_obj_ = from_edge_obj;
     int N = from.size( );
     ForceAssertEq( N, to_edge_obj.isize( ) );
     ForceAssertEq( N, from_edge_obj.isize( ) );
     vec<int> used( edges.size( ), 0 );
     for ( int i = 0; i < N; i++ ) 
     {    ForceAssertEq( to_edge_obj[i].size( ), to[i].size( ) );
          ForceAssertEq( from_edge_obj[i].size( ), from[i].size( ) );
          for ( int j = 0; j < to_edge_obj[i].isize( ); j++ )
          {    int o = to_edge_obj[i][j];
               ForceAssertGe( o, 0 );
               ForceAssertLt( o, edges.isize( ) );
               ++used[o];
               int w = i, v = to_[i][j];
               int wf = BinPosition( from[v], w );
               // The following assert won't do what we want if there are multiple
               // edges between two given vertices (in which case wf doesn't
               // make sense).
               // ForceAssertEq( o, from_edge_obj[v][wf] );    
                    }    }
     for ( int i = 0; i < used.isize( ); i++ )
     {    if ( used[i] != 1 )
          {    cout << "Edge " << i << " is used " << used[i]
                    << " times, whereas it should be used exactly once.\n";    }
          ForceAssertEq( used[i], 1 );    }    }

template<class E>
digraphE<E>::digraphE( const vec< vec<int> >& from, const vec< vec<int> >& to,
          const vec<E>& edges, const vec< vec<int> >& to_edge_obj,
          const vec< vec<int> >& from_edge_obj ) 
      : digraph(from, to) // redundant with initialize?
{    Initialize( from, to, edges, to_edge_obj, from_edge_obj );    }

template<class V>
void digraphV<V>::Initialize( const vec< vec<int> >& from, const vec< vec<int> >& to,
          const vec<V>& verts )
{    digraph::Initialize( from, to );
     verts_ = verts;
     ForceAssertEq( N( ), verts.isize( ) );    }

template<class V>
digraphV<V>::digraphV( const vec< vec<int> >& from, const vec< vec<int> >& to,
          const vec<V>& verts ) 
     : digraph(from, to) // redundant with initialize?
{    Initialize( from, to, verts );    }

template<class V, class E>
void digraphVE<V,E>::Initialize( const vec< vec<int> >& from, 
     const vec< vec<int> >& to, const vec<V>& verts, const vec<E>& edges, 
     const vec< vec<int> >& to_edge_obj, const vec< vec<int> >& from_edge_obj )
{    digraphE<E>::Initialize( from, to, edges, to_edge_obj, from_edge_obj );
     verts_ = verts;
     ForceAssertEq( from.size( ), verts.size( ) );    }

template<class V, class E>
digraphVE<V,E>::digraphVE( const vec< vec<int> >& from, 
     const vec< vec<int> >& to, const vec<V>& verts, const vec<E>& edges, 
     const vec< vec<int> >& to_edge_obj, const vec< vec<int> >& from_edge_obj )
     : digraphE<E>( from, to, edges, to_edge_obj, from_edge_obj ) // redundant??
{    Initialize( from, to, verts, edges, to_edge_obj, from_edge_obj );    }

template<class V, class E>
digraphVE<V,E>::digraphVE( const digraphE<E>& G, const vec<V>& verts )
     : digraphE<E>(G)
{    verts_ = verts;
     ForceAssertEq( G.N( ), verts.isize( ) );    }

template<class E> Bool digraphE<E>::IsComplete( 
     const vec<int>& vertices, const vec<int>& edges ) const
{    ForceAssert( vertices.UniqueOrdered( ) );
     ForceAssert( edges.UniqueOrdered( ) );
     for ( int u = 0; u < vertices.isize( ); u++ )
     {    int v = vertices[u];
          for ( int j = 0; j < From(v).isize( ); j++ )
          {    int w = From(v)[j];
               if ( !BinMember( vertices, w ) ) return False;
               int e = EdgeObjectIndexByIndexFrom( v, j );
               if ( !BinMember( edges, e ) ) return False;    }
          for ( int j = 0; j < To(v).isize( ); j++ )
          {    int w = To(v)[j];
               if ( !BinMember( vertices, w ) ) return False;
               int e = EdgeObjectIndexByIndexTo( v, j );
               if ( !BinMember( edges, e ) ) return False;    }    }
     return True;    }

template<class E> void digraphE<E>::DualComponentRelation( 
     equiv_rel& e, const vec<Bool>& exclude ) const
{    e.Initialize( EdgeObjectCount( ) );
     for ( int v = 0; v < N( ); v++ )
     {    for ( int j1 = 0; j1 < To(v).isize( ); j1++ )
          {    int e1 = EdgeObjectIndexByIndexTo( v, j1 );
               if ( exclude.nonempty( ) && exclude[e1] ) continue;
               for ( int j2 = 0; j2 < From(v).isize( ); j2++ )
               {    int e2 = EdgeObjectIndexByIndexFrom( v, j2 );
                    if ( exclude.nonempty( ) && exclude[e2] ) continue;
                    e.Join(e1, e2);    }     }    }    }

template<class E> void digraphE<E>::Initialize( 
     const digraphE& g, const equiv_rel& e )
{    edges_ = g.edges_;
     vec<int> reps;
     e.OrbitRepsAlt(reps);
     int nreps = reps.size( );
     vec<int> to_reps( g.N( ) );
     for ( int i = 0; i < nreps; i++ )
     {    vec<int> o;
          e.Orbit( reps[i], o );
          for ( int j = 0; j < o.isize( ); j++ )
               to_reps[ o[j] ] = i;    }
     from_.resize(nreps), to_.resize(nreps);
     from_edge_obj_.resize(nreps), to_edge_obj_.resize(nreps);
     int nedges = g.EdgeObjectCount( );
     vec<int> to_left_vertex(nedges, -1), to_right_vertex(nedges, -1);
     for ( int w = 0; w < g.N( ); w++ )
     {    for ( int j = 0; j < g.To(w).isize( ); j++ )
          {    int m = g.EdgeObjectIndexByIndexTo( w, j );
               int v = g.To(w)[j];
               to_left_vertex[m] = v, to_right_vertex[m] = w;    }    }
     for ( int m = 0; m < nedges; m++ )
     {    if ( to_left_vertex[m] < 0 || to_right_vertex[m] < 0 ) continue;
          int v = to_reps[ to_left_vertex[m] ];
          int w = to_reps[ to_right_vertex[m] ];
          from_[v].push_back(w), to_[w].push_back(v);
          from_edge_obj_[v].push_back(m), to_edge_obj_[w].push_back(m);    }
     for ( int v = 0; v < N( ); v++ )
     {    SortSync( from_[v], from_edge_obj_[v] );
          SortSync( to_[v], to_edge_obj_[v] );    }    }

template<class E> digraphE<E>::digraphE( 
     const digraphE& g, const equiv_rel& e ) : edges_( g.edges_ )
{    Initialize( g, e );    }

template<class E> digraphE<E>::digraphE( const vec<digraphE>& g )
{    Initialize(g);    }

template<class E> void digraphE<E>::Initialize( const vec<digraphE>& g )
{    for ( int i = 0; i < g.isize( ); i++ )
          Append( g[i] );    }

template<class E> void digraphE<E>::Initialize( const vec<digraphE>& g,
     const vec< pair< pair<int,int>, pair<int,int> > >& joins )
{    digraphE<E> G(g);
     equiv_rel e( G.N( ) );
     vec<int> start( g.isize( ) );
     start[0] = 0;
     for ( int i = 1; i < g.isize( ); i++ )
          start[i] = start[i-1] + g[i-1].N( );
     for ( int i = 0; i < joins.isize( ); i++ )
     {    int v = start[ joins[i].first.first ] + joins[i].first.second;
          int w = start[ joins[i].second.first ] + joins[i].second.second;
          e.Join( v, w );    }
     Initialize( G, e );    }

template<class E> digraphE<E>::digraphE( const vec<digraphE>& g,
     const vec< pair< pair<int,int>, pair<int,int> > >& joins )
{    Initialize( g, joins );    }


template<class E> void digraphE<E>::Initialize( const digraph& g, const vec<E>& edges ){    
  int nedges = g.N();
  ForceAssertEq( nedges, edges.isize() );
  equiv_rel e( 2*nedges );
  for ( int v = 0; v < nedges; v++ ){    
    for ( size_t j = 0; j < g.From(v).size( ); j++ ){    
      int w = g.From(v)[j];
      e.Join( 2*v + 1, 2*w );    
    }    
  }
  vec<int> reps;
  e.OrbitRepsAlt(reps);
 
  int N = reps.size( );
  vec< vec<int> > from(N), to(N);
  vec< vec<int> > from_edge_obj(N), to_edge_obj(N);
  for ( int i = 0; i < edges.isize( ); i++ ){    
    int x = BinPosition( reps, e.ClassId( 2*i ) );
    int y = BinPosition( reps, e.ClassId( 2*i + 1 ) );
    from[x].push_back(y), to[y].push_back(x);    
    from_edge_obj[x].push_back(i), to_edge_obj[y].push_back(i);    
  }
  
  for ( int i = 0; i < N; i++ ){
    SortSync( from[i], from_edge_obj[i] );
    SortSync( to[i], to_edge_obj[i] );
  }    

  Initialize( from, to, edges, to_edge_obj, from_edge_obj );
}

template<class E> digraphE<E>::digraphE( const digraph& g, const vec<E>& edges ){
  Initialize( g, edges );
}


template<class E> digraphE<E>::digraphE( const digraph& g ){
  vec<int> edges( g.N(), vec<int>::IDENTITY );
  Initialize( g, edges );
}


template<class E> Bool digraphE<E>::ThisClose( int v, int w, E d ) const
{    if ( d < 0 ) return False;
     if ( v == w ) return True;
     set< pair<int,E> > unprocessed, processed;

     unprocessed.insert( make_pair( v, 0 ) );
     while( !unprocessed.empty( ) )
     {    int x = unprocessed.begin( )->first;
          E dx = unprocessed.begin( )->second;
          typename set< pair<int,E> >::iterator u 
               = processed.lower_bound( make_pair( x, 0 ) );
          unprocessed.erase( unprocessed.begin( ) );
          if ( u != processed.end( ) && u->first == x )
          {    if ( u->second <= dx ) continue;
               processed.erase(u);    }
          processed.insert( make_pair( x, dx ) );
          for ( int j = 0; j < From(x).isize( ); j++ )
          {    int y = From(x)[j];
               E dy = dx + EdgeObjectByIndexFrom( x, j );
               if ( dy > d ) continue;
               if ( y == w ) return True;
               typename set< pair<int,E> >::iterator p 
                    = processed.lower_bound( make_pair( y, 0 ) );
               if ( p != processed.end( ) && p->first == y )
               {    if ( p->second <= dy ) continue;
                    processed.erase(p);    }
               typename set< pair<int,E> >::iterator u 
                    = unprocessed.lower_bound( make_pair( y, 0 ) );
               if ( u != unprocessed.end( ) && u->first == y )
               {    if ( u->second <= dy ) continue;
                    unprocessed.erase(u);    }
               unprocessed.insert( make_pair( y, dy ) );    }    }
     return False;    }

template<class E> void digraphE<E>::ToLeft( vec<int>& to_left ) const
{    to_left.resize( EdgeObjectCount( ) );
     for ( int i = 0; i < N( ); i++ )
     {    for ( int j = 0; j < From(i).isize( ); j++ )
          {    int e = EdgeObjectIndexByIndexFrom( i, j );
               to_left[e] = i;    }    }    }

template<class E> void digraphE<E>::ToRight( vec<int>& to_right ) const
{    to_right.resize( EdgeObjectCount( ) );
     for ( int i = 0; i < N( ); i++ )
     {    for ( int j = 0; j < To(i).isize( ); j++ )
          {    int e = EdgeObjectIndexByIndexTo( i, j );
               to_right[e] = i;    }    }    }

template<class E>
void digraphE<E>::GetSuccessors( const vec<int>& v, vec< pair<int,E> >& from_v )
{    set< pair<int,E> > check, fromv;

     for ( int i = 0; i < v.isize( ); i++ )
          check.insert( make_pair( v[i], 0 ) );
     while( !check.empty( ) )
     {    int x = check.begin( )->first;
          E dx = check.begin( )->second;
          typename set< pair<int,E> >::iterator u 
               = fromv.lower_bound( make_pair( x, 0 ) );
          check.erase( check.begin( ) );
          if ( u != fromv.end( ) && u->first == x )
          {    if ( u->second <= x ) continue;
               fromv.erase(u);    }
          fromv.insert( make_pair( x, dx ) );
          for ( int i = 0; i < From(x).isize( ); i++ )
          {    int y = From(x)[i]; 
               E dy = dx + EdgeObjectByIndexFrom( x, i );
               typename set< pair<int,E> >::iterator a 
                    = check.lower_bound( make_pair( y, 0 ) );
               if ( a != check.end( ) && a->first == y )
               {    if ( a->second <= dy ) continue;
                    check.erase(a);    }
               typename set< pair<int,E> >::iterator b 
                    = fromv.lower_bound( make_pair( y, 0 ) );
               if ( b != fromv.end( ) && b->first == y )
               {    if ( b->second <= dy ) continue;
                    fromv.erase(b);    }
               check.insert( make_pair( y, dy ) );    }    }
     from_v.clear( );
     for ( typename set< pair<int,E> >::iterator i = fromv.begin( ); 
          i != fromv.end( ); ++i )
     {    from_v.push_back(*i);    }    }

template<class E>
void
digraphE<E>::PrettyDOT( ostream& out, const vec<double>& lengths,
			Bool label_contigs, Bool label_vertices,
			Bool label_edges,
			const vec<int>* componentsToPrint,
                        const Bool edge_labels_base_alpha,
			const vec<String> *label_edges_extra,
			const vec<String> *label_contigs_extra,
			const vec<int> *verticesToPrint ) const
{
  // Set up output.
  out << "digraph G {\n\n";
  // out << "orientation=landscape;\n";
  if (label_vertices)
    out << "node [width=0.1,height=0.1,fontsize=12,shape=plaintext];\n";
  else out << "node [width=0.1,height=0.1,fontsize=10,shape=point];\n";
  out << "edge [fontsize=12];\n";
  if (label_contigs) out << "margin=1.0;\n";
  out << "rankdir=LR;\n";
  out << "labeljust=l;\n";
  
  // Define components.
  vec< vec<int> > components;
  Components( components );
  
  // Contig labels, and estimate space after label.
  vec<int> label_distance;
  vec<String> strContigLabel;
  if ( label_contigs ) {
    if ( label_contigs_extra )
      strContigLabel = *label_contigs_extra;
    else {
      strContigLabel.resize( components.size( ) );
      for (size_t ii=0; ii<components.size( ); ii++) {
	strContigLabel[ii] = "contig " + ToString( ii );
      }
    }
    
    label_distance.resize( strContigLabel.size( ), 0 );
    for (int ii=0; ii<(int)strContigLabel.size( ); ii++)
      label_distance[ii] = 1 + (int)( strContigLabel[ii].size( ) / 2 );
  }
  
  // Selected components.
  vec<int> select;
  if ( componentsToPrint ) select = *componentsToPrint;
  else {
    select.reserve( components.size( ) );
    for (int ii=0; ii<(int)components.size( ); ii++) select.push_back( ii );
  }

  // Vertices to be skipped.
  vec<bool> skip_vtx;
  if ( verticesToPrint ) {
    skip_vtx.resize( this->N( ), true );
    for (size_t ii=0; ii<verticesToPrint->size( ); ii++)
      skip_vtx[ (*verticesToPrint)[ii] ] = false;
  }
  
  // Print the contigs.  We put each contig in its own cluster (the
  // subgraph's name MUST start with "cluster" for this to have any effect).
  for ( int sel_id = select.isize( ) - 1; sel_id >= 0; sel_id-- ) {
    int i = select[sel_id];

    out << "\nsubgraph cluster" << i << " {\n";
    out << "color=white;\n";
    if ( label_contigs && label_contigs_extra )
      out << "label=\"" << strContigLabel[i]
	  << "\","
	  << "fontsize=18,"
	  << "fontname=\"Times-Bold\"\n";
    
    vec<int> &o = components[i];
    
    // Find "leftmost" vertex.
    Sort(o);
    vec<float> pos( o.size( ) );
    vec<Bool> placed( o.size( ), False );
    pos[0] = 0.0, placed[0] = True;
    while( Sum(placed) < o.isize( ) ) {
      for ( int i1 = 0; i1 < o.isize( ); i1++ ) {
	int v = o[i1];
	for ( int j = 0; j < From(v).isize( ); j++ ) {
	  int w = From(v)[j];
	  int i2 = BinPosition( o, w );
	  if ( !( placed[i1] ^ placed[i2] ) ) continue;
	  edge_t e = EdgeObjectIndexByIndexFrom( v, j );
	  if ( placed[i1] ) pos[i2] = pos[i1] + lengths[e];
	  else pos[i1] = pos[i2] - lengths[e];
	  placed[i1] = placed[i2] = True;
	}
      }
    }
    
    float left = Min(pos);
    int leftj = 0;
    for ( leftj = 0; leftj < pos.isize( ); leftj++ )
      if ( pos[leftj] == left ) break;
    int leftv = o[leftj];
    
    // Print component.
    
    for ( int vi = 0; vi < o.isize( ); vi++ ) {
      int v = o[vi];
      if ( verticesToPrint && skip_vtx[v] ) continue;
      if (label_vertices)
	out << v << " [label=" << "\"" << v << "\"" 
	    << ",fontcolor=black];\n";
      
      for ( int j = 0; j < From(v).isize( ); j++ ) {
	int ei = EdgeObjectIndexByIndexFrom( v, j );
	int w = From(v)[j];
	float wd = 0.1; // this value not used
	String color, label;
	Bool bold = False;
	double len = lengths[ei];
	if ( len < 100.0 ) {
	  color = "gray";
	  if ( v == w ) label = ToString( len, 0 );
	  wd = 1.0;
	}
	else if ( len >= 100.0 && len < 1000.0 ) {
	  color = "black";
	  wd = 2.0;
	}
	else if ( len >= 1000.0 && len < 10000.0 ) {
	  color = "red";
	  wd = 4.0;
	  label = ToString( len/1000.0, 1 ) + " kb";
	}
	else {
	  color = "magenta";
	  bold = True;
	  wd = 8.0;
	  label = ToString( len/1000.0, 0 ) + " kb";
	}
	out << v << " -> " << w
	    << " [minlen=" << wd << ",color=" << color;
	if (bold) out << ",style=bold";
	if (label_edges) {
	  if ( label == "" ) 
	    label = ( edge_labels_base_alpha ? BaseAlpha(ei) : ToString(ei) );
	  else
	    label = ( edge_labels_base_alpha ? BaseAlpha(ei) : ToString(ei) ) 
	      + " (" + label + ")"; }
        if ( label_edges_extra ) label += " " + (*label_edges_extra)[ei];
	
	if ( label != "" ) out << ",label=\"" << label << "\"";
	if ( label_contigs && v == leftv && j == 0 && ! label_contigs_extra )
	  out << ",taillabel=\"" << strContigLabel[i] 
	      << "\",labelangle=180,"
	      << "weight=10000,"
	      << "labeldistance=" << label_distance[i] << ",labelfontsize=18,"
	      << "labelfontname=\"Times-Bold\"";
	out << "];\n";
      }
    }
    out << "}\n";
  }
  
  out << "\n}" << endl;
}

// Method: DumpGraphML
// Output the digraph structure in a textual format that can be easily
// read without reference to our code base.
template<class E>
void
digraphE<E>::DumpGraphML( const String& graphMLFileName ) const
{
  vec< vec< String > > edgeLabels( N() );
  for ( int v = 0; v < N( ); v++ ) {
    for ( int j = 0; j < From(v).isize( ); j++ ) {
      int w = From(v)[j];
      edgeLabels[ v ].push_back( BaseAlpha( EdgeObjectIndexByIndexFrom( v, j ) ) );
    }
  }
  
  Ofstream( grml, graphMLFileName );
  WriteGraphML( grml, edgeLabels );
}

template<class E> void digraphE<E>::ComponentsE( vec< vec<int> >& comp ) const
{    comp.clear( );
     equiv_rel e( N( ) );
     for ( int v = 0; v < N( ); v++ )
     {    for ( int j = 0; j < From(v).isize( ); j++ )
               e.Join( v, From(v)[j] );    }
     for ( int x = 0; x < N( ); x++ )
     {    if ( e.Representative(x) )
          {    vec<int> o;
               e.Orbit( x, o );
               Sort(o);
               vec<int> C;
               for ( int i = 0; i < o.isize( ); i++ )
               {    int v = o[i];
                    for ( int j = 0; j < From(v).isize( ); j++ )
                         C.push_back( EdgeObjectIndexByIndexFrom( v, j ) );    }
               comp.push_back(C);    }    }    }

template<class E> void LongestPath( const digraphE<E>& G, int (E::*len)( ) const,
     vec<int>& a_longest_path )
{    vec<int> D;
     const int infinity = 2000000000;
     DistancesToEnd( G, len, infinity, True, D );
     int M = 0, v = 0;
     for ( int x = 0; x < G.N( ); x++ )
          if ( D[x] > M ) { M = D[x], v = x; }
     a_longest_path.clear( );
     while( G.From(v).nonempty( ) )
     {    for ( int j = 0; j < G.From(v).isize( ); j++ )
          {    int w = G.From(v)[j];
               if ( D[w] == D[v] - ((G.EdgeObjectByIndexFrom( v, j )).*len)( ) )
               {    a_longest_path.push_back( G.EdgeObjectIndexByIndexFrom( v, j ) );
                    v = w;
                    break;    }    }    }    }

template<class E> void DistancesToEnd( const digraphE<E>& G,
     int (E::*len)( ) const, const int max_dist, const Bool fw, vec<int>& D )
{
     // Let D(v) be the maximum length of a path starting at v, to be computed.
     // Define initial values for D(v) to be 'undefined', except for sinks, 
     // which are zero.

     D.resize_and_set( G.N( ), -1 );
     for ( int v = 0; v < G.N( ); v++ )
     {    if ( fw && G.Sink(v) ) D[v] = 0;
          if ( !fw && G.Source(v) ) D[v] = 0;    }

     // Initialize vertices to process.

     vec<Bool> to_process( G.N( ), False );
     vec<int> to_processx;
     for ( int v = 0; v < G.N( ); v++ )
     {    if ( (fw && G.Sink(v)) || (!fw && G.Source(v)) ) 
          {    to_process[v] = True, to_processx.push_back(v);    }    }

     // Now compute D.  Uncomputed values are set to 'infinity'.

     while( to_processx.nonempty( ) )
     {    int v = to_processx.back( );
          to_processx.pop_back( );
          to_process[v] = False;
          for ( int j = 0; j < (fw ? G.To(v) : G.From(v) ).isize( ); j++ )
          {    int w = ( fw ? G.To(v) : G.From(v) )[j];
               if ( D[w] >= max_dist ) continue;
               const E& e = ( fw ? G.EdgeObjectByIndexTo(v, j)
                    : G.EdgeObjectByIndexFrom(v, j) );
               int Dw_new = (e.*len)( ) + D[v];
               if ( Dw_new > D[w] ) 
               {    D[w] = Dw_new;
                    if ( !to_process[w] )
                    {    to_process[w] = True; 
                         to_processx.push_back(w);    }    }    }    }
     for ( int v = 0; v < G.N( ); v++ )
          if ( D[v] < 0 ) D[v] = max_dist;    }

template<class E> void RemoveHangingEnds( digraphE<E>& G, 
     int (E::*len)( ) const, const int max_del, const double min_ratio )
{
     // Track hanging ends.

     vec<Bool> hanging( G.EdgeObjectCount( ), False );

     // Define the maximum length that we care about.

     const int max_dist = int( ceil( double(max_del) * min_ratio ) );

     // Go through two passes (forward and reverse).

     for ( int pass = 1; pass <= 2; pass++ )
     {
          // Compute distances to end.

          vec<int> D;
          DistancesToEnd( G, len, max_dist, pass == 1, D );

          // Identify hanging ends.

          for ( int v = 0; v < G.N( ); v++ )
          {    const vec<int>& V = ( pass == 1 ? G.From(v) : G.To(v) );
               vec<int> d( V.size( ) ); 
               vec<int> id( V.size( ), vec<int>::IDENTITY );
               for ( int j = 0; j < V.isize( ); j++ )
               {    d[j] = ((pass == 1 
                         ? G.EdgeObjectByIndexFrom(v,j) : G.EdgeObjectByIndexTo(v,j))
                         .*len)( ) + D[ V[j] ];    }
               ReverseSortSync( d, id );
               for ( int j = 1; j < d.isize( ); j++ )
               {    if ( d[j] <= max_del && d[0] >= d[j] * min_ratio )
                    {    hanging[ ( pass == 1 
                              ? G.EdgeObjectIndexByIndexFrom( v, id[j] )
                              : G.EdgeObjectIndexByIndexTo( v, id[j] ) ) ] 
                              = True;    }    }    }    }

     // Remove hanging ends.

     vec<int> to_delete;
     for ( int i = 0; i < G.EdgeObjectCount( ); i++ )
          if ( hanging[i] ) to_delete.push_back(i);
     G.DeleteEdges(to_delete);    }


// Remove short hanging ends.  Look for
//
//                 x
//                 |
//                 e
//                 |
//        u --c--> v --d--> w
//
// where x is a source or sink, e is short (and can go either way), whereas
// c and d are long.  Works for T = HyperKmerPath and T = HyperFastavector.

template<class T> void RemoveHangingEnds2( T& h ) {

  for ( int x = 0; x < h.N( ); x++ ) {
    
    // Check that basic assumptions are satisfied, including length(e) <= 5kb.
    
    int v, c, d, e;
    if ( h.Source(x) && h.From(x).size( ) == 1 ) {
      v = h.From(x)[0];
      e = h.EdgeObjectIndexByIndexFrom( x, 0 );
    } else if ( h.Sink(x) && h.To(x).size( ) == 1 ) {
      v = h.To(x)[0];
      e = h.EdgeObjectIndexByIndexTo( x, 0 );
    } else 
      continue;

    if ( h.EdgeLengthKmers(e) > 5000 ) continue;

    if ( h.Source(x) ) {
      if ( !( h.From(v).size( ) == 1 && h.To(v).size( ) == 2 ) ) continue;
      d = h.EdgeObjectIndexByIndexFrom( v, 0 );
      c = h.EdgeObjectIndexByIndexTo( v, 0 );
      if ( c == e ) c = h.EdgeObjectIndexByIndexTo( v, 1 );
    } else {
      if ( !( h.From(v).size( ) == 2 && h.To(v).size( ) == 1 ) ) continue;
      c = h.EdgeObjectIndexByIndexTo( v, 0 );
      d = h.EdgeObjectIndexByIndexFrom( v, 0 );
      if ( d == e ) d = h.EdgeObjectIndexByIndexFrom( v, 1 );
    }

    // We require that there is an edge "competing with e", that is at least
    // 20 times longer.
    
    static vec<int> v_only(1), to_v, from_v;
    v_only[0] = v;
    int max_competitor = 0;
    if ( h.Source(x) ) {
      h.digraph::GetPredecessors( v_only, to_v );
      for ( int j = 0; j < to_v.isize( ); j++ ) {
	int z = to_v[j];
	for ( int i = 0; i < h.To(z).isize( ); i++ ) {
	  int e = h.EdgeObjectIndexByIndexTo( z, i );
	  max_competitor = Max( max_competitor, h.EdgeLengthKmers(e) );
	}
      }
    } else {
      h.digraph::GetSuccessors( v_only, from_v );
      for ( int j = 0; j < from_v.isize( ); j++ ) {
	int z = from_v[j];
	for ( int i = 0; i < h.From(z).isize( ); i++ ) {
	  int e = h.EdgeObjectIndexByIndexFrom( z, i );
	  max_competitor = Max( max_competitor, h.EdgeLengthKmers(e) );
	}
      }
    }

    if ( 20 * h.EdgeLengthKmers(e) > max_competitor ) continue;
    
    // Edit the graph.
    
    if ( h.Source(x) ) h.DeleteEdgeFrom( x, 0 );
    else h.DeleteEdgeTo( x, 0 );

  }
}


// Find the indices of all edges e that form self-loops, i.e., e goes from v -> v.
template<class E> vec<int> digraphE<E>::SelfLoops( ) const
{
  vec<int> to_left, to_right;
  ToLeft( to_left );
  ToRight( to_right );
  vec<Bool> used;
  Used( used );
  
  vec<int> self_loops;
  for ( int i = 0; i < EdgeObjectCount(); i++ )
    if ( to_left[i] == to_right[i] && used[i] )
      self_loops.push_back( i );
  
  return self_loops;
}




template<class E> void digraphE<E>::LoopSubgraph( vec<int>& loop_edges ) const
{    loop_edges.clear( );
     vec< vec<int> > SCC;
     StronglyConnectedComponents(SCC);
     for ( int i = 0; i < SCC.isize( ); i++ )
     {    const vec<int>& V = SCC[i];
          for ( int r = 0; r < V.isize( ); r++ )
          {    int v = V[r];
               for ( int j = 0; j < From(v).isize( ); j++ )
               {    if ( BinMember( V, From(v)[j] ) )
                    {    loop_edges.push_back( EdgeObjectIndexByIndexFrom( 
                              v, j ) );    }    }    }    }
     Sort(loop_edges);    }

template<class E> void digraphE<E>::SplayVertex( const int v )
{    int n = N( );
     AddVertices( To(v).size( ) );
     for ( int j = To(v).isize( ) - 1; j >= 0; j-- )
          GiveEdgeNewToVx( EdgeObjectIndexByIndexTo( v, j ), v, n + j );
     n = N( );
     AddVertices( From(v).size( ) );
     for ( int j = From(v).isize( ) - 1; j >= 0; j-- )
     {    GiveEdgeNewFromVx( EdgeObjectIndexByIndexFrom( v, j ), 
               v, n + j );    }    }

template<class E> void digraphE<E>::LiberateEdge( 
     const int e, const int v, const int w )
{    int j = EdgeObjectIndexToFromIndex( v, e );
     DeleteEdgeFrom( v, j );
     SplayVertex(v), SplayVertex(w);    }

template<class E> void digraphE<E>::GiveEdgeNewFromVx
( int edge_id, int old_from_v, int new_from_v ) {
       int i = Position( from_edge_obj_[old_from_v], edge_id );
       ForceAssert( i != -1 );
       int w = from_[old_from_v][i];
       int j = Position( to_edge_obj_[w],edge_id );
       ForceAssert( j != -1 );
       to_[w][j] = new_from_v;
       from_[old_from_v].erase( from_[old_from_v].begin() + i );
       from_edge_obj_[old_from_v].erase( from_edge_obj_[old_from_v].begin() + i );
       from_[new_from_v].push_back(w);
       from_edge_obj_[new_from_v].push_back(edge_id);
       SortSync( to_[w], to_edge_obj_[w] );
       SortSync( from_[new_from_v], from_edge_obj_[new_from_v] );
     }

template<class E> void digraphE<E>::GiveEdgeNewToVx
( int edge_id, int old_to_w, int new_to_w ) {
       int j = Position( to_edge_obj_[old_to_w], edge_id );
       ForceAssert( j != -1 );
       int v = to_[old_to_w][j];
       int i = Position( from_edge_obj_[v],edge_id );
       ForceAssert( i != -1 );
       from_[v][i] = new_to_w;
       to_[old_to_w].erase( to_[old_to_w].begin() + j );
       to_edge_obj_[old_to_w].erase( to_edge_obj_[old_to_w].begin() + j );
       to_[new_to_w].push_back(v);
       to_edge_obj_[new_to_w].push_back(edge_id);
       SortSync( from_[v], from_edge_obj_[v] );
       SortSync( to_[new_to_w], to_edge_obj_[new_to_w] );
     }

template<class E> void digraphE<E>::AddEdge( int v, int w, const E& e )
{    int n = EdgeObjectCount( );
     edges_.push_back(e);
     int i = upper_bound( from_[v].begin(), from_[v].end(), w ) - from_[v].begin();
     from_[v].insert( from_[v].begin()+i, w );
     from_edge_obj_[v].insert( from_edge_obj_[v].begin()+i, n );
     int j = upper_bound( to_[w].begin(), to_[w].end(), v ) - to_[w].begin();
     to_[w].insert( to_[w].begin()+j, v );
     to_edge_obj_[w].insert( to_edge_obj_[w].begin()+j, n );    }

template<class E>
Bool digraphE<E>::EdgePaths( const int v, const int w, vec< vec<int> >& paths,
     const int max_copies, const int max_paths, const int max_iterations )
{    vec<int> left, right;
     ToLeft(left), ToRight(right);

     // Pretest to determine if the computation will explode.  This only works if
     // max_copies is not set.

     if ( max_copies < 0 && ( max_paths >= 0 || max_iterations >= 0 ) )
     {    vec<int> subs;
          int path_count = 0;
          for ( int i = 0; i < From(v).isize( ); i++ )
          {    int e = EdgeObjectIndexByIndexFrom( v, i );
               subs.push_back(e);    }
          int iterations = 0;
          while( subs.nonempty( ) )
          {    if ( max_iterations > 0 && ++iterations > max_iterations ) 
                    return False;
               int p = subs.back( );
               subs.pop_back( );
               int x = right[p];
               if ( x == w ) 
               {    if ( max_paths >= 0 && ++path_count > max_paths ) 
                         return False;    }
               else
               {    for ( int j = 0; j < From(x).isize( ); j++ )
                    {    int e = EdgeObjectIndexByIndexFrom( x, j );
                         subs.push_back(e);    }    }    }    }

     // Now do the computation for real.

     vec< vec<int> > subs;
     paths.clear( );
     for ( int i = 0; i < From(v).isize( ); i++ )
     {    int e = EdgeObjectIndexByIndexFrom( v, i );
          vec<int> one;
          one.push_back(e);
          subs.push_back(one);    }
     int iterations = 0;
     while( subs.nonempty( ) )
     {    if ( max_iterations > 0 && ++iterations > max_iterations ) return False;
          vec<int> p = subs.back( );
          subs.resize( subs.isize( ) - 1 );
          int x = right[ p.back( ) ];
          if ( x == w ) 
          {    paths.push_back(p);
               if ( max_paths >= 0 && paths.isize( ) > max_paths ) return False;    }
          else
          {    for ( int j = 0; j < From(x).isize( ); j++ )
               {    int e = EdgeObjectIndexByIndexFrom( x, j );
                    vec<int> pp(p);
                    pp.push_back(e);
                    if ( max_copies >= 0 )
                    {    vec<int> pps(pp);
                         Sort(pps);
                         Bool fail = False;
                         for ( int r = 0; r < pps.isize( ); r++ )
                         {    int s = pps.NextDiff(r);
                              if ( s - r > max_copies )
                              {    fail = True;
                                   break;    }
                              r = s - 1;    }
                         if (fail) continue;    }
                    subs.push_back(pp);    }    }    }
     return True;    }

template<class E> void digraphE<E>::DeleteEdgeTo( int w, int j )
     {    int v = to_[w][j];
          int i = InputToOutputFrom( w, j );
          to_[w].erase( to_[w].begin( ) + j );
          to_edge_obj_[w].erase( to_edge_obj_[w].begin( ) + j );
          from_[v].erase( from_[v].begin( ) + i );
          from_edge_obj_[v].erase( from_edge_obj_[v].begin( ) + i );    }

template<class E> void digraphE<E>::DeleteEdgeFrom( int v, int j )
     {    int w = from_[v][j];
          int i = InputFromOutputTo( v, j );
          from_[v].erase( from_[v].begin( ) + j );
          from_edge_obj_[v].erase( from_edge_obj_[v].begin( ) + j );
          to_[w].erase( to_[w].begin( ) + i );
          to_edge_obj_[w].erase( to_edge_obj_[w].begin( ) + i );    }

template<class E> 
vec<int> digraphE<E>::EdgesBetween( const int v, const int w ) const
{    vec<int> b;
     for ( int i = 0; i < From(v).isize( ); i++ )
     {    if ( From(v)[i] == w )
               b.push_back( EdgeObjectIndexByIndexFrom( v, i ) );    }
     return b;    }

template<class E> vec<int> digraphE<E>::EdgesBetween( const vec<int>& v ) const
{    vec<int> b;
     for ( int j = 0; j < v.isize( ); j++ )
     {    for ( int i = 0; i < From(v[j]).isize( ); i++ )
          {    if ( BinMember( v, From(v[j])[i] ) )
                    b.push_back( EdgeObjectIndexByIndexFrom( v[j], i ) );    }    }
     Sort(b);
     return b;    }

template<class E> 
vec<E> digraphE<E>::EdgeObjectsBetween( const int v, const int w ) const
{    vec<E> b;
     for ( int i = 0; i < From(v).isize( ); i++ )
     {    if ( From(v)[i] == w )
               b.push_back( EdgeObjectByIndexFrom( v, i ) );    }
     return b;    }

template<class E> int digraphE<E>::InputToOutputFrom( int w, int i ) const
{    int v = to_[w][i];
     int ei = to_edge_obj_[w][i];
     for ( int j = 0; j < from_[v].isize( ); j++ )
          if ( from_edge_obj_[v][j] == ei ) return j;
     ForceAssert( 0 == 1 );
     return -1;    }

template<class E> int digraphE<E>::InputFromOutputTo( int w, int i ) const
{    int v = from_[w][i];
     int ei = from_edge_obj_[w][i];
     for ( int j = 0; j < to_[v].isize( ); j++ )
          if ( to_edge_obj_[v][j] == ei ) return j;
     ForceAssert( 0 == 1 );
     return -1;    }

template<class E> void digraphE<E>::ChangeEdgeObjectFrom( int v, int i, const E& e )
{    int ne = edges_.size( );
     edges_.push_back(e);
     int w = From(v)[i];
     int j = InputFromOutputTo( v, i );
     from_edge_obj_[v][i] = ne;
     to_edge_obj_[w][j] = ne;    }

template<class E> E digraphE<E>::MinEdge( int v, int w )
{    E m = 0;
     Bool first = True;
     for ( int j = 0; j < From(v).isize( ); j++ )
     {    if ( From(v)[j] != w ) continue;
          if (first) m = EdgeObjectByIndexFrom( v, j );
          else m = Min( m, EdgeObjectByIndexFrom( v, j ) );
          first = False;    }
     ForceAssert( !first );
     return m;    }

template<class E> E digraphE<E>::MaxEdge( int v, int w )
{    E M = 0;
     Bool first = True;
     for ( int j = 0; j < From(v).isize( ); j++ )
     {    if ( From(v)[j] != w ) continue;
          if (first) M = EdgeObjectByIndexFrom( v, j );
          else M = Max( M, EdgeObjectByIndexFrom( v, j ) );
          first = False;    }
     ForceAssert( !first );
     return M;    }

template<class E> void digraphE<E>::AddVertices( int nadd ) 
{    int nvert = N( );
     from_.resize( nvert + nadd );
     to_.resize( nvert + nadd );
     from_edge_obj_.resize( nvert + nadd );
     to_edge_obj_.resize( nvert + nadd );    }

template<class E> void digraphE<E>::DeleteEdges( const vec<int>& to_delete )
{    vec<int> to_delete_local;
     if ( !to_delete.UniqueOrdered( ) )
     {    to_delete_local = to_delete;
          UniqueSort(to_delete_local);    }
     const vec<int>& tod 
          = ( to_delete_local.nonempty( ) ? to_delete_local: to_delete );
     for ( int v = 0; v < N( ); v++ )
     {    for ( int j = From(v).isize( ) - 1; j >= 0; j-- )
          {    int e = EdgeObjectIndexByIndexFrom( v, j );
               if ( BinMember( tod, e ) ) DeleteEdgeFrom( v, j );    }    }    }

template<class E> void digraphE<E>::DeleteEdges( const vec<int>& to_delete,
     const vec<int>& to_left )
{    vec<int> to_delete_local;
     if ( !to_delete.UniqueOrdered( ) )
     {    to_delete_local = to_delete;
          UniqueSort(to_delete_local);    }
     const vec<int>& tod 
          = ( to_delete_local.nonempty( ) ? to_delete_local: to_delete );
     vec<int> vs;
     for ( int i = 0; i < to_delete.isize( ); i++ )
          vs.push_back( to_left[ to_delete[i] ] );
     UniqueSort(vs);
     for ( int i = 0; i < vs.isize( ); i++ )
     {    int v = vs[i];
          for ( int j = From(v).isize( ) - 1; j >= 0; j-- )
          {    int e = EdgeObjectIndexByIndexFrom( v, j );
               if ( BinMember( tod, e ) ) DeleteEdgeFrom( v, j );    }    }    }

template<class E> int digraphE<E>::EdgeObjectIndexByIndexTo( int v, int j ) const
{    CheckGoodVertex(v);
     AssertGe( j, 0 );
     AssertLt( j, to_edge_obj_[v].isize( ) );
     return to_edge_obj_[v][j];    }
   
template<class E> int digraphE<E>::EdgeObjectIndexToFromIndex( int v, int e ) const
{    AssertGe( v, 0 );
     AssertLt( v, from_edge_obj_.isize( ) );
     for ( int i = 0; i < from_edge_obj_[v].isize( ); i++ )
          if ( from_edge_obj_[v][i] == e ) return i;
     return -1;    }

template<class E> int digraphE<E>::EdgeObjectIndexToToIndex( int v, int e ) const
{    AssertGe( v, 0 );
     AssertLt( v, to_edge_obj_.isize( ) );
     for ( int i = 0; i < to_edge_obj_[v].isize( ); i++ )
          if ( to_edge_obj_[v][i] == e ) return i;
     return -1;    }

template<class E> bool operator!=( const digraphE<E>& g1, const digraphE<E>& g2 )
{ return !(g1==g2); }

template<class E> bool operator==( const digraphE<E>& g1, const digraphE<E>& g2 )
{
    if ( static_cast<digraph const&>(g1) != static_cast<digraph const&>(g2) )
        return false;

    // digraphs are the same, now check edge objects
    typedef vec<int> V;
    typedef V::const_iterator VI;
    typedef vec<V> VV;
    typedef VV::const_iterator VVI;
    VV const& vv1 = g1.FromEdgeObj();
    VV const& vv2 = g2.FromEdgeObj();
    if ( vv1.size() != vv2.size() )
        return false;

    VVI oE(vv1.end());
    for ( VVI o1(vv1.begin()), o2(vv2.begin()); o1 != oE; ++o1, ++o2 )
    {
        if ( o1->size() != o2->size() )
            return false;

        VI iE(o1->end());
        for ( VI i1(o1->begin()), i2(o2->begin()); i1 != iE; ++i1, ++i2 )
            if ( !(g1.EdgeObject(*i1) == g2.EdgeObject(*i2)) )
                return false;
    }
    return true;
}

template<class E> 
     void Compare( ostream& out, const digraphE<E>& g1, const digraphE<E>& g2 )
{    if ( g1.N( ) != g2.N( ) )
          cout << "first graph has " << g1.N( ) << " vertices but "
               << "second graph has " << g2.N( ) << "\n";
     if ( g1.From( ) != g2.From( ) ) cout << "from_ not the same\n";
     if ( g1.To( ) != g2.To( ) ) cout << "to_ not the same\n";
     if ( g1.Edges( ) != g2.Edges( ) ) cout << "edges_ not the same\n";
     if ( g1.ToEdgeObj( ) != g2.ToEdgeObj( ) )
          cout << "to_edge_obj_ not the same\n";
     if ( g1.FromEdgeObj( ) != g2.FromEdgeObj( ) )
          cout << "from_edge_obj_ not the same\n";
     if ( g1 != g2 ) cout << "DIGRAPHS ARE NOT EQUAL\n";
     return;    }

template<class E> void digraphE<E>::Clear( )
{    from_.clear( ), to_.clear( );
     from_edge_obj_.clear( ), to_edge_obj_.clear( );
     edges_.clear( );    }

template<class E> const E& digraphE<E>::EdgeObject( int i ) const
{    AssertGe( i, 0 );
     AssertLt( i, edges_.isize( ) );
     return edges_[i];    }

template<class E> E& digraphE<E>::EdgeObjectMutable( int i )
{    AssertGe( i, 0 );
     AssertLt( i, edges_.isize( ) );
     return edges_[i];    }

template<class V> const V& digraphV<V>::Vert( int v ) const
{    AssertGe( v, 0 );
     AssertLt( v, N( ) );
     return verts_[v];    }

template<class V> V& digraphV<V>::VertMutable( int v )
{    AssertGe( v, 0 );
     AssertLt( v, N( ) );
     return verts_[v];    }

template<class V, class E> const V& digraphVE<V,E>::Vert( int v ) const
{    AssertGe( v, 0 );
     AssertLt( v, N( ) );
     return verts_[v];    }

template<class V, class E> V& digraphVE<V,E>::VertMutable( int v )
{    AssertGe( v, 0 );
     AssertLt( v, N( ) );
     return verts_[v];    }

template<class V> void digraphV<V>::DeleteVertex( const int v )
{    int n = N( );
     AssertGe( v, 0 );
     AssertLt( v, n );
     DeleteEdgesAtVertex(v);
     verts_.erase( verts_.begin( ) + v );
     from_.erase( from_.begin( ) + v );
     to_.erase( to_.begin( ) + v );
     for ( int x = 0; x < n - 1; x++ )
     {    for ( int j = 0; j < From(x).isize( ); j++ )
               if ( From(x)[j] >= v ) FromMutable(x)[j]--;
          for ( int j = 0; j < To(x).isize( ); j++ )
               if ( To(x)[j] >= v ) ToMutable(x)[j]--;    }    }

template<class V> void digraphV<V>::DeleteVertices( const vec<int>& v )
{    for ( int m = v.isize( ) - 1; m >= 0; m-- )
          DeleteVertex( v[m] );    }

template<class V> void digraphV<V>::AddVertex( const V& v )
{    verts_.push_back(v);
     from_.resize( from_.size( ) + 1 );
     to_.resize( to_.size( ) + 1 );    }

template<class E> 
vec<int> digraphE<E>::EdgesSomewhereBetween( const int v, const int w ) const
{    vec<int> answer, after_v, before_w, both;
     GetSuccessors1( v, after_v ), GetPredecessors1( w, before_w );
     Intersection( after_v, before_w, both );
     for ( int l = 0; l < both.isize( ); l++ )
     {    int s = both[l];
          for ( int j = 0; j < From(s).isize( ); j++ )
          {    int t = From(s)[j];
               if ( BinMember( both, t ) ) 
                    answer.append( EdgesBetween( s, t ) );    }    }
     UniqueSort(answer);
     return answer;    }

template<class E>
size_t digraphE<E>::writeBinary( BinaryWriter& writer ) const
{
    size_t len = digraph::writeBinary(writer);

    typedef vec< vec<int> >::const_iterator OItr;
    typedef vec<int>::const_iterator Itr;
    OItr oend(from_edge_obj_.end());
    typename vec<E>::size_type totSize = 0;
    for ( OItr oitr(from_edge_obj_.begin()); oitr != oend; ++oitr )
        totSize += oitr->size();
    len += writer.write(totSize);
    for ( OItr oitr(from_edge_obj_.begin()); oitr != oend; ++oitr )
    {
        vec<int> const& fff = *oitr;
        for ( Itr itr(fff.begin()), end(fff.end()); itr != end; ++itr )
            len += writer.write(edges_[*itr]);
    }

    return len;
}

template<class E>
void digraphE<E>::readBinary( BinaryReader& reader )
{
    digraph::readBinary(reader);
    reader.read(&edges_);

    Mimic(from_,from_edge_obj_);
    typedef vec< vec<int> >::iterator OItr;
    typedef vec<int>::iterator Itr;
    int count = 0;
    OItr oend(from_edge_obj_.end());
    for ( OItr oitr(from_edge_obj_.begin()); oitr != oend; ++oitr )
    {
        vec<int>& fff = *oitr;
        for ( Itr itr(fff.begin()), end(fff.end()); itr != end; ++itr )
            *itr = count++;
    }

    MimicReserve(to_,to_edge_obj_);
    for ( size_t iii = 0; iii < from_.size(); ++iii )
    {
        vec<int> const& fff = from_[iii];
        vec<int> const& fe = from_edge_obj_[iii];
        for ( size_t jjj = 0; jjj < fe.size(); ++jjj )
            to_edge_obj_[fff[jjj]].push_back(fe[jjj]);
    }
}

template<class V, class E>
size_t digraphVE<V,E>::writeBinary( BinaryWriter& writer ) const
{    size_t len = digraphE<E>::writeBinary(writer);
     len += writer.write(verts_);
     return len;    }

template<class V, class E>
void digraphVE<V,E>::readBinary( BinaryReader& reader )
{    digraphE<E>::readBinary(reader);
     reader.read( &verts_ );    }

template<class E> void EmbeddedSubPath<E>::TestValid( ) const
{    ForceAssertEq( e_.isize( ), a_.isize( ) - 1 );
     for ( int u = 0; u < a_.isize( ) - 1; u++ )
     {    const vec<int>& fr = D_->From( a_[u] );
          ForceAssertGe( e_[u], 0 );
          ForceAssertLt( e_[u], fr.isize( ) );
          ForceAssertEq( fr[ e_[u] ], a_[u+1] );
          ForceAssertEq( D_->EdgeObjectIndexByIndexFrom( a_[u], e_[u] ),
               esafe_[u] );    }    }

template<class E> void EmbeddedSubPath<E>::Repair( )
{    for ( int u = 0; u < e_.isize( ); u++ )
     {    if ( D_->EdgeObjectIndexByIndexFrom( a_[u], e_[u] ) != esafe_[u] )
               e_[u] = D_->EdgeObjectIndexToFromIndex(
                    a_[u], esafe_[u] );    }    }

#endif

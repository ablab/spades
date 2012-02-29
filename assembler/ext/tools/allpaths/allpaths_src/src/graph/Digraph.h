///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// =================================================================================
// PLEASE READ: THERE ARE OTHER FILES YOU MAY WANT TO LOOK AT:
// - DigraphPaths.h
// - FindCells.h
// Please add to this list if you add new files having digraph functions.
// =================================================================================

// This file contains the definition of class digraph, which represents a finite
// directed graph, digraphV<V>, which has vertex objects from class V, and
// digraphE<E>, which has edge objects from class E.  Ultimately we may also
// want a class digraphVE<V,E>.
//
// The file Digraph.cc contains some of the definitions for digraph functions.
//
// The file DigraphTemplate.h contains some of the templatized definitions for
// digraphE functions.  Putting these in a separate file makes this file more
// readable, but has the additional advantage of reducing the inlining of these
// functions, thereby reducing compilation time and executable size (in principle).
// Generally, the less places DigraphTemplate.h is included, the better.

#ifndef DIGRAPH_H
#define DIGRAPH_H

#include "Bitvector.h"
#include "CoreTools.h"
#include "Equiv.h"
#include "VecUtilities.h"
#include "math/Functions.h"
#include "math/Matrix.h"
#include "system/TraceVal.h"
#include "feudal/BinaryStream.h"
#include <cstddef>

typedef pair<int, int> VertexPair;

typedef int vrtx_t;
typedef int edge_t;

/**
   Class: digraph

   A digraph may have multiple edges between two given vertices, but no edges
   from a vertex to itself.

   A digraph is stored as two vec< vec<int> >'s, "to" and "from".  For each 
   vertex v, to[v] is the *sorted* list of vertices which have an edge to v, and 
   from[v] is the *sorted* list of vertices which have an edge from v.
*/
// TODO: potentially dangerous truncation of index by this class's use of ints

class digraph
{
     public:

     digraph( ) { }
     digraph( const vec< vec<int> >& from, const vec< vec<int> >& to );
     digraph( const matrix< Bool >& );
     explicit digraph( int n ) { Initialize( n ); }

     // "Reducing" constructors
     // Reduce g to the vertices in v
     digraph( const digraph& g, const vec<int>& v );
     // Reduce g to the vertex/vertices in component #i
     digraph( const digraph& g, const int i );
     
     void Initialize( const vec< vec<int> >& from, const vec< vec<int> >& to );
     void Initialize( const digraph& g, const vec<int>& v );
     void Initialize( const matrix< Bool >& );

     // create an empty graph with n nodes

     void Initialize( int n ) {
       Clear();
       from_.resize( n );
       to_.resize( n );
     }

     void TestValid( ) const; // incomplete test

     void Clear( )
     {    from_.clear( ), to_.clear( );    }

     int N( ) const { return from_.size( ); } // number of vertices

     void CheckGoodVertex( int v ) const
     {    AssertGe( v, 0 );
          AssertLt( v, N( ) );    }

     const vec< vec<int> >& From( ) const { return from_; }
     const vec< vec<int> >& To( ) const { return to_; }
     vec< vec<int> >& FromMutable( ) { return from_; }
     vec< vec<int> >& ToMutable( ) { return to_; }

     const vec<int>& From( int v ) const { CheckGoodVertex(v); return from_[v]; }
     const vec<int>& To( int v ) const { CheckGoodVertex(v); return to_[v]; }
     int FromSize( int v ) const { CheckGoodVertex(v); return from_[v].size( ); }
     int ToSize( int v ) const { CheckGoodVertex(v); return to_[v].size( ); }

     vec<int>& FromMutable( int v ) { CheckGoodVertex(v); return from_[v]; }
     vec<int>& ToMutable( int v ) { CheckGoodVertex(v); return to_[v]; }

     // Create new vertices.

     void AddVertices( int nadd ) {
       int nvert = N( );
       from_.resize( nvert + nadd );
       to_  .resize( nvert + nadd );
     }

     // Add an edge from v to w.

     void AddEdge( int v, int w );
  
     Bool HasEdge( int v, int w ) const;

     Bool Source( int v ) const { return To(v).empty( ); }
     Bool Sink( int v ) const { return From(v).empty( ); }

     void Sources( vec<int>& v ) const;
     void Sinks( vec<int>& v ) const;

     // Given the complete subgraph defined by a set of vertices S, find its sources
     // and sinks.  The set S must be sorted.

     void SubgraphSources( const vec<int>& S, vec<int>& v ) const;
     void SubgraphSinks( const vec<int>& S, vec<int>& v ) const;

     // GetPredecessors: find all vertices which have a directed path to a vertex
     // in v.  This includes the vertices in v by definition.  Return a sorted list 
     // to_v.  GetSuccessors: go the other way.  
     // See also another version of GetSuccessors for class digraphE.

     void GetPredecessors( const vec<int>& v, vec<int>& to_v ) const;
     void GetSuccessors( const vec<int>& v, vec<int>& from_v ) const;
     void GetPredecessors1( const int v, vec<int>& to_v ) const;
     void GetSuccessors1( const int v, vec<int>& from_v ) const;

     // Reverse the graph.

     void Reverse( );

     // Components: find the connected components.  Each component is a sorted list
     // of vertices.  ComponentsAlt orders the components differently.
     
     void Components( vec< vec<int> >& comp ) const;
     void ComponentsAlt( vec< vec<int> >& comp ) const;
     size_t NComponents() const;

     // Determine if a graph has a path of nonzero length from v to v.

     Bool LoopAt( const int v ) const;

     // Return the strongly connected components of a graph.  The code is an
     // algorithm of Tarjan.  We borrowed pseudocode from Wikipedia.  The answer
     // is a sorted vector of sorted vectors.

     void StronglyConnectedComponents( vec< vec<int> >& SCC ) const;

     // Determine if the connected component defined by vertices sub has a directed
     // cycle.  Don't call this with anything other than a connected component.

     Bool HasCycle( const vec<int>& sub ) const;

     // Determine if graph is acyclic.

     Bool Acyclic( ) const;

     // ComponentRelation: return equivalence relation corresponding to connected 
     // components.  Note that this can be very slow on a large graph.  The function
     // "Components" can return the same information, faster.

     void ComponentRelation( equiv_rel& e ) const;

     // Return number of connected components in a graph.

     int ConnectedComponents( ) const;

     // Return the connected component containing a given vertex.

     vec<int> ComponentOf( const int v );

     // Return a sorted list of all vertices whose removal would increase the 
     // number of components in the graph.  Complexity = O( sum( n^2 log(n) ) ),
     // where n ranges over the component sizes.

     void CutPoints( vec<int>& cuts ) const;

     // Given an edge referred to by To, find the corresponding edge in From,
     // and the other way.

     int InputToOutputFrom( int w, int i ) const;
     int InputFromOutputTo( int w, int i ) const;

     // Delete an edge.  Indices of edges <j are untouched, while those >j decrease 
     // by 1.  In particular, adding an edge and then deleting it is guaranteed to 
     // leave the indices of other edges unchanged.

     void DeleteEdgeTo( int w, int j );
     void DeleteEdgeFrom( int v, int j );

     // Remove all edges entering or exiting a given vertex.

     void DeleteEdgesAtVertex( int v );

     // AllPaths: find all paths between two vertices which contain no duplicated
     // vertices.  Zero length paths (i.e. paths from a vertex to the same vertex)
     // are included.  Return answer as a list of lists, with each inner list a 
     // sequence of vertices.  This last can then be expanded out to reflect
     // alternative edges between adjacent vertices.
     //
     // If w = -1, instead find all paths from v to a sink.
     // If v = -1, instead find all paths from a source to w.
     //
     // If maxpaths >= 0, fail, returning False if greater than that many paths or 
     // partial paths found.
     //
     // If allow_self_loop and v = w (and >= 0), then loops from v to v are
     // allowed (but v cannot appear more than twice).
     //
     // Note that the behavior of this code may not be as expected.  Thus in this 
     // case
     //                    y        x
     //      * ------> * <----- * -----> *
     //                  -----> 
     //
     // AllPaths( -1, -1, ... ) will not use edge y.

     Bool AllPaths( int v, int w, vec< vec<int> >& paths, int maxpaths = -1,
          const Bool allow_self_loop = False, const int maxpushes = -1 ) const;

     friend void BinaryWrite( int fd, const digraph& g );
     friend void BinaryRead( int fd, digraph& g );

     // Create a representation of a given graph in DOT:
     // http://www.graphviz.org/doc/info/lang.html
     // The edge_colors and edge_labels arguments should be in bijective 
     // correspondence with from_.

     void DOT( ostream& out ) const;
     void DOT( ostream& out, const vec<String>& vertex_colors ) const;
     void DOT( ostream& out, const vec<String>& vertex_colors,
          const vec< vec<String> >& edge_colors ) const;
     void DOT( ostream& out, const vec< vec<String> >& edge_labels ) const;
     void PrettyDOT( ostream& out, Bool label_contigs = True, 
		     Bool label_vertices = False,
		     const vec<int>* componentsToPrint = NULL, 
		     const vec<String> *label_contigs_extra = NULL,
		     const vec<int> *verticesToPrint = NULL ) const;

     // DOT_vl: generate a graph with vertex labels in circles.  There is an 
     // optional argument layout, which could be set to circo if you want 
     // circular layout.  There is another optional argument legend, which is 
     // rendered inside a square box, but at a random place (not necessarily at
     // the top).  Lines within the box are left-right centered.  The version with
     // legends is similar: you can have more than one.

     void DOT_vl( ostream& out, const vec<String> & vertex_labels,
          const String& layout = "", 
          const vec<String>& legend = vec<String>( ), 
          const String& color = "" ) const;

     void DOT_vl( ostream& out, const vec<String> & vertex_labels,
          const String& layout = "", 
          const vec< vec<String> >& legends = vec< vec<String> >( ),
          const vec<String>& colors = vec<String>( ),
          const vec< vec<String> >& edge_labels = vec< vec<String> >( ) ) const;

     // Print the complete subgraph defined by the given vertices vs, which need
     // not be sorted.

     void DOT_vl( ostream& out, const vec<String> & vertex_labels, 
          vec<int> vs ) const;

     // Create a representation of a given graph in GraphML:
     // http://graphml.graphdrawing.org
     // The edge_labels arguments should be in bijective 
     // correspondence with from_.

     void WriteGraphML( ostream& out, const vec< vec<String> >& edge_labels ) const;

     /// ReplaceWithTransitiveClosure()
     /// Replace this graph with its transitive closure.  That is, if there is a
     /// directed path from u to w in the original graph, there is a direct edge
     /// (u,w) in the new graph.

     void ReplaceWithTransitiveClosure();

     size_t writeBinary( BinaryWriter& writer ) const
     { return writer.write(from_); }

     void readBinary( BinaryReader& reader );

     static size_t externalSizeof() { return 0; }

     friend bool operator==( digraph const& g1, digraph const& g2 )
     { return g1.from_ == g2.from_; }

     friend bool operator!=( digraph const& g1, digraph const& g2 )
     { return g1.from_ != g2.from_; }

     protected:

     vec< vec<int> > from_, to_;

};

SELF_SERIALIZABLE(digraph);

void BinaryWrite( int fd, const digraph& g );
void BinaryRead( int fd, digraph& g );

// =================================================================================

// Initial minimal implementation of digraphV.

template<class V> class digraphV : public digraph {

     public:

     void TestValid( ) const; // incomplete test

     digraphV( ) { }

     digraphV( const vec< vec<int> >& from, const vec< vec<int> >& to,
          const vec<V>& verts );
     void Initialize( const vec< vec<int> >& from, const vec< vec<int> >& to,
          const vec<V>& verts );

     const V& Vert( int v ) const;
     V& VertMutable( int v );

     void AddVertex( const V& v );

     // DeleteVertex: delete a vertex and all edges incident upon it.
     // Deleting a vertex renumbers everything.

     void DeleteVertex( const int v );

     // DeleteVertices take a unique-sorted vector as input.

     void DeleteVertices( const vec<int>& v );

     private:

     vec<V> verts_;

};

template <class V>
struct Serializability< digraphV<V> > : public SelfSerializable {};

// =================================================================================

template<class E> class EmbeddedSubPath; // forward declaration

template<class E> class digraphE;

template<class E> void BinaryWrite( int fd, const digraphE<E>& g );
template<class E> void BinaryRead( int fd, digraphE<E>& g );

/**
   Class Template: digraphE

   A digraphE<E> is stored as a digraph, plus a vec<E> edges_, plus indexing vectors
   (vec< vec<int> >'s to_edge_obj_ and from_edge_obj_) such that 
   if to_[w][i] = v, then to_edge_obj_[w][i] is the index in edges_ of the edge
   object for v --> w, and similarly for from_edge_obj_.

   The vec<E> edges_ is allowed to have unused entries (although not initially) and 
   there is a RemoveDeadEdgeObjects( ) function to eliminate them.

   Note that there are a bunch of member functions given here for digraphE<E>,
   that could also be implemented for class digraph.
*/

template<class E> class digraphE : public digraph {

     public:

     void TestValid( ) const; // incomplete test

     // Constructor 1: build the empty digraph.

     digraphE( ) { }

     // Constructor 2: build digraph from arbitrary digraph data.
     // (Note the INSANE order of the arguments: first from before to, then
     // to_edge_obj before from_edge_obj.)

     digraphE( const vec< vec<int> >& from, const vec< vec<int> >& to,
          const vec<E>& edges, const vec< vec<int> >& to_edge_obj,
          const vec< vec<int> >& from_edge_obj );
     void Initialize( const vec< vec<int> >& from, const vec< vec<int> >& to,
          const vec<E>& edges, const vec< vec<int> >& to_edge_obj,
          const vec< vec<int> >& from_edge_obj );

     // Constructor 3a: given a collection of edge objects, create a graph having
     // one edge and two vertices for each of the edge objects.
     // [ invoked by: digraphE( const vec<E>& edges, EDGES_SEPARATE ); ]
     // Constructor 3b: given a collection of edge objects, create a line graph out
     // of the edges.
     // [ invoked by: digraphE( const vec<E>& edges, EDGES_IN_LINE ); ]

     enum ConstructorBehavior {
       EDGES_SEPARATE,
       EDGES_IN_LINE 
     };

     digraphE( const vec<E>& edges, const ConstructorBehavior constructor_type  );

     // Constructor 4: given a collection of edge objects, and an equivalence 
     // relation on them, build a graph having two vertices per equivalence class,
     // with one edge between those two vertices for each member of the equivalence
     // class.

     digraphE( const vec<E>& edges, const equiv_rel& e );

     // Same as above, but not a constructor:

     void EdgeEquivConstructor( const vec<E>& edges, const equiv_rel& e );

     // Constructor 5: given a digraph, and given a list of vertex indices,
     // create the digraph having those vertices (with indices starting at 0,
     // but in the given order), and having all the edges that were between those
     // vertices.  Thus this is a "complete subgraph" constructor.  The order of the
     // edges is first by vertex v, and then by the order within From(v).

     digraphE( const digraphE& g, const vec<int>& v );
     void Initialize( const digraphE& g, const vec<int>& v );

     // Constructor 6: extract the nth connected component from another digraph.

     digraphE( const digraphE& g, int n );

     // Constructor 7: from another digraph and replacement edge objects.

     template<class F> digraphE( const digraphE<F>& g, const vec<E>& edges )
     {    from_ = g.From( );
          to_ = g.To( );
          from_edge_obj_ = g.FromEdgeObj( );
          to_edge_obj_ = g.ToEdgeObj( );
          edges_ = edges;    }
     template<class F> void Initialize ( const digraphE<F>& g, const vec<E>& edges )
     {    from_ = g.From( );
          to_ = g.To( );
          from_edge_obj_ = g.FromEdgeObj( );
          to_edge_obj_ = g.ToEdgeObj( );
          edges_ = edges;    }

     // Constructor 8: from a given digraph and an equivalence relation on the
     // vertices.  The vertices of the constructed digraph are the equivalence
     // classes of vertices of the given graph.
     
     digraphE( const digraphE& g, const equiv_rel& e );
     void Initialize( const digraphE& g, const equiv_rel& e );

     // Constructor 9: form the disjoint union of a collection of digraphs.

     digraphE( const vec<digraphE>& g );
     void Initialize( const vec<digraphE>& g );

     // Constructor 10: from a collection of digraphs and a set of identifications
     // between vertices in their disjoint union, each of which is specified as
     // ( (g1,v1), (g2,v2) ) where g1, g2 refer to graphs and v1, v2 refer to
     // vertices on those graphs.

     digraphE( const vec<digraphE>& g, 
          const vec< pair< pair<int,int>, pair<int,int> > >& joins );
     void Initialize( const vec<digraphE>& g, 
          const vec< pair< pair<int,int>, pair<int,int> > >& joins );

     // Constructor 11: from a given digraph and a collection of subsets of its
     // edges, each of which is given the induced subgraph structure, when are then
     // merged into a disjoint union.  The numbering of edge objects in the new
     // graph is the obvious order C[0][0], C[0][1], ..., C[1][0], C[1][1], ... .

     digraphE( const digraphE& g, const vec< vec<int> >& C );


     // Constructor 12: from a digraph and the collection of edges corresponding
     // to digraph vertices

     digraphE( const digraph& g, const vec<E>& edges );
     void Initialize( const digraph& g, const vec<E>& edges ); 
 
     // Constructor 13: from a digraph (edges will correspond to vertex numbers in digraph) 
     digraphE( const digraph& g );
     

     // Delete the entire graph.

     void Clear( );

     // Return number of edge objects.

     int EdgeObjectCount( ) const { return edges_.size( ); }

     // A bunch of ways to access edge objects:

     const E& EdgeObject( int i ) const;

     E& EdgeObjectMutable( int i );

     const E& EdgeObjectByIndexFrom( int v, int j ) const
     {    CheckGoodVertex(v);
          AssertGe( j, 0 );
          AssertLt( j, from_edge_obj_[v].isize( ) );
          return edges_[ from_edge_obj_[v][j] ];    }

     E& EdgeObjectByIndexFromMutable( int v, int j )
     {    CheckGoodVertex(v);
          AssertGe( j, 0 );
          AssertLt( j, from_edge_obj_[v].isize( ) );
          return edges_[ from_edge_obj_[v][j] ];    }

     int EdgeObjectIndexByIndexFrom( int v, int j ) const
     {    CheckGoodVertex(v);
          AssertGe( j, 0 );
          AssertLt( j, from_edge_obj_[v].isize( ) );
          return from_edge_obj_[v][j];    }

     const E& EdgeObjectByIndexTo( int v, int j ) const
     {    CheckGoodVertex(v);
          AssertGe( j, 0 );
          AssertLt( j, to_edge_obj_[v].isize( ) );
          return edges_[ to_edge_obj_[v][j] ];    }

     E& EdgeObjectByIndexToMutable( int v, int j )
     {    CheckGoodVertex(v);
          AssertGe( j, 0 );
          AssertLt( j, to_edge_obj_[v].isize( ) );
          return edges_[ to_edge_obj_[v][j] ];    }

     int EdgeObjectIndexByIndexTo( int v, int j ) const;

     // The following two functions return -1 if they can't find e.
   
     int EdgeObjectIndexToFromIndex( int v, int e ) const;
     int EdgeObjectIndexToToIndex( int v, int e ) const;

     const vec<int>& FromEdgeObj( int v ) const { return from_edge_obj_[v]; }
     const vec<int>& ToEdgeObj( int v ) const { return to_edge_obj_[v]; }
     vec<int>& FromEdgeObjMutable( int v ) { return from_edge_obj_[v]; }
     vec<int>& ToEdgeObjMutable( int v ) { return to_edge_obj_[v]; }

     const vec< vec<int> >& FromEdgeObj( ) const { return from_edge_obj_; }
     const vec< vec<int> >& ToEdgeObj( ) const { return to_edge_obj_; }
     vec< vec<int> >& FromEdgeObjMutable( ) { return from_edge_obj_; }
     vec< vec<int> >& ToEdgeObjMutable( ) { return to_edge_obj_; }

     const vec<E>& Edges( ) const { return edges_; }
     vec<E>& EdgesMutable( ) { return edges_; }

     // Find indices of all edges between two vertices.

     vec<int> EdgesBetween( const int v, const int w ) const;

     // Return sorted list of indices of all edges that are between two vertices, 
     // possibly indirectly.

     vec<int> EdgesSomewhereBetween( const int v, const int w ) const;

     // Find indices of all edges between a set v of vertices.  We assume that
     // v is sorted.

     vec<int> EdgesBetween( const vec<int>& v ) const;

     // Find all edges between two vertices.

     vec<E> EdgeObjectsBetween( const int v, const int w ) const;

     // Given an edge referred to by To, find the corresponding edge in From,
     // and the other way.

     int InputToOutputFrom( int w, int i ) const;
     int InputFromOutputTo( int w, int i ) const;

     // Change an edge object.

     void SetEdgeObject( int i, const E& e )
     {    edges_[i] = e;    }

     void ChangeEdgeObjectFrom( int v, int i, const E& e );

     // Move an edge's endpoint, without changing the edge object itself.

     void GiveEdgeNewFromVx( int edge_id, int old_from_v, int new_from_v );
     void GiveEdgeNewToVx( int edge_id, int old_to_w, int new_to_w );

     // MinEdge: find the minimum length of an edge from v to w.  Assert if there is
     // no edge.  This only makes sense if Min(E,E) is defined.

     E MinEdge( int v, int w );

     // MaxEdge: find the maximum length of an edge from v to w.  Assert if there is
     // no edge.  This only makes sense if Max(E,E) is defined.

     E MaxEdge( int v, int w );

     // GetSuccessors: find all vertices which have a directed path from a vertex
     // in v, and return the minimum distance to each, as a sorted list.

     void GetSuccessors( const vec<int>& v, vec< pair<int,E> >& from_v );

     // Return equivalence relation e on the edges generated by the rule that if 
     // edge A follows edge B, then they are equivalent.  The classes of e are thus
     // the components of the directed line graph associated to the given graph.
     // If "exclude" is specified, it should have one entry per edge.  The excluded
     // edges are not used in creating the equivalence relation.

     void DualComponentRelation( equiv_rel& e, const vec<Bool>& exclude ) const;

     // ToLeft, ToRight: create vectors that map edge indices to vertex indices.

     void ToLeft( vec<int>& to_left ) const;
     void ToRight( vec<int>& to_right ) const;

     // Use Dijkstra's algorithm or the Bellman-Ford algorithm to find the shortest 
     // path in a digraphE<int> between the 'start' vertex and the 'stop' vertex, 
     // where the length of the path is defined to be the sum of the edge values.  
     // Negative edge values are allowed but negative cycles will cause the code
     // to assert.  The handling of ties is not specified here.  Return the list 
     // of vertices that defines the path.  In the event of failure, return an 
     // empty path.
     //
     // Implementation stolen from pseudocode in Wikipedia.  Dijkstra's algorithm
     // is Implemented as O(N^2), which is stupid.  Should be reimplemented using 
     // priority queue.
     //
     // Note that there is a faster algorithm for acyclic graphs, see 
     // Cormen et al. p. 536.
     // 
     // Compare DistancesToEnd.

     void ShortestPath( const int start, const int stop, vec<int>& path ) const;

     // Add vertices.

     void AddVertices( int nadd );

     // Add an edge from v to w.  Indices of edges from v to any x<=w are untouched,
     //  while ones from v to x>w are incremented by 1; likewise for w's edges.

     void AddEdge( const int v, const int w, const E& e );

     // Delete an edge.  Indices of edges <j are untouched, while those >j decrease 
     // by 1.  In particular, adding an edge and then deleting it is guaranteed to 
     // leave the indices of other edges unchanged.

     void DeleteEdgeTo( int w, int j );
     void DeleteEdgeFrom( int v, int j );

     // Two functions to delete a bunch of edges.  For both functions, the input
     // list does not have to be sorted, and may contain duplicates.  The first
     // Delete function is O(N), where N is the number of vertices in the graph.
     // The second Delete function is O(to_delete).

     void DeleteEdges( const vec<int>& to_delete );
     void DeleteEdges( const vec<int>& to_delete, const vec<int>& to_left );
     
     // SplayVertex.  Replace a given vertex v by To(v).size( ) + From(v).size( )
     // vertices and move the edges entering and exiting v to these vertices.
     // This does not actually touch the edge objects.

     void SplayVertex( const int v );

     // LiberateEdge.  Remove a given edge e: v --> w, and spread out all the
     // edges that touch v or w, so that the two vertices v and w are replaced
     // by From(v) + To(v) + From(w) + To(w) - 2 edges.  (The - 2 is there to
     // account for e's contribution.)

     void LiberateEdge( const int e, const int v, const int w );

     // ComponentsE: find the connected components.  Each component is a list 
     // of edges.
     
     void ComponentsE( vec< vec<int> >& comp ) const;

     // ThisClose: determine if there is a path from v to w whose sum of edge
     // objects is <= d.  This assumes that edge objects can be added and that
     // <= makes sense on them.  In the special case v = w, return True (so long
     // as d >= 0).

     Bool ThisClose( int v, int w, E d ) const;

     // Find the indices of all edges e that form self-loops, i.e.,
     // e goes from v -> v.
     vec<int> SelfLoops( ) const;
     
     // Find the 'loop subgraph' of a given graph.  This consists of all edges
     // are part of a loop.  We return a sorted list of the edge indices.  

     void LoopSubgraph( vec<int>& loop_edges ) const;

     // Determine if the given vertices and edges (sorted lists) comprise the 
     // union of zero or more connected components.

     Bool IsComplete( const vec<int>& vertices, const vec<int>& edges ) const;

     // Remove all edges entering or exiting a given vertex.

     void DeleteEdgesAtVertex( int v );

     // If two edges have the same stop and start, and are same as E-objects,
     // delete one.  This does not actually delete the objects.

     void RemoveDuplicateEdges( );

     // Determine which edge objects are used.

     void Used( vec<Bool>& used ) const;

     // Eliminate unused edge objects.

     void RemoveDeadEdgeObjects( );

     // Remove vertices having no edges coming in or going out, or only specified
     // vertices having this same property.

     void RemoveEdgelessVertices( );
     void RemoveEdgelessVertices( const vec<int>& to_remove );

     // Reverse the graph or the connected component containing a given vertex of 
     // it.  The current version does nothing to the edge objects.

     void Reverse( );
     void ReverseComponent( int v );

     // Change order of vertices.

     void ReorderVertices( const vec<int>& new_order );

     // Change order of components.

     void ReorderComponents( const vec<int>& new_order );

     // Find the edges in each connected component.

     void ComponentEdges( vec< vec<int> >& edges ) const;

     // Distance: determine the lengths of all directed paths from vertex v to
     // vertex w, which are no longer than max_dist.  Return answer as D.
     // This is only implemented for the case where E = int.

     void Distance( int v, int w, int max_dist, vec<int>& D ) const;

     // EdgePaths: find all paths from vertex v to w, as lists of edges.

     Bool EdgePaths( const int v, const int w, vec< vec<int> >& paths,
          const int max_copies = -1, const int max_paths = -1,
          const int max_iterations = -1 );

     // AllPathsFixedLength: find all paths from v to w having length L.  This
     // is only implemented for the case where E = int.

     void AllPathsFixedLength( int v, int w, int L, 
          vec< vec<int> >& paths ) const;

     // AllPathsLengthRange: find all paths from v to w having length between L1
     // and L2, as lists of edges.  This is only implemented for the case where 
     // E = int.  If max_paths is specified and more than that many paths are found,
     // return False.  If max_loops is specified and the code loops more than that 
     // number of times, return False.

     Bool AllPathsLengthRange( int v, int w, int L1, int L2, 
          const vec<int>& to_right, vec< vec<int> >& paths, 
          int max_paths = 0, int max_loops = 0, const Bool no_dups = False ) const;

     // AllPathsLengthRangeAlt: same as AllPathsLengthRange, but processes as FIFO
     // instead of LIFO.

     Bool AllPathsLengthRangeAlt( int v, int w, int L1, int L2, 
          const vec<int>& to_right, vec< vec<int> >& paths, 
          int max_paths = 0, int max_loops = 0, const Bool no_dups = False,
          const Bool eq_ok = False, const int max_partials = 0 ) const;

     // SplitEdge( v, j, e1, e2 ): Replace edge from_[v][j] from v to w by a pair 
     // of edges 
     //
     //   e1    e2
     // v --> n --> w,
     //
     // where n is a new vertex [numbered N( ), if N( ) is called before SplitEdge].
     // This pushes e1 and then e2 onto edges_, and leaves an unused entry in 
     // edges_ where the original edge was.

     void SplitEdge( int v, int j, const E& e1, const E& e2 );

     // JoinEdges( x, e ): Suppose vertex x has exactly one edge entering it (from
     // a vertex v) and exactly one edge exiting it (to a vertex w, v != w).  
     // Remove both edges and substitute a single edge from v to w, leaving x as a 
     // vertex having no edges coming in or going out, and also leaving unused 
     // entries in edges_.

     void JoinEdges( int x, const E& e );

     // RemoveUnneededVertices: Apply JoinEdges to all applicable vertices x.  This 
     // assumes that class E has a member function Append( const E& d ), used here 
     // to compute the edge e in JoinEdges( x, e ).  This function removes edgeless 
     // vertices, but will leave unused entries in edges_.  If you want to remove
     // them, call RemoveDeadEdgeObjects.

     void RemoveUnneededVertices( );

     // PopBubbles.
     // Input is a set of vertices v.  Each v must be located at the opening of
     // a bubble, with exactly two edges that lead to the same successor w:
     //        _-_
     //  --> v     w -->
     //        -_-
     // Re-route all edges entering v so that they enter w instead.  Then
     // remove v and the bubble edges from the graph.
     // Last, contract any redundant edges by a call to RemoveUnneededVertices.
     void PopBubbles( const vec<int> & bubble_vs );

     // Append(D): append another digraph, forming the disjoint union.  The
     // new vertices from D are numbered after the existing vertices.

     void Append( const digraphE<E>& D );

     // TransferEdges( v, w ): move all edges entering and exiting v, so that they
     // enter and exit w instead, leaving v without any edges entering and exiting 
     // it.  If enter_only, only transfer entering edges.

     void TransferEdges( int v, int w, const Bool enter_only = False );

     // ContractEdgeTo, ContractEdgeFrom: contract the specified edge.
     // Afterwards that edge no longer exists, and all other edges of
     // the far-end vertex are transfered to v.

     void ContractEdgeTo( int w, int j ) 
     {    int v = to_[w][j];
          DeleteEdgeTo( w, j );
	  if ( v != w ) TransferEdges( v, w );    }
     void ContractEdgeFrom( int v, int j ) 
     {    int w = from_[v][j];
          DeleteEdgeFrom( v, j );
	  if ( v != w ) TransferEdges( w, v );    }

     // Glue.  Given two sequences of edges in the graph, given by subpaths a and b,
     //
     //    e0       em-1            f0      fn-1
     // a0 --> ...  --> am,      b0 --> ... --> bn
     //
     // such that {e0,...,em-1} does not meet {f0,...,fn-1},
     //
     // plus a third sequence of edges not in the graph
     //
     //    f0       fp-1  
     // c0 --> ...  --> cp,  
     // 
     // and two injective maps:
     //
     // EE: {0,...m} --> {0,...p},     FF: {0,...,n} --> {0,...,p},
     //
     // such that EE(0) = 0, EE(m) = p, FF(0) = 0, and FF(n) = p,
     //
     // modify the graph by replacing the "a" and "b" sequences by the "c"
     // sequence, using the maps EE and FF to transfer edges involving any of
     // a0,...,am,b0,...,bn.  However, these vertices are left in the graph as
     // vertices having no edges entering or exiting them.
     //
     // Notes on inputs:
     // The digraph c gives c1,...,cp and the intervening edges.

     void Glue( const EmbeddedSubPath<E>& a, const EmbeddedSubPath<E>& b,
          const vec<int>& EE, const vec<int>& FF, const digraphE<E>& c );

     friend void BinaryWrite<E>( int fd, const digraphE<E>& g );
     friend void BinaryRead<E>( int fd, digraphE<E>& g );

     // DOT_vel. This is for a digraphE<int>.  Print in DOT with given vertex 
     // labels in circles and edges labeled with the integers they carry.

     void DOT_vel( ostream& out, const vec<String> & vertex_labels );

     // PrettyDOT: Like digraph::DOT, but prettier...
     // -- Components are labeled (as "contigs", since this is a holdover from
     //    HyperKmerPath and HyperBasevector)
     // -- Longer edges (as indicated by the lengths vector) are drawn longer
     // -- Edges are colored according to their length:
     //    < 100: gray
     //    100-1000: black
     //    1000-10000: red
     //    > 10000: magenta
   
     void PrettyDOT( ostream& out, const vec<double>& lengths,
          Bool label_contigs = True, Bool label_vertices = False,
	  Bool label_edges = False, const vec<int>* componentsToPrint = NULL,
          const Bool edge_labels_base_alpha = True,
	  const vec<String> *label_edges_extra = 0,
	  const vec<String> *label_contigs_extra = 0,
	  const vec<int>* verticesToPrint = NULL ) const;
  
     // Method: DumpGraphML
     // Output the digraph structure in a textual format that can be easily
     // read without reference to our code base.

     void DumpGraphML( const String& graphMLFileName ) const;

     size_t writeBinary( BinaryWriter& writer ) const;
     void readBinary( BinaryReader& reader );

     static size_t externalSizeof() { return 0; }

     private:

     vec<E> edges_;
     vec< vec<int> > to_edge_obj_, from_edge_obj_;

};

template <class E>
struct Serializability< digraphE<E> > : public SelfSerializable {};

// Routines to compare two digraphs.

template<class E> bool operator!=( const digraphE<E>& g1, const digraphE<E>& g2 );
template<class E> bool operator==( const digraphE<E>& g1, const digraphE<E>& g2 );
template<class E> 
     void Compare( ostream& out, const digraphE<E>& g1, const digraphE<E>& g2 );

// DistancesToEnd.  For each vertex v, let D(v) be the maximum (possibly infinity)
// of the lengths over all paths that start with v, where length is defined by
// the function len.  Stop counting as soon as length has reached max_dist.  If
// fw = True, do as just stated.  If fw = False, go in reverse direction.

template<class E> void DistancesToEnd( const digraphE<E>& G,
     int (E::*len)( ) const, const int max_dist, const Bool fw, vec<int>& D );

// LongestPath.  For a given nonempty acyclic graph, find a longest path from
// a source to a sink.  Do not call this on a graph with cycles, as it may
// either assert or run forever.  The answer is a sequence of edge indices.

template<class E> void LongestPath( const digraphE<E>& G, int (E::*len)( ) const, 
     vec<int>& a_longest_path );

// RemoveHangingEnds.  For each edge e, let M(e) be the maximum (possibly
// infinity) of the lengths over all paths that start with e, where length
// is defined by the function len.  If edges ei, ej, ... emanate from a vertex 
// and M(ei) <= max_del and M(ej)/M(ei) >= min_ratio, delete ei.
//
// Unfortunately this can do "bad" things.  For example in the following graph
//
//                    y        x
//        ------> * <----- * -----> *
//                  -----> 
//
// edge x can be deleted using edge y.

template<class E> void RemoveHangingEnds( digraphE<E>& G,
     int (E::*len)( ) const, const int max_del, const double min_ratio );

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

template<class T> void RemoveHangingEnds2( T& h );

// =================================================================================

// Initial minimal implementation of digraphVE.

template<class V, class E> class digraphVE : public digraphE<E> {

     public:

     digraphVE( ) { }

     digraphVE( const vec< vec<int> >& from, const vec< vec<int> >& to,
          const vec<V>& verts, const vec<E>& edges, 
          const vec< vec<int> >& to_edge_obj,
          const vec< vec<int> >& from_edge_obj );

     void Initialize( const vec< vec<int> >& from, const vec< vec<int> >& to,
          const vec<V>& verts, const vec<E>& edges, 
          const vec< vec<int> >& to_edge_obj,
          const vec< vec<int> >& from_edge_obj );

     digraphVE( const digraphE<E>& G, const vec<V>& verts );

     int N( ) const { return verts_.size( ); }

     const V& Vert( int v ) const;
     V& VertMutable( int v );

     size_t writeBinary( BinaryWriter& writer ) const;
     void readBinary( BinaryReader& reader );

     private:

     vec<V> verts_;

};

template <class V, class E>
struct Serializability< digraphVE<V,E> > : public SelfSerializable {};

// =================================================================================

// Class Template: EmbeddedSubPath
//
// An EmbeddedSubPath is a directed path within a given <digraphE>.  It keeps the
// address of the digraph, so should not be used in situations where that address
// might change.
//
// Constructor: takes as input a sequence of consecutive edges
//
//    e0       em-1
// a0 --> ...  --> am,
//
// where the edge from ai to ai+1 is given by from_[ a[i] ][ e[i] ].  

template<class E> class EmbeddedSubPath {

     public:

     void TestValid( ) const;

     // Repair e_ entries which are wrong because digraph has been edited.

     void Repair( );

     EmbeddedSubPath( ) { D_ = 0; }

     EmbeddedSubPath( const digraphE<E>& D, const vec<int>& a, const vec<int>& e )
          : D_(&D), a_(a), e_(e)
     {    esafe_.resize( e.size( ) );
          for ( int i = 0; i < e.isize( ); i++ )
               esafe_[i] = D.EdgeObjectIndexByIndexFrom( a[i], e[i] );
          TestValid( );    }

     // The following constructor takes a vector of vertices.  There must be a
     // unique edge from each vertex to the next.

     EmbeddedSubPath( const digraphE<E>& D, const vec<int>& a ) : D_(&D), a_(a)
     {    ForceAssertGt( a.size( ), 0 );
          e_.resize( a.isize( ) - 1 ), esafe_.resize( a.isize( ) - 1 );
          for ( int i = 0; i < e_.isize( ); i++ )
          {    int v = a[i], w = a[i+1];
               Bool found = False;
               for ( int j = 0; j < D.From(v).isize( ); j++ )
               {    if ( D.From(v)[j] == w )
                    {    ForceAssert( !found );
                         found = True;
                         e_[i] = j;
                         esafe_[i] = D.EdgeObjectIndexByIndexFrom( v, j );    }    }
               if ( !found )
               {    cout << "No edge found from vertex " << v << " to vertex "
                         << w << ".\n";
                    ForceAssert( 0 == 1 );    }    }    }

     int NVertices( ) const { return a_.size( ); }
     int NEdges( ) const { return e_.size( ); }

     int Vertex( int i ) const { return a_[i]; }
     int FirstVertex( ) const { return a_.front( ); }
     int LastVertex( ) const { return a_.back( ); }

     bool IsFromSource( ) const { return D_->Source(FirstVertex()); }
     bool IsToSink( ) const { return D_->Sink(LastVertex()); }

     const E& EdgeObject( int i ) const
     {    AssertGe( i, 0 );
          AssertLt( i, e_.isize( ) );
          return D_->EdgeObjectByIndexFrom( a_[i], e_[i] );    }

     int EdgeObjectIndexAbs( int i ) const
     {    return esafe_[i];    }

     int EdgeObjectIndex( int i ) const
     {    AssertGe( i, 0 );
          AssertLt( i, e_.isize( ) );
          return D_->EdgeObjectIndexByIndexFrom( a_[i], e_[i] );    }

     int EdgeObjectFromIndex( int i ) const
     {    AssertGe( i, 0 );
          AssertLt( i, e_.isize( ) );
          return e_[i];    }

     void SetVertex( int i, int v ) { a_[i] = v; }
     void SetEdge( int i, int e ) 
     {    e_[i] = e; 
          esafe_[i] = D_->EdgeObjectIndexByIndexFrom( a_[i], e );    }

     // Add vertex and edge to left or right.

     void Prepend( int a, int e )
     {    a_.push_front(a);
          e_.push_front(e);
          esafe_.push_front( D_->EdgeObjectIndexByIndexFrom(a,e) );    }
     void Append( int a, int e )
     {    esafe_.push_back( D_->EdgeObjectIndexByIndexFrom(a_.back(),e) );
          a_.push_back(a);
          e_.push_back(e);    }
  
     // Return true if this edge already exists in this EmbeddedSubPath.
     // If true, calling Prepend( a, e ) or Append( a, e ) would produce a loop
     // and/or a duplicated edge.
 
     bool Contains( int edge_id ) 
     {    for ( int i = 0; i < e_.isize( ); i++ )
               if ( edge_id == esafe_[i] ) return true;
          return false;     }

     // Return true if there are any edges in common between p1 and p2.

     friend Bool HasSharedEdge( const EmbeddedSubPath& p1, 
          const EmbeddedSubPath& p2 )
     {    vec<int> edges1, edges2;
          edges1.reserve( p1.NEdges( ) );
          edges2.reserve( p2.NEdges( ) );
          for ( int i = 0; i < p1.NEdges( ); i++ )
          {    edges1.push_back( p1.D_-> 
                    EdgeObjectIndexByIndexFrom( p1.a_[i], p1.e_[i] ) );    }
          for ( int i = 0; i < p2.NEdges( ); i++ )
          {    edges2.push_back( p2.D_-> 
                    EdgeObjectIndexByIndexFrom( p2.a_[i], p2.e_[i] ) );    }
          Sort(edges1), Sort(edges2);
          return Meet( edges1, edges2 );    }

     friend bool operator== ( const EmbeddedSubPath& lhs, 
          const EmbeddedSubPath& rhs ) 
     {    return ( lhs.D_ == rhs.D_ && lhs.a_ == rhs.a_ && lhs.e_ == rhs.e_ &&
                lhs.esafe_ == rhs.esafe_ );    }

     private:

     const digraphE<E>* D_;
     vec<int> a_;
     vec<int> e_;
     vec<int> esafe_;

};

#endif

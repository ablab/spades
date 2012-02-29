///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Define class HyperKmerPath, which is a directed graph whose edges are
// KmerPaths, and other ancillary classes.

#ifndef HYPER_FASTAVECTOR_H
#define HYPER_FASTAVECTOR_H

#include "Fastavector.h"
#include "Superb.h"
#include "graph/Digraph.h"
#include "paths/HyperBasevector.h"

// HyperFastavector.  Like a HyperBasevector, but allows ambiguous base codes.

class HyperFastavector : public digraphE<fastavector> {

     public:

     HyperFastavector( ) { }
     HyperFastavector( int K ) : K_(K) { }

     // Constructor: from K, a digraphE g, and some edge objects.  This ignores
     // the edge objects of g and puts in the new edge objects.

     template<class T> HyperFastavector( int K, const digraphE<T>& g, 
          const vec<fastavector>& edges ) : digraphE<fastavector>( g, edges )
     {    K_ = K;    }

     // Constructor from a HyperBasevector:

     HyperFastavector( const HyperBasevector& h );

     // Constructor: given a collection of fastavectors, create a graph having
     // one edge and two vertices for each of the edge objects

     HyperFastavector( int K, const vec<fastavector>& p ) 
          : digraphE<fastavector>( p, EDGES_SEPARATE )
     {    K_ = K;    }

     // Constructor: form the disjoint union of some HyperFastavectors.
     
     HyperFastavector( int K, const vec<HyperFastavector>& v )
     {    SetK(K);   
          SetToDisjointUnionOf(v);    }
  
     // Constructor: from an HyperFastavector and a list of vertices.

     HyperFastavector( const HyperFastavector& h, const vec<int>& v )
          : digraphE<fastavector>( (const digraphE<fastavector>& ) h, v )
     {    K_ = h.K( );    }

     // Constructor: load from file
     HyperFastavector( const String& filename );

     void Initialize( int K, const vec< vec<int> >& from, 
          const vec< vec<int> >& to, const vec<fastavector>& edges, 
          const vec< vec<int> >& from_edge_obj, const vec< vec<int> >& to_edge_obj )
     {    K_ = K;
          from_ = from;
          to_ = to;
          FromEdgeObjMutable( ) = from_edge_obj;
          ToEdgeObjMutable( ) = to_edge_obj;
          EdgesMutable( ) = edges;    }

     // SetToDisjointUnionOf: clear a given HyperFasta and set it to the 
     // disjoint union of a given collection of HyperFastavector.

     void SetToDisjointUnionOf( const vec<HyperFastavector>& v );

     int K( ) const { return K_; }
     void SetK( int K ) { K_ = K; }

     void LowerK( int newK );

     void RemoveUnneededVertices( );

     // Reverse entire graph.

     void Reverse( );

     int EdgeLength( int e ) const { return EdgeObject(e).size( ); }
     int EdgeLengthKmers( int e ) const { return EdgeObject(e).size( ) - K( ) + 1; }

     // Attempt to linearize this HyperFastavector into a set of superb objects.
     // This process involves many heuristics which are not thoroughly optimized.

     void ConvertToSuperbs( vec<superb> & superbs, const int max_threads,
          const int MIN_EDGE_TO_SAVE ) const;
  
  // Write a fasta-format file which contains the bases in this HyperFastavector
  // scaffolded together (with gaps, etc.) as defined by the scaffolds.
  // If supplied, rc defines which contigs are reverse-complement.
  void WriteScaffoldsFasta( const String & fasta_out,
			    const vec<superb> & scaffolds,
			    const vec<Bool> & rc = vec<Bool>() ) const;
  
     void PrintSummaryPlus( ostream& out, const Bool show_bases = False ) const;

     void PrintSummaryDOT0w( ostream& out,
			     Bool label_contigs = True,
			     Bool label_vertices = False,
			     Bool label_edges = False,
			     const vec<int>* componentsToPrint = NULL,
			     const Bool edge_labels_base_alpha = False,
			     const vec<String> *label_edges_extra = NULL,
			     const vec<String> *label_contigs_extra = NULL,
			     const vec<int> *verticesToPrint = NULL ) const;
     friend void BinaryWrite( int fd, const HyperFastavector& h );
     friend void BinaryRead( int fd, HyperFastavector& h );

     friend void BinaryWrite( int fd, const vec<HyperFastavector>& h )
     {    BinaryWriteComplex( fd, h );    }
     friend void BinaryRead( int fd, vec<HyperFastavector>& h )
     {    BinaryReadComplex( fd, h );    }

     friend Bool operator==( const HyperFastavector& h1, 
          const HyperFastavector& h2 );

     friend ostream& operator<<( ostream& out, const HyperFastavector& h )
     {    for ( int v = 0; v < h.N( ); v++ )
          {    for ( size_t j = 0; j < h.From(v).size(); j++ )
               {    int e = h.EdgeObjectIndexByIndexFrom( v, j );
                    int w = h.From(v)[j];
                    h.EdgeObject(e).Print( out, "edge_" + ToString(e) + "[vert_" 
                         + ToString(v) + "-->vert_" + ToString(w) 
                         + "]" );    }    }
          return out;    }

     private:

     int K_;

};

// Semantic type: edgeName_t
// Alphanumeric name of an assembly edge, obtained as
// BaseAlpha( edge_id ).
SemanticType( String, edgeName_t );

#endif

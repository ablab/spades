///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Define class HyperBasevector, which is a directed graph whose edges are
// basevectors, and other ancillary classes.

#ifndef HYPER_BASEVECTOR_H
#define HYPER_BASEVECTOR_H

#include "Basevector.h"
#include "feudal/BinaryStream.h"
#include "graph/Digraph.h"
#include "paths/HyperKmerPath.h"
#include "paths/KmerBaseBroker.h"
#include "paths/KmerPath.h"

// Class: HyperBasevector
// 
// A HyperBasevector is a <kmer numbering>-independent representation of 
// a HyperKmerPath.  However, it is not independent of K.

class HyperBasevector : public digraphE<basevector> {

     public:

     HyperBasevector( ) { }
     HyperBasevector( int K ) : K_(K) { }

     // Constructor from a HyperKmerPath having no gaps:

     HyperBasevector( const HyperKmerPath& h, const KmerBaseBroker& kbb );

     void Initialize( const HyperKmerPath& h, const KmerBaseBroker& kbb );

     // Constructor: given a collection of basevectors, create a graph having
     // one edge and two vertices for each of the edge objects

     HyperBasevector( int K, const vec<basevector>& p ) 
          : digraphE<basevector>( p, EDGES_SEPARATE )
     {    K_ = K;    }

     // Constructor: form the disjoint union of some HyperBasevectors.
     
     HyperBasevector( int K, const vec<HyperBasevector>& v )
     {    SetK(K);   
          SetToDisjointUnionOf(v);    }

     void Initialize( int K, const vec< vec<int> >& from, 
          const vec< vec<int> >& to, const vec<basevector>& edges, 
          const vec< vec<int> >& from_edge_obj, const vec< vec<int> >& to_edge_obj )
     {    K_ = K;
          from_ = from;
          to_ = to;
          FromEdgeObjMutable( ) = from_edge_obj;
          ToEdgeObjMutable( ) = to_edge_obj;
          EdgesMutable( ) = edges;    }

     // Constructor: given a HyperBasevector, and given a list of vertex indices,
     // create the HyperBasevector having those vertices (with indices starting at 0,
     // but in the given order), and having all the edges that were between those
     // vertices.  Thus this is a "complete subgraph" constructor.

     HyperBasevector( const HyperBasevector& h, const vec<int>& v )
          : digraphE<basevector>( (const digraphE<basevector>& ) h, v )
     {    K_ = h.K( );    }

     explicit HyperBasevector( const String& filename );

     // Check that the edge adjacencies in this HyperBasevector make sense.
     // If they don't, this causes a FatalErr.

     void TestValid( ) const;
  
     // SetToDisjointUnionOf: clear a given HyperBasevector and set it to the 
     // disjoint union of a given collection of HyperBasevectors.

     void SetToDisjointUnionOf( const vec<HyperBasevector>& v );

     // EdgePathToBases.  Given a sequence of edge ids, return the associated
     // basevector.

     basevector EdgePathToBases( const vec<int>& e );

     int K( ) const { return K_; }
     void SetK( int K ) { K_ = K; }

     void LowerK( int newK );

     void RemoveUnneededVertices( );

     int EdgeLength( int e ) const { return EdgeObject(e).size( ); }
     int EdgeLengthBases( int e ) const { return EdgeObject(e).size( ); }
     int EdgeLengthKmers( int e ) const { return EdgeObject(e).size( ) - K( ) + 1; }

     // TotalEdgeLength: return the total number of bases.

     int TotalEdgeLength() const {
       longlong length = 0;
       for ( int e = 0; e < EdgeObjectCount(); e++ )
	 length += EdgeObject(e).size();
       return length;
     }

     // Reverse entire graph.

     void Reverse( );

     void GenerateVecbasevector( vecbvec &bases ) const;

     void RemoveSmallComponents( int min_kmers );

     // Write this HyperBasevector to file.

     void DumpFastb( const String& fn ) const;
     void DumpFasta( const String& fn ) const;

     void PrintSummaryDOT0w( ostream& out, Bool label_contigs = True,
			     Bool label_vertices = False, Bool label_edges = False,
			     const vec<int>* componentsToPrint = NULL,
                             const Bool edge_labels_base_alpha = False,
			     const vec<String> *label_edges_extra = NULL,
			     const vec<String> *label_contigs_extra = NULL,
			     const vec<int> *verticesToPrint = NULL ) const;

     friend void BinaryWrite( int fd, const HyperBasevector& h );
     friend void BinaryRead( int fd, HyperBasevector& h );

     friend void BinaryWrite( int fd, const vec<HyperBasevector>& h )
     {    BinaryWriteComplex( fd, h );    }
     friend void BinaryRead( int fd, vec<HyperBasevector>& h )
     {    BinaryReadComplex( fd, h );    }

     size_t writeBinary( BinaryWriter& writer ) const;
     void readBinary( BinaryReader& reader );
     static size_t externalSizeof();

     friend Bool operator==( const HyperBasevector& h1, const HyperBasevector& h2 );

     friend ostream& operator<<( ostream& out, const HyperBasevector& h )
     {    for ( int v = 0; v < h.N( ); v++ )
          {    for ( size_t j = 0; j < h.From(v).size(); j++ )
               {    int e = h.EdgeObjectIndexByIndexFrom( v, j );
                    int w = h.From(v)[j];
                    h.EdgeObject(e).Print( out, BaseAlpha(e) + " [vert_" 
                         + ToString(v) + "-->vert_" + ToString(w) 
                         + "]" );    }    }
          return out;    }

     private:

     int K_;

};

SELF_SERIALIZABLE(HyperBasevector);

#endif

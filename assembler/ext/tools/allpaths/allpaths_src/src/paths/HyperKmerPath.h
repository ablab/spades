/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// Define class HyperKmerPath, which is a directed graph whose edges are
// KmerPaths, and other ancillary classes.

#ifndef HYPER_KMER_PATH_H
#define HYPER_KMER_PATH_H

#include "Equiv.h"
#include "graph/Digraph.h"
#include "paths/KmerBaseBroker.h"
#include "paths/KmerPath.h"


/**
   Class: HyperKmerPath

   A graph in which each edge is a <KmerPath> ( and each node is just something that
   holds adjacent edges together ).  The graph is a factored representation of a set
   of <KmerPaths>.  The graph *may be disconnected* (i.e. have several connected 
   components).  Note that in this graph, there may be many different edges between 
   a given pair of nodes ( or many different self-loops on a given node ).

   The nodes of this graph have no interpretation of their own; they're just places
   where the edges join.  Each edge, on the other hand, is a KmerPath.

   Some main uses of HyperKmerPaths:

   For each <long-insert pair>, a HyperKmerPath is built representing the possible
   <closures> of the pair.
   
   From each <neighborhood>, <LocalizeReads> constructs a HyperKmerPath representing
   the possible actual sequences of that neighborhood in the genome.  These are then
   joined into increasingly larger graphs until a graph representation of the entire
   assembly is built.

   See also <HyperBasevector>, which is a sequence-space equivalent of HyperKmerPath.
*/

class HyperKmerPath : public digraphE<KmerPath> {

     public:

     // ===========================================================================
     // ======================= INTEGRITY TESTS (partial) =========================
     // ===========================================================================

     void TestValid( ) const;

     // ===========================================================================
     // ====================== CONSTRUCTORS AND THEIR KIN =========================
     // ===========================================================================

     // Constructor 1: empty HyperKmerPath

     HyperKmerPath( ) 
     { K_ = 0; }

     // Constructors 2ab: given a collection of KmerPaths, create a graph having
     // one edge and two vertices for each of the edge objects

     HyperKmerPath( int K, const vec<KmerPath>& p ) 
          : digraphE<KmerPath>( p, EDGES_SEPARATE )
     {    K_ = K;    }

     HyperKmerPath( int K, const vecKmerPath& p )
          : digraphE<KmerPath>( VecOfKmerPath( p ), EDGES_SEPARATE )
     {    K_ = K;    }

     // Constructor 3: given an equivalence relation on a bunch of KmerPaths, build
     // the HyperKmerPath having two vertices per equivalence class, with one edge 
     // between those two vertices for each member of the equivalence class.

     HyperKmerPath( int K, const vecKmerPath& p, const equiv_rel& e )
          : digraphE<KmerPath>( VecOfKmerPath( p ), e )
     {    K_ = K;    }

     // Constructor 4: given a HyperKmerPath, and given a list of vertex indices,
     // create the HyperKmerPath having those vertices (with indices starting at 0,
     // but in the given order), and having all the edges that were between those
     // vertices.  Thus this is a "complete subgraph" constructor.

     HyperKmerPath( const HyperKmerPath& h, const vec<int>& v )
          : digraphE<KmerPath>( (const digraphE<KmerPath>& ) h, v )
     {    K_ = h.K( );    }



     // Constructor 6: from a file.

     explicit HyperKmerPath( const String& filename );

     // Constructor 7: extract a given component from another HyperKmerPath.

     HyperKmerPath( const HyperKmerPath& h, int n )
          : digraphE<KmerPath>( (const digraphE<KmerPath>& ) h, n )
     {    K_ = h.K( );    }

     // Constructor 8: from a HyperKmerPath h and an equivalence relation on its
     // vertices.  This identifies vertices according to the given equivalence
     // relation.

     HyperKmerPath( const HyperKmerPath& h, const equiv_rel& e )
          : digraphE<KmerPath>( (const digraphE<KmerPath>& ) h, e )
     {    K_ = h.K( );    }

     // Constructor 9: from the disjoint union of some HyperKmerPaths.

     HyperKmerPath( int K, const vec<HyperKmerPath>& v )
     {    SetK(K);
          SetToDisjointUnionOf(v);    }

     // SetToDisjointUnionOf: clear a given HyperKmerPath and set it to the disjoint
     // union of a given collection of HyperKmerPaths.

     void SetToDisjointUnionOf( const vec<HyperKmerPath>& v );

     // Constructor 10: from K, a digraphE g, and some edge objects.  This ignores 
     // the edge objects of g and puts in the new edge objects.

     template<class T>
     HyperKmerPath( int K, const digraphE<T>& g, const vec<KmerPath>& edges )
          : digraphE<KmerPath>( g, edges )
     {    K_ = K;    }

     // Constructor 11: from raw data.

     HyperKmerPath( int K, const vec< vec<int> >& from, 
          const vec< vec<int> >& to, const vec<KmerPath>& edges, 
          const vec< vec<int> >& from_edge_obj, const vec< vec<int> >& to_edge_obj )
          : digraphE<KmerPath>( from, to, edges, to_edge_obj, from_edge_obj )
     {    K_ = K;    }

     void Initialize( int K, const vec< vec<int> >& from, 
          const vec< vec<int> >& to, const vec<KmerPath>& edges, 
          const vec< vec<int> >& from_edge_obj, const vec< vec<int> >& to_edge_obj )
     {    K_ = K;
          from_ = from;
          to_ = to;
          FromEdgeObjMutable( ) = from_edge_obj;
          ToEdgeObjMutable( ) = to_edge_obj;
          EdgesMutable( ) = edges;    }

     // Constructor 12: from a collection of HyperKmerPaths and a set of 
     // identifications between vertices in their disjoint union, each of which is 
     // specified as ( (g1,v1), (g2,v2) ) where g1, g2 refer to HyperKmerPaths and 
     // v1, v2 refer to vertices on those HyperKmerPaths.

     HyperKmerPath( int K, const vec<HyperKmerPath>& g,
          const vec< pair< pair<int,int>, pair<int,int> > >& joins )
     {    vec< digraphE<KmerPath> > gg( g.size( ) );
          for ( size_t i = 0; i < g.size(); i++ )
               gg[i] = digraphE<KmerPath>( g[i] );
          digraphE<KmerPath>::Initialize( gg, joins );
          SetK(K);    }

     // Constructor 13: from a given HyperKmerPath and a collection of subsets of
     // its edges, each of which is given the induced subgraph structure, when are
     // then merged into a disjoint union.

     HyperKmerPath( const HyperKmerPath& h, const vec< vec<int> >& C )
          : digraphE<KmerPath>( (const digraphE<KmerPath>& ) h, C )
     {    SetK( h.K( ) );    }

     // Constructor 14: from a given digraphE<KmerPath>.

     HyperKmerPath( const int K, const digraphE<KmerPath>& g )
          : digraphE<KmerPath>(g)
     {    SetK(K);    }

     // Constructor 15: from given digraph, vec<KmerPath>.
  
     HyperKmerPath( const int K, const digraph& g, const vec<KmerPath> &edges )
          : digraphE<KmerPath>(g, edges)
     {    SetK(K);    }

     void Initialize( const int K, const digraph& g, const vec<KmerPath> &edges )
     {    digraphE<KmerPath>::Initialize(g, edges);
          SetK(K);    }

     int K( ) const { return K_; }
     void SetK( int K ) { K_ = K; }

     // ===========================================================================
     // ============================= EDITORS =====================================
     // ===========================================================================

     // Reverse the component containing a given vertex.

     void ReverseComponent( int v );

     // Remove components having less than a specified number of kmers.  For this,
     // the number of kmers in a component is defined to be the sum of the lengths
     // of its edges (which may not be the most useful definition).

     void RemoveSmallComponents( int min_kmers );

     // Remove the loop subgraph from this HKP, leaivng only non-loop edges.
     void RemoveLoopSubgraph();

     // Remove the loop subgraph, and also remove branches as much as possible.
     // The resulting HyperKmerPath should be clean and acyclic.
     void MakeAcyclic();
  
     // Reverse entire graph.

     void Reverse( );

     // If two components are reverse complements of each other,
     // delete one of them.

     void DeleteReverseComplementComponents( );

     // ContractEmptyEdges: For each edge  v ------> w labelled with the
     // empty KmerPath, delete the edge and pull all edges of w into v.

     void ContractEmptyEdges( );

     // ReduceLoops: Wherever we have u ------> v <------> w ------> x, change 
     // v <------> w into a self-loop at v and remove w.

     void ReduceLoops( );


     // MethodDecl: CanonicalizeEdges
     // Canonicalize the KmerPaths on all the edges.
     void CanonicalizeEdges();

     // Zipper: look for two edges that start at the same vertex and go to a 
     // different vertex, and such that the edges agree at the beginning.
     // Merge up to the point where they disagree.  Ditto for reverse.

     void Zipper( );

     // Compress edge objects, so that any adjacent and mergeable
     // KmerPathIntervals are merged.

     void CompressEdgeObjects( );


     // ===========================================================================
     // ============================ PRINTERS =====================================
     // ===========================================================================

     // PrintSummary: generate one line per edge, e.g.
     //     35 --- 523 +/- 8 --> 16
     // would be outputted for an edge from vertex 35 to vertex 16 that has a
     // mean length in kmers of 523 and a variability of 8 (resulting from gaps).

     void PrintSummary( ostream& out ) const;

     // PrintSummaryPlus: like PrintSummary, but also
     // 1. Print kmer ranges for each edge.
     // 2. Organize by graph component and within component, by rough order.

     void PrintSummaryPlus( ostream& out, 
          const void * scratch1 = 0, KmerBaseBroker* kbb = 0, 
	  void * scratch2 = 0, void * scratch3 = 0,
          int scratch4 = 0, Bool print_kmer_ranges = False,
          const vec<String>* edge_remarks = 0, 
          Bool print_component_id_line = True, 
          const vec<String>* component_remarks = 0,
          const vec<Bool>* component_remarks_only = 0 ) const;

     // PrintSummaryPlusAlt: same as PrintSummaryPlus, but more efficient because
     // it doesn't build an equivalence relation.  Takes less arguments.  
     // Numbers components slightly differently.

     void PrintSummaryPlusAlt( ostream& out, const vec<String>* edge_remarks = 0 ) 
          const;

     // Create a fastb file having one record per edge, no gaps allowed.

     void DumpFastb( const String& fn, const KmerBaseBroker& kbb ) const;

     // Create a fasta file having one record per edge, with gaps replaced by Ns.

     void DumpFasta( const String& fn, const KmerBaseBroker& kbb ) const;

     // PrintSummaryDOT: similar to PrintSummary but generate DOT output.
     // PrintSummaryDOT0: similar but don't label edges
     // PrintSummaryDOT0w: similar but don't label edges, do weight them, and 
     // color-code:
     // < 100 kmers: gray
     // 100-1000 kmers: black
     // 1000-10000 kmers: red
     // > 10000 kmers: magenta

     void PrintSummaryDOT( ostream& out ) const;
     void PrintSummaryDOTAlt( ostream& out ) const;
     void PrintSummaryDOT0( ostream& out ) const;
     void PrintSummaryDOT0w( ostream& out,
			     Bool label_contigs = True,
			     Bool label_vertices = False,
			     Bool label_edges = False,
			     const vec<int>* componentsToPrint = NULL,
                             const Bool edge_labels_base_alpha = False,
                             const vec<String> *edge_labels_extra = NULL,
			     const vec<String> *label_contigs_extra = NULL,
			     const vec<int> *verticesToPrint = NULL ) const;



     // ===========================================================================
     // =========================== BINARY I/O ====================================
     // ===========================================================================

     friend void BinaryWrite( int fd, const HyperKmerPath& h );
     friend void BinaryRead( int fd, HyperKmerPath& h );

     friend void BinaryWrite( int fd, const vec<HyperKmerPath>& h )
     {    BinaryWriteComplex( fd, h );    }
     friend void BinaryRead( int fd, vec<HyperKmerPath>& h )
     {    BinaryReadComplex( fd, h );    }

     // ===========================================================================
     // ============================== OTHER ======================================
     // ===========================================================================


     void MakeEdgeDatabase( vec<tagged_rpint>& edgedb ) const;

     int EdgeLength( int e ) const { return EdgeObject(e).KmerCount( ); }
     int EdgeLengthKmers( int e ) const { return EdgeObject(e).KmerCount( ); }

     // TotalEdgeLength: return the total number of kmers in this HyperKmerPath.
     int TotalEdgeLength( ) const {
       int length = 0;
       for ( int e = 0; e < EdgeObjectCount( ); e++ )
	 length += EdgeObject(e).KmerCount( );
       return length;
     }

     // EdgeN50: compute the N50 edge size.

     int EdgeN50( ) const;



     // GetPath(EmbeddedSubPath): return the concatenation of the paths in the
     // subpath.

     KmerPath GetPath( const EmbeddedSubPath<KmerPath>& s )
     {    KmerPath answer;
          for ( int i = 0; i < s.NEdges( ); i++ )
               answer.Append( s.EdgeObject(i) );
          return answer;    }


     // MethodDecl: FindIsomorphicComponents
     // Find HyperKmerPath connected components isomorphic to at least one other
     // connected component.
     // Return the equivalence relation of the components, and the list of ids of those components
     // that have an isomorphic partner in the graph.
     // If you just want to know whether the graph has _any_ isomorphic components,
     // call with stopIfFoundOne==True.
     Bool FindIsomorphicComponents( equiv_rel& componentRelation, vec<int>& isomorphicComponentReps,
				    Bool stopIfFoundOne = False ) const;
     

     // ===========================================================================
     // ============================= PRIVATE =====================================
     // ===========================================================================

     private:

     int K_;

};  // class HyperKmerPath




/// ComponentsAreIsomorphic:
///
/// Given two HyperKmerPaths and one vertex on each, check
/// whether the connected components containing those vertices
/// are isomorphic, with the vertices in correspondence.
/// That is, check whether there is a map from vertices of hkp1
/// to vertices of hkp2 so that corresponding vertices are joined
/// by edges labelled with the same KmerPaths.  (Only the vertices
/// in the connected components of the seed vertices, though.)
///
/// For this to be feasible (in polynomial time :-), we require
/// that the edges in/out of each vertex all bear distinct labels.
/// The function asserts if this is not the case.
///
/// There is an optional vec<int>* argument.  If the components
/// are indeed isomorphic, then this vector will be filled with
/// the isomorphism: v[i]=j means vertex i in hkp1 matches up
/// with vertex j in hkp2, and v[i]=-1 means vertex i was not
/// in the connected component of the seed vertex.
///
/// This is a global function, not a member function.

bool ComponentsAreIsomorphic( const HyperKmerPath& hkp1, int seed_vx_1,
			      const HyperKmerPath& hkp2, int seed_vx_2,
			      vec<int>* p_isomorphism = NULL );




// EdgeLabel generates a label corresponding to a KmerPath.

String EdgeLabel( int v, int w, const KmerPath& p, int K );







// Semantic type: edgeName_t
// Alphanumeric name of an assembly edge, obtained as
// BaseAlpha( edge_id ).
SemanticType( String, edgeName_t );

#endif

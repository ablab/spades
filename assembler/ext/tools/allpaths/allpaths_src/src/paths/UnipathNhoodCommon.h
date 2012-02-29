///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef UNIPATH_NHOOD_COMMON_H
#define UNIPATH_NHOOD_COMMON_H

// File: UnipathNhoodCommon.h
//
// This file defines a toolkit for building a neighborhood (abbreviated "nhood")
// of unipaths around a given seed unipath, and identifying the reads that came from
// this neighborhood.

#include "CoreTools.h"
#include "Basevector.h"
#include "VecUtilities.h"
#include "graph/Digraph.h"
#include "paths/Sepdev.h"
#include "paths/simulation/Placement.h"

// Logical type: pos_rel_to_seed_t
// Position relative to <neighborhood seed>.
SemanticTypeStd( int, pos_rel_to_seed_t );


// Class template: Tedge
//
// An edge defines a separation between two unipaths, along with a deviation value
// for it.
template<typename T>
class Tedge {

     public:

     typedef T value_type;

     int uid1;
     int uid2;
     Tsepdev<T> sepdev;

     T sep() const { return sepdev.Sep(); }
     T dev() const { return sepdev.Dev(); }

     Tedge( ) { }
     template <typename U>
     Tedge( int u1, int u2, U s, U d ) : uid1(u1), uid2(u2), sepdev(s,d) { }

     friend Bool operator<( const Tedge<T>& e1, const Tedge<T>& e2 )
     {    if ( e1.uid1 < e2.uid1 ) return True;
          if ( e1.uid1 > e2.uid1 ) return False;
          if ( e1.uid2 < e2.uid2 ) return True;
          if ( e1.uid2 > e2.uid2 ) return False;
          if ( e1.dev() < e2.dev() ) return True;
          if ( e1.dev() > e2.dev() ) return False;
          return False;    }
  
      friend Bool operator== ( const Tedge<T> &e1, const Tedge<T> &e2 ) {
	if ( e1.uid1 != e2.uid1 ) return false;
	if ( e1.uid2 != e2.uid2 ) return false;
	if ( e1.sep( ) != e2.sep( ) ) return false;
	return ( e1.dev( ) == e2.dev( ) );
      }
};

typedef Tedge<int> edge;


// Given a vec of edges, we can build a sepdev graph:
template<typename T>
void BuildGraphFromEdges( const vec< Tedge<T> >& given_edges,
			  int nuni, // number of unipaths
			  digraphE< Tsepdev<T> >& G ) {
     vec< vec<int> > from(nuni), to(nuni);
     vec< vec<int> > from_edge_obj(nuni), to_edge_obj(nuni);
     vec< Tsepdev<T> > edges;
     for ( int i = 0; i < given_edges.isize( ); i++ )
     {    const Tedge<T>& e = given_edges[i];
          from[e.uid1].push_back(e.uid2);
          to[e.uid2].push_back(e.uid1);
          from_edge_obj[e.uid1].push_back(i);
          to_edge_obj[e.uid2].push_back(i);
          edges.push_back(  Tsepdev<T>( e.sep(), e.dev() ) );    }
     for ( int i = 0; i < nuni; i++ )
     {    SortSync( from[i], from_edge_obj[i] );
          SortSync( to[i], to_edge_obj[i] );    }
     G.Initialize( from, to, edges, to_edge_obj, from_edge_obj );
}










// Class: ustart
//
// A ustart defines a start point for a unipath, relative to a given fixed
// position (that is, relative to the <seed> of the unipath's <neighborhood>).
class ustart {

     public:

     ustart( ) { }
     ustart( int uid, int start, const vec<int>& dev )
          : uid_(uid), start_(start), cached_mean_dev_(-1), dev_(dev)  { }

     int Uid( ) const { return uid_; }

     int Start( ) const { return start_; }

     const vec<int>& Dev( ) const { return dev_; }

     int MeanDev( ) const
     { 
       if ( cached_mean_dev_ < 0 ) {
         double sqsum = 0.0;
         for ( int i = 0; i < dev_.isize( ); i++ )
           sqsum += double( dev_[i] ) * double( dev_[i] );
         cached_mean_dev_ = int(round(sqrt(sqsum)));    
       }
       return cached_mean_dev_;
     }

     friend Bool operator<( const ustart& e1, const ustart& e2 )
     {    return e1.start_ < e2.start_;    }

     struct OrderByDescendingMeanDev 
       : public binary_function<ustart,ustart,bool> {
       bool operator() ( const ustart& lhs, const ustart& rhs ) const {
         return ( lhs.MeanDev() > rhs.MeanDev() );
       }
     };

     private:

     int uid_;
     int start_;
     mutable int cached_mean_dev_;
     vec<int> dev_;

};


class edgeplus : public edge {

public:

     int off1, off2;

     edgeplus( ) { }
     edgeplus( int u1, int u2, int s, int d, int off1, int off2 ) 
          : edge( u1, u2, s, d ), off1(off1), off2(off2) { }

     friend Bool operator<( const edgeplus& e1, const edgeplus& e2 )
     {    if ( e1.uid1 < e2.uid1 ) return True;
          if ( e1.uid1 > e2.uid1 ) return False;
          if ( e1.uid2 < e2.uid2 ) return True;
          if ( e1.uid2 > e2.uid2 ) return False;
          if ( e1.sep( ) < e2.sep( ) ) return True;
          if ( e1.sep( ) > e2.sep( ) ) return False;
          if ( e1.dev( ) < e2.dev( ) ) return True;
          if ( e1.dev( ) > e2.dev( ) ) return False;
          if ( e1.off1 < e2.off1 ) return True;
          if ( e1.off1 > e2.off1 ) return False;
          if ( e1.off2 < e2.off2 ) return True;
          return False;    }

     friend Bool operator==( const edgeplus& e1, const edgeplus& e2 )
     {    return e1.uid1 == e2.uid1 && e1.uid2 == e2.uid2
               && e1.sep( ) == e2.sep( ) && e1.dev( ) == e2.dev( )
               && e1.off1 == e2.off1 && e1.off2 == e2.off2;    }

};


// Print a unipath neighborhood.  Argument "locs" should be valid if USE_TRUTH is true.

void PrintNhood( int v, const vec<ustart>& processed, const vec<int>& ulen,
     Bool USE_TRUTH, const VecPlacementVec* locs );

#endif

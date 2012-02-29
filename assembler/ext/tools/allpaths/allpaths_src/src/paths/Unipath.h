///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef UNIPATH_H
#define UNIPATH_H

// This define Unipath(...).  See the documentation in the main program
// Unipather.cc.  This file defines a callable module.  Output is:
// - unipaths (a vecKmerPath);
// - unipathsdb (a vec<TAGGED_RPINT>, which is some sort of *_tagged_rpint).
// - unipaths_file (optional file of unipaths);
// If unipaths_file is unspecified, then it is not generated.  It is also used as a 
// scratch file, and if specified, performance may be better.

#include "Basevector.h"
#include "CoreTools.h"
#include "paths/HyperBasevector.h"
#include "paths/HyperKmerPath.h"
#include "paths/KmerPath.h"

template<class TAGGED_RPINT>
void Unipath( 
     /* inputs: */  const vecKmerPath& paths, const vecKmerPath& paths_rc,
                    const vec<TAGGED_RPINT>& pathsdb, 
     /* outputs: */ vecKmerPath& unipaths, vec<TAGGED_RPINT>& unipathsdb, 
                    Bool log_progress = False, const String& unipaths_file = "",
                    bool verbose = false );

// The implementation for Unipath is in Unipath.cc, and it must be
// explicitly instaintiated there for each TAGGED_RPINT type.



// DecomposePath: given a KmerPath x, break it up into unipaths, and express x
// in terms of those.

void DecomposePath( const KmerPath& x, Bool brief = False );

// BuildUnipathAdjacencyGraph.  Given a collection of paths p, and a collection of
// unipaths, based on the same kmer numbering, build the graph whose nodes are the
// unipaths and for which there is an edge from u1 to u2 if there exists a path p
// for which the last kmer of u1 is the predecessor of the first kmer of u2.

template<class TAGGED_RPINT>
void BuildUnipathAdjacencyGraph(
     /* inputs: */ const vecKmerPath& paths, 
                   const vecKmerPath& paths_rc,
                   const vec<TAGGED_RPINT>& pathsdb, 
                   const vecKmerPath& unipaths,
                   const vec<TAGGED_RPINT>& unipathsdb, 
     /* output: */ digraph& A );

// BuildUnipathAdjacencyHyperKmerPath.  Take the output of 
// BuildUnipathAdjacencyGraph, and turn it into a HyperKmerPath (changing vertices
// to edges).

void BuildUnipathAdjacencyHyperKmerPath(
     int K, const digraph& A, const vecKmerPath& unipaths, HyperKmerPath& h );

void BuildUnibaseAdjacencyHyperBasevector(
     int K, const digraph& A, const vecbasevector& unibases, HyperBasevector& h );

// UnipathInvolution.  Define map that takes a unipath to its reverse complement.

template<class TAGGED_RPINT>
void UnipathInvolution( const vecKmerPath& unipaths, 
     const vec<TAGGED_RPINT>& unipathsdb, vec<unipath_id_t>& to_rc );

void PrintInColumns( const vec<String>& s );




///////////////////////////////////////////////////////////////////////
//
// Implementations for functions templatized over TAGGED_RPINT classes
//
///////////////////////////////////////////////////////////////////////

/**
   Function: BuildUnipathAdjacencyGraph

   Build a <digraph> whose nodes are unipaths, and there is an edge
   (u1,u2) if there is a (k+1)-mer of the form
   (last-kmer-of(u1), first-kmer-of(u2)) in some read.

   See also: BuildUnipathAdjacencyHyperKmerPath().
*/
template<class TAGGED_RPINT>
void BuildUnipathAdjacencyGraph( const vecKmerPath& paths, 
     const vecKmerPath& paths_rc, const vec<TAGGED_RPINT>& pathsdb, 
     const vecKmerPath& unipaths, const vec<TAGGED_RPINT>& unipathsdb, 
     digraph& A )
{    int nuni = unipaths.size( );
     A.Clear( );
     A.ToMutable( ).resize(nuni), A.FromMutable( ).resize(nuni);
     vec<path_interval_id_t> places;
     vec<kmer_id_t> before, after;
     for ( unipath_id_t i = 0; i < nuni; i++ )
     {    const KmerPath& u = unipaths[i];

          // Find predecessor kmers of the first kmer
          // of this unipath, and successor kmers of
          // the last kmer of this unipath.
     
          kmer_id_t left = u.FirstSegment( ).Start( );
          kmer_id_t right = u.LastSegment( ).Stop( );
          Contains( pathsdb, left, places );
          before.clear( ), after.clear( );
          for ( size_t j = 0; j < places.size(); j++ )
          {    const TAGGED_RPINT& t = pathsdb[ places[j] ];
	    longlong id = t.PathId( );
               const KmerPath& p = ( id >= 0 ? paths[id] : paths_rc[-id-1] );
               int seg = t.PathPos( ), pos = left - t.Start( );
               if ( pos == 0 && seg == 0 ) continue;
               if ( pos > 0 ) --pos;
               else
               {    --seg;
                    pos = p.Length(seg) - 1;    }
               before.push_back( p.Segment(seg).Start( ) + pos );    }
	  UniqueSort( before );
	  if ( before.size( ) > 4 )
	    FatalErr( "BuildUnipathAdjacencyGraph: kmer_id " << left << " seems to have " << before.size( ) << " distinct predecessors in the paths.\nLogically, it cannot have more than 4." );
	  
          Contains( pathsdb, right, places );
          for ( size_t j = 0; j < places.size(); j++ )
          {    const TAGGED_RPINT& t = pathsdb[ places[j] ];
	    longlong id = t.PathId( );
               const KmerPath& p = ( id >= 0 ? paths[id] : paths_rc[-id-1] );
               int seg = t.PathPos( ), pos = right - t.Start( );
               if ( pos == p.Length(seg) - 1 && seg == p.NSegments( ) - 1 ) continue;
               if ( pos < p.Length(seg) - 1 ) ++pos;
               else
               {    ++seg;
                    pos = 0;    }
               after.push_back( p.Segment(seg).Start( ) + pos );    }
          UniqueSort(after);
	  if ( after.size( ) > 4 )
	    FatalErr( "BuildUnipathAdjacencyGraph: kmer_id " << right << " seems to have " << after.size( ) << " distinct successors in the paths.\nLogically, it cannot have more than 4." );

          for ( size_t j = 0; j < before.size(); j++ )
          {    Contains( unipathsdb, before[j], places );
	       ForceAssertLe( places.size(), 1u );  // each kmer is in exactly one unipath
	       if ( places.size( ) != 1 ) continue;
               A.ToMutable( )[i].push_back( 
                    unipathsdb[ places.front() ].PathId( ) );    }
          for ( size_t j = 0; j < after.size(); j++ )
          {    Contains( unipathsdb, after[j], places );
	       ForceAssertLe( places.size(), 1u );  // each kmer is in exactly one unipath
	       if ( places.size( ) != 1 ) continue;
               A.FromMutable( )[i].push_back( 
                    unipathsdb[ places.front() ].PathId( ) );    }    }    
     for ( unipath_id_t i = 0; i < nuni; i++ )
     {    Sort( A.FromMutable(i) );
          Sort( A.ToMutable(i) );    }    }

/**
   Function: UnipathInvolution

   Create a mapping from each unipath to its reverse complement.
   Unlike <read paths> and <genome part paths>, where we normally
   create two arrays -- one for the paths and one for their reverse
   complements, unipaths are normally kept all in one array
   (for historical reasons?) -- hence, the need for this function.
*/
template<class TAGGED_RPINT>
void UnipathInvolution( const vecKmerPath& unipaths, 
     const vec<TAGGED_RPINT>& unipathsdb, vec<unipath_id_t>& to_rc )
{    to_rc.resize( unipaths.size( ) );
     for ( size_t i = 0; i < unipaths.size( ); i++ )
     {    KmerPath u = unipaths[i];
          u.Reverse( );
          vec<unipath_interval_id_t> locs;
	  // each kmer occurs in exactly one unipath;
	  // and so each kmer range (a KmerPathInterval of a unipath)
	  // occurs in exactly one unipath.
          Contains( unipathsdb, u.Segment(0), locs );
          to_rc[i] = unipathsdb[ locs[0] ].PathId( );    }
     for ( size_t i = 0; i < unipaths.size( ); i++ )
          ForceAssertEq( (unipath_id_t) i, to_rc[ to_rc[i] ] );    }

/**
   Function: CheckUnipathSanity

   Check that the unipaths are consistent with their path interval database.
*/
template <class TAG>
void CheckUnipathSanity( const vecKmerPath& unipaths, const vec<TAG>& unipathsdb ) {
  for ( size_t i = 0; i < unipathsdb.size(); i++ ) {
    const TAG& rpi = unipathsdb[i];
    if ( i+1 < unipathsdb.size() )
      ForceAssertLt( rpi.Stop(), unipathsdb[i+1].Start() );
    ForceAssertLe( rpi.Start(), rpi.Stop() );
    for ( kmer_id_t k = rpi.Start(); k <= rpi.Stop(); k++ ) {
      ForceAssertLe( 0, rpi.ReadId() );
      ForceAssertLt( static_cast<size_t>(rpi.ReadId()), unipaths.size() );
      const KmerPath& up = unipaths[ unipath_id_t(rpi.ReadId()) ];
      ForceAssertLt( rpi.PathPos(), up.NSegments() );
      ForceAssertEq( up.Start( rpi.PathPos() ), rpi.Start() );
      ForceAssertEq( up.Stop( rpi.PathPos() ), rpi.Stop() );
    }
  }
}


#endif

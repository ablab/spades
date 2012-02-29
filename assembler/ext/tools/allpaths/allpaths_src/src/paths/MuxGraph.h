// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

/// A MuxGraph tracks the minimal extensions of each path.  We expect
/// this graph to be sparse, but well-connected, i.e. most nodes will
/// have one edge, some will have a few edges, and few will have many.
///
/// For an overview of Mux-based insert walking, see paths/doc/mux.pdf

#ifndef PATHS_MUXGRAPH_H
#define PATHS_MUXGRAPH_H

#include "paths/Mux.h"

#include <map>

class MuxGraph {
 public:
  // Constructors.

  MuxGraph()
  {}

  MuxGraph( int numReads )
  {
    this->resize( numReads );
  }

  int size() const;

  // Resize the graph to accomodate the forward and reverse versions of numReads paths.
  void resize( const int numReads );


  // Methods to get/set edges.

  void GetMuxesOf( const OrientedKmerPathId& pathId, vec<Mux>& muxes ) const;

  int NumMuxesOf( const OrientedKmerPathId& pathId ) const;

  void SetMuxOf( const OrientedKmerPathId& pathId, const Mux& theMux );

  void SetMuxesOf( const OrientedKmerPathId& pathId, const vec<Mux>& muxes );

  // Methods for I/O.

  void Write( const String& filename ) const;

  void Read( const String& filename );
  bool FilesExist( const String& filename ) const; // do all files exist?
  
  bool VerifySameAs( const MuxGraph& other ) const;

  void PrintDot( const String& filename, int partition = 0 ) const;

 private:  

  bool IsSingle( const Mux& aMux ) const { return aMux.GetSegment() >= 0; }
  bool IsSpecial( const Mux& aMux ) const { return aMux.GetSegment() < 0; }

  // Special nodes in the singleEdgeNodes point to nodes in
  // multiEdgeNodes as follows:
  //   start index of nodes = aMux.GetNumKmers()
  //   number of nodes      = -aMux.GetSegment()-1;
  //
  // By this convention, empty nodes are those with aMux.GetSegment() == -1.

  static const int s_empty = -1;

  void SetToEmpty( Mux& aMux ) const { aMux.SetSegment( s_empty ); }
  void SetToMulti( Mux& aMux, int startIndex, int numMuxes ) const
  {
    aMux.SetNumKmers( startIndex );
    aMux.SetSegment( -(numMuxes)-1 );
  }

  int GetStartIndexFromMulti( const Mux& aMux ) const 
  {
    return aMux.GetNumKmers();
  }

  int GetNumMuxesFromSpecial( const Mux& aMux ) const
  {
    return -(aMux.GetSegment())-1;
  }

  void GetFromMulti( const Mux& aMux, int& startIndex, int& numMuxes ) const
  {
    startIndex = GetStartIndexFromMulti( aMux );
    numMuxes   = GetNumMuxesFromSpecial( aMux );
  }

  bool IsEmpty( const Mux& aMux ) const { return aMux.GetSegment() == s_empty; }
  bool IsMulti( const Mux& aMux ) const { return aMux.GetSegment() < s_empty; }

  vec<Mux> m_singleEdgeNodes;
  vec<Mux> m_multiEdgeNodes;
};

#endif

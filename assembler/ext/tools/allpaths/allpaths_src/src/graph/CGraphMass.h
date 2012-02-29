/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// Since this has KmerPath stuff in it, it should be in directory 'paths'
// rather than this directory.

#ifndef C_GRAPH_MASS_H
#define C_GRAPH_MASS_H

#include "String.h"
#include "Vec.h"
#include "paths/KmerPath.h"
#include "paths/Unipath.h"

/**
 * class CGraphMass
 *
 * The kmer-mass of an edge is defined as the recursive total
 * kmer-length of all the edges following (or preceding) the given
 * edge (edge included).  For efficiency reasons, the mass is
 * computed up to a fixed cap.
 */
class CGraphMass {

public:

  CGraphMass( const HyperKmerPath *hkp = 0, int cap = 3000 );
  
  void SetPointers( const HyperKmerPath *hkp );
  
  void SetCap( int cap ) { cap_ = cap; }

  int Cap( ) const { return cap_; }
  
  // Find (unique) vertex at the left of the given edge.
  int LeftVertexIdx( int edge_id ) const { return toleft_idx_[edge_id]; }

  // Find (unique) vertex at the right of the given edge.
  int RightVertexIdx( int edge_id ) const { return toright_idx_[edge_id]; }

  // Forward (or reverse) mass of given edge (mass of edge is included).
  int DirectedKMass( int edge_id, bool fw = true ) const;
  
  
private:

  // Build indexes (toleft_idx_ and toright_idx_).
  void BuildIndexes( );

  // Core recursive method (return true if mass >= cap_).
  bool KMass( int edge_id, int &mass, vec<int> &seen ) const;
  
  // Add mass of edge (skip already seen edges, return true if mass >= cap_).
  bool AddEdge( int edge_id, int &mass, vec<int> &seen ) const;

  // Find all edges from and to vertex_id.
  vec<int> FindAllEdges( int vertex_id ) const;
  
  
private:
  
  const HyperKmerPath *hkp_;    // the HyperKmerPath
  int cap_;                     // kmer-mass cap for the recursive algorithm

  vec<int> toleft_idx_;         // index of vertices at left of edges
  vec<int> toright_idx_;        // index of vertices at right of edges
  
};

#endif

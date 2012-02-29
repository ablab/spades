///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef C_SUPER_LINKS_H
#define C_SUPER_LINKS_H

#include "Intvector.h"
#include "PairsManager.h"
#include "Superb.h"
#include "VecTemplate.h"
#include "math/Functions.h"
#include "paths/Alignlet.h"
#include "paths/reporting/COffset.h"

/**
 * class CSuperLinks
 */
class CSuperLinks {
  
public:
  
  CSuperLinks( ) { }

  CSuperLinks( const PairsManager *pairs,
	       const vec<superb> *supers,
	       const vec<alignlet> *aligns,
	       const vec<int> *index );
  
  void SetPointers( const PairsManager *pairs,
		    const vec<superb> *supers,
		    const vec<alignlet> *aligns,
		    const vec<int> *index );
  
  // Find all links from super1 to (rc2 of) super2.
  COffset AllLinks( int super1, int super2, bool rc2 ) const;

  // Find all links from super_id to other scaffolds.
  void AllLinks( int super_id, vec<COffset>& out,
		 float slop = -1, float *stretch = 0 ) const;
  
  // Print all links from super_id (very large output if full = true).
  void PrintAllLinks( ostream &out, int super_id, bool full = false ) const;

  // Returns ids of all pairs linking this super with a different super.
  vec<int> AllPairs( int super_id ) const;
  
  // Add link from pair id to the given COffset.
  void AddLink( int pair_id, COffset &offset ) const;
  
  
private:

  // Centralized factory for the maps in the class.
  void GenerateMaps( );
  
  // Window on super (these assert if read does not align).
  int SuperId( int read_id ) const;
  bool FwOnSuper( int read_id ) const;
  pair<int,int> WinOnSuper( int read_id ) const;
  
  
private:

  const PairsManager *pairs_;
  const vec<superb> *supers_;
  const vec<alignlet> *aligns_;
  const vec<int> *index_;

  vec<int> super_id_;     // id of super containing this edge
  vec<int> super_pos_;    // position in super for this edge
  vec<int> super_begin_;  // begin of edge on super
  vec<int> super_end_;    // end of edge on super
  VecIntVec read_ids_;    // reads in this edge that own an align

  vec<int> true_begin_;   // cached true begin for this super
  vec<int> true_end_;     // cached true end for this super
  
};

#endif

///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef _INSERT_WALKER
#define _INSERT_WALKER

/*******************************************************************************
 *
 *        InsertWalker
 *
 * The InsertWalker module is a tool for "walking" long inserts along a unipath
 * graph, finding the set of all paths from one end to another.  This module is
 * designed for use with LocalizeReadsLG, as a replacement for the old Mux-based
 * insert walking machinery.
 *
 * This module contains two classes: InsertWalker and InsertWalkerLoc (which is
 * only in InsertWalker.cc.)  The main function is InsertWalker::WalkUnipaths.
 *
 *
 * Josh Burton
 * December 2009
 *
 ******************************************************************************/

#include <set>
#include "graph/Digraph.h"
#include "math/HoInterval.h"
#include "paths/HyperKmerPath.h"
#include "paths/KmerPath.h"
#include "TaskTimer.h" // TaskTimer






/*******************************************************************************
 *
 * The InsertWalker is a simple container class.  It loads in a set of data
 * structures that it needs for insert walking, and then it performs the insert
 * walking with the function WalkUnipaths().  The function WalkUnipaths() uses
 * the local InsertWalkerLoc class.
 *
 ******************************************************************************/
class InsertWalker {
  
public:
  
  // Constructor: Supplies pointers to the needed data structures.
  InsertWalker( const int K, 
		const vecKmerPath * unipaths,
		const vec<tagged_rpint> * unipathsdb,
		const digraph * unipath_adjs );
  
  // PlaceReadsOnUnipaths: Locate a pair of reads on the unipath graph.
  void PlaceReadsOnUnipaths( const KmerPath & read1, const KmerPath & read2,
			     int & u1, int & u2, int & sep_offset ) const;
  
  // WalkUnipaths:  Attempt to bridge the separation between two unipaths
  // by walking along the unipath graph.
  // If the insert cannot be walked for any reason, this returns an empty set.
  // This is the main workhorse function in the InsertWalker class.
  // For a more thorough documentation of this algorithm, see InsertWalker.cc.
  set<int> WalkUnipaths( const int & u1, const int & u2, const ho_interval & sep_range,
			 TaskTimer & timer,
			 const bool verbose = false ) const;
  
  
private:
  
  // Pointers to outside objects
  const int _K;
  const vecKmerPath * _unipaths;
  const vec<tagged_rpint> * _unipathsdb;
  const digraph * _unipath_adjs;
  
  // Derived data
  vec<int> _unipath_component;
};



#endif

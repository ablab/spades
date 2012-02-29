// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

#ifndef PATHS_GAP_SKIP
#define PATHS_GAP_SKIP

#include "Vec.h"
#include "paths/KmerPath.h"
// just need a forward declaration here; the .cc needs the full include
class NegativeGapValidator;


// The gap skipper returns a FarEnd for each possible skip.
//
// If the reads keep going beyond the merger, then:
//  .DONE is False
//  .loc1 and .loc2 point to the common k-mers of the given KmerPaths
//       which align on the far side of the gap
//  .merged is the merged path in between the orig locs and the returned ones;
//       for the rest of the path, recurse on .loc1 and .loc2.
//
// If one read ends before the two reads next share a kmer, then:
//  .DONE is True
//  .merged is the entire path from the orig locs to the end; don't recurse
//  .stop{1,2} indicate where path{1,2} ended in .merged -- measured in
//      segments away from the end, so one of them is 0, "last segment".
//      (This is the index if at the left edge, else NSegments-1-index.)


struct FarEnd {
  // Want a vec of these, so I need a default constructor
  FarEnd( ) : DONE( False ) { }
  FarEnd( const KmerPathLoc& _loc1, const KmerPathLoc& _loc2 ) :
    loc1(_loc1), loc2(_loc2), DONE( False ) { }

  KmerPathLoc loc1;
  KmerPathLoc loc2;
  KmerPath merged;
  Bool DONE; // True: stop recursing; "merged" goes all the way to the end.
  int stop1;
  int stop2;

  void flip() {
    swap(loc1,loc2);
    swap(stop1,stop2);
  }

};

// NOTE: in FarEnds return by the gap skipping routines, we require that
// loc1 and loc2 point to the same KmerPaths as loc1orig and loc2orig,
// respectively.  If the order of arguments is switched in a recursive
// call, then the caller should call the .flip() method on the return
// values to get things in the right order.


void GapSkipLeft( KmerPathLoc loc1orig,
		  KmerPathLoc loc2orig,
		  vec<FarEnd>& far_ends,
		  const NegativeGapValidator* ngv = NULL );
void GapSkipRight( KmerPathLoc loc1orig,
		   KmerPathLoc loc2orig,
		   vec<FarEnd>& far_ends, 
		   const NegativeGapValidator* ngv = NULL );


void GapOnGapRight( KmerPathLoc loc1orig,
		    KmerPathLoc loc2orig,
		    int min_gap_used, int max_gap_used,
		    vec<FarEnd>& far_ends,
		    vec< pair<int,int> >& real_gap_used,
		    const NegativeGapValidator* ngv = NULL );
void GapOnGapLeft( KmerPathLoc loc1orig,
		   KmerPathLoc loc2orig,
		   int min_gap_used, int max_gap_used,
		   vec<FarEnd>& far_ends,
		   vec< pair<int,int> >& real_gap_used,
		   const NegativeGapValidator* ngv = NULL );



#endif


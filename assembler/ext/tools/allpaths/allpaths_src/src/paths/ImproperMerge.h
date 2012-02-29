// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

#ifndef IMPROPER_MERGE
#define IMPROPER_MERGE

#include "CoreTools.h"
#include "Vec.h"

#include "paths/KmerPath.h"
#include "paths/NegativeGapValidator.h"

/// For each alignment, the aligner return an ImproperMerger struct,
/// containing a KmerPath and four KmerPathLoc's:
///
/// .merged is the merged agreeing section of the two reads.  
///   This extends to the last matching kmer in both directions.
///   NOTE that this may not be the end of the original paths, even
///   if the alignment is genuinely proper(!), since the end of one
///   path may land in a gap on the other path.
///
/// .left_end1, .left_end2, .right_end1, .right_end2
///   Each is a KmerPathLoc, pointing to the leftmost/rightmost
///   matched kmer on the first/second aligned read.
///
/// The struct also has two member functions, PushLoc{1,2}(loc),
/// which can take a KmerPathLoc on the matched section of one
/// of the two original paths, and return the corresponding
/// loc on .merged.
///
/// If you plan to do negative gap validation on a result, only the
/// .merged member need be validated.

struct ImproperMerger {

  KmerPath merged;
  
  KmerPathLoc left_end1;
  KmerPathLoc left_end2;
  KmerPathLoc right_end1;
  KmerPathLoc right_end2;

  KmerPathLoc PushLoc1( const KmerPathLoc& loc1 ) const;
  KmerPathLoc PushLoc2( const KmerPathLoc& loc2 ) const;

  void flip() {
    swap( left_end1, left_end2 );
    swap( right_end1, right_end2 );
  }

};



/// The improper path merger tool takes as arguments two paths
/// which contain a shared kmer, and finds the stretches of those
/// two reads which can be aligned and merged.  The set of such
/// alignments is always nonempty, since at the very least the
/// matching stretch contains the original shared kmer.  (But with
/// the min_perfect_match argument you can demand more than that.)
/// As usual, it is possible to get more than one alignment, in the
/// presence of stretchy gaps.
///
/// The aligned stretch in each direction reaches either to the end of 
/// one of the reads ("is proper") or else extends to the last matched 
/// kmer on the two given reads.  There is some subtlety in this def'n,
/// since we could just as well consider the aligned portion to include
/// gaps just past the last matched kmer, and only end at the first
/// unmatchable kmer instead.  This implements the first version because
/// it is much easier, not because it is more correct.

// p1.Segment(ind1) and p2.Segment(ind2) are overlapping sequence.

void ImproperMergePaths( const KmerPath& p1, const KmerPath& p2, 
			 int ind1, int ind2, vec<ImproperMerger>& ans,
			 int min_perfect_match=1,
			 const NegativeGapValidator* ngv = NULL );


#endif

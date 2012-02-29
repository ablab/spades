// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

#include "CoreTools.h"
#include "Vec.h"

#include "paths/KmerPath.h"
#include "paths/AlignAndMerge.h"
#include "paths/GapSkip.h"

#include "paths/VerifyGapLengths.h"

// I doubt anyone but MergePaths will ever call these helpers,
// so I'm making them private to this file.

void RecursiveAlignAndMergeLeft( const KmerPathLoc& loc1orig,
				 const KmerPathLoc& loc2orig,
				 vec<MergedKmerPath>& ans,
				 const NegativeGapValidator* ngv = NULL );
void RecursiveAlignAndMergeRight( const KmerPathLoc& loc1orig,
				  const KmerPathLoc& loc2orig,
				  vec<MergedKmerPath>& ans,
				  const NegativeGapValidator* ngv = NULL );
// return value: length of the longest perfect match to the left/right.
// MergePaths uses this to implement min_perfect_match where the long 
// match region doesn't need to contain the original anchor point



// I believe this aligner finds all possible mergers.
// Warning: the number of mergers can be very large; it may grow
// exponentially with the length of the input.  Don't use blindly.

// I probably made a bad design decision in writing everything
// twice, once for moving left and once for moving right.  This
// is intended to increase running efficiency, but requires writing and
// maintaining two only slightly different versions of everything.

// If it turns out that MergePaths is being used in cases where
// the answer involves millions of possible mergers, then
// the answer should be returned in a vecKmerPath (mastervec),
// changing the code below so that things don't all end up self-owned.
// Not dealing with this for now; I hope it won't be an issue.


/// MergePaths() takes two KmerPaths and tries to merge them, starting from
/// the overlapping KmerPathIntervals p1.Interval(ind1) and p2.Interval(ind2).
/// For each merger it finds, it pushes a MergedKmerPath object onto ans.

// VERSION 2: accepts an optional argument of a NegativeGapValidtor*.
// If non-NULL, this NGV is used when a new gap is created, to make sure
// the gap's minimum length is NGV-approved.  The point is to avoid, eg:
//
// path1: ...---------------- (   g a p  7 5 - 7 5   ) ---------------
// path2: ...------- (very stretchy gap) --
// merger: ..---------------- (gap 0-74) -- (gap 0-74) --------
// ngv: .....---------------- (gap 47-74)--(gap 47-74) --------
//
// In this example, the two kmers added into the middle don't belong at all,
// but the NGV is too late to stop them.  Note that now the original path
// doesn't even align to the post-NGV version of the merger!

void MergePaths( const KmerPath& p1, const KmerPath& p2, 
		 int ind1, int ind2, 
		 vec<MergedKmerPath>& ans,
		 int min_perfect_match, 
		 const NegativeGapValidator* ngv,
		 bool DEBUG_GAP_SIZES ) {

  Assert( min_perfect_match >= 1 );  // Anything 0 or negative is nonsense

  ans.clear();

  // The given intervals must be overlapping sequence.
  if ( !p1.isSeq(ind1) || !p2.isSeq(ind2) 
       || !p1.Segment(ind1).Overlaps(p2.Segment(ind2)) )
    return;  // bad anchor point

  // Create KmerPathLocs pointing to an aligning k-mer in each path.
  // If the alignment isn't proper, return immediately (with ans still empty)
  pair<KmerPathLoc,KmerPathLoc> lp = CreateAlignedLocs( p1, p2, ind1, ind2);

  // Look left:
  KmerPathLoc left1 = lp.first;
  KmerPathLoc left2 = lp.second;
  if ( ! ScanLeftPerfectMatch(left1,left2) )
    return;

  // Look right:
  KmerPathLoc right1 = lp.first;
  KmerPathLoc right2 = lp.second;
  if ( ! ScanRightPerfectMatch(right1,right2) )
    return;

  int perfect_match_len = right1 - left1 + 1;
  // but there may be longer matches on the far sides of gaps!

  // Finish the no-gap case ASAP:
  if(  (left1.atBegin() || left2.atBegin()) &&
       (right1.atEnd() || right2.atEnd())  ) {

    if( perfect_match_len < min_perfect_match ) 
      return;

    ans.resize(1);  // now ans[0] holds an empty MergedKmerPath
    ans[0].longest_perfect_match = perfect_match_len;

    // Copy first half of path, and set leftmost and given:
    if( left1.atBegin() ) {
      p2.CopySubpath( p2.Begin(), right2, ans[0].path );
      ans[0].left_end = make_pair( left2.GetIndex(), 0 );
      ans[0].given = ind2;
    } else {
      p1.CopySubpath( p1.Begin(), right1, ans[0].path );
      ans[0].left_end = make_pair( 0, left1.GetIndex() );
      ans[0].given = ind1;
    }

    // Copy second half of path, and set right_end:
    if( right1.atEnd() ) {
      p2.CopySubpathNoFirstKmer( right2, p2.End(), ans[0].path );
      ans[0].right_end = make_pair( p2.NSegments()-1 - right2.GetIndex(), 0);
    } else {
      p1.CopySubpathNoFirstKmer( right1, p1.End(), ans[0].path );
      ans[0].right_end = make_pair( 0, p1.NSegments()-1 - right1.GetIndex() );
    }
    // We're done
    return;
  }


  // Otherwise, we hit a gap, on either the left or right.
  // Begin the recursive aligning in both directions, and
  // stitch together all the answers:

  // Hold on to the perfectly-matched middle:
  // (and how far the matching segment is from its right-hand side)
  KmerPath middle;
  p1.CopySubpathNoFirstKmer( left1, right1, middle );
  int match_offset = right1.GetIndex() - lp.first.GetIndex();

  // containers for the left and right possibilities:
  vec<MergedKmerPath> left_ans, right_ans;

  // Look left:
  RecursiveAlignAndMergeLeft( left1, left2, left_ans, ngv );
  // Don't take the time to search right if there were no left alignments:
  if( left_ans.empty() ) 
    return;   // no left alignments

  // Look right:
  RecursiveAlignAndMergeRight( right1, right2, right_ans, ngv );

  if( right_ans.empty() )
    return;   // no right alignments;

  // Glue together the lefts, the middle, and the rights, in all ways:
  ans.reserve( left_ans.size() * right_ans.size() );
  int max_match_len;

  for(uint i=0; i<left_ans.size(); i++)
    for(uint j=0; j<right_ans.size(); j++) {

      // Check that there's a long-enough perfect match:
      max_match_len = max( perfect_match_len, 
			   max( left_ans[i].longest_perfect_match,
				right_ans[j].longest_perfect_match ) );
      if( max_match_len < min_perfect_match )
	continue;

      // left:
      ans.push_back( left_ans[i] );  // sets .path and .left_end
      // middle (which was already NoFirstKmer'ed):
      ans.back().path.Append(middle);
      // This is the only moment when we can set .given easily:
      ans.back().given = ans.back().path.NSegments() - match_offset - 1;
      // right:
      ans.back().path.AppendNoFirstKmer( right_ans[j].path );
      ans.back().right_end = right_ans[j].right_end;
      // match len:
      ans.back().longest_perfect_match = max_match_len;
    }

  if( DEBUG_GAP_SIZES )
    for( vec<MergedKmerPath>::iterator merger = ans.begin();
	 merger != ans.end(); merger++ )
      VerifyGapLengths(*merger, p1, p2);

}



void RecursiveAlignAndMergeLeft( const KmerPathLoc& loc1orig,
				 const KmerPathLoc& loc2orig,
				 vec<MergedKmerPath>& ans,
				 const NegativeGapValidator* ngv ) {

  KmerPathLoc loc1=loc1orig, loc2=loc2orig;
  // Scan perfect match; return nothing if improper.
  if ( ! ScanLeftPerfectMatch(loc1,loc2) )
    return;

  int match_len = loc1orig - loc1 + 1;

  // Base case: We got to the end of one of the reads
  if ( loc1.atBegin() ) {
    ans.resize(1);
    loc2.GetPath().CopyHead( loc2orig, ans[0].path );
    ans[0].left_end = 
      make_pair( loc2.GetIndex(), 0);
    ans[0].longest_perfect_match = match_len;
    return;
  }
  if ( loc2.atBegin() ) {
    ans.resize(1);
    loc1.GetPath().CopyHead( loc1orig, ans[0].path );
    ans[0].left_end = 
      make_pair( 0, loc1.GetIndex() );
    ans[0].longest_perfect_match = match_len;
    return;
  }

  // Otherwise, we must have reached a gap.  Call the gap skipper.

  KmerPath initial_perfect_match;
  loc1.GetPath().CopySubpath( loc1, loc1orig, initial_perfect_match );

  vec<FarEnd> far_ends;

  GapSkipLeft( loc1, loc2, far_ends, ngv );

  // Recursively call self on each far_end; build answers in ans.
  vec<MergedKmerPath> sub_answers;

  for( vec<FarEnd>::iterator far = far_ends.begin();
       far != far_ends.end(); far++ ) {
    if( far->DONE ) {   // Don't recurse: far->merged goes to end
      ans.push_back( MergedKmerPath(far->merged) );
      ans.back().path.Append( initial_perfect_match );
      ans.back().left_end = make_pair(far->stop1, far->stop2);
      ans.back().longest_perfect_match = match_len;
    } else {
      sub_answers.clear();
      RecursiveAlignAndMergeLeft( far->loc1, far->loc2, sub_answers, ngv );
      for( vec<MergedKmerPath>::iterator sub_ans = sub_answers.begin();
	   sub_ans != sub_answers.end(); sub_ans++ ) {
	// recursive answer + merged gap section + perfect match
	ans.push_back( sub_ans->path );
	ans.back().path.AppendNoFirstKmer( far->merged );
	ans.back().path.Append( initial_perfect_match );
	ans.back().left_end = sub_ans->left_end;
	ans.back().longest_perfect_match = 
	  max( match_len, sub_ans->longest_perfect_match );
      }
    }
  }
}



void RecursiveAlignAndMergeRight( const KmerPathLoc& loc1orig,
				  const KmerPathLoc& loc2orig,
				  vec<MergedKmerPath>& ans,
				  const NegativeGapValidator* ngv ) {

  KmerPathLoc loc1=loc1orig, loc2=loc2orig;
  // Scan perfect match; return nothing if improper.
  if ( ! ScanRightPerfectMatch(loc1,loc2) )
    return;

  int match_len = loc1 - loc1orig + 1;

  // Base case: We got to the end of one of the reads
  if ( loc1.atEnd() ) {
    ans.resize(1);
    loc2.GetPath().CopyTail( loc2orig, ans[0].path );
    ans[0].right_end = 
      make_pair( loc2.GetPath().NSegments()-1 - loc2.GetIndex(), 0);
    ans[0].longest_perfect_match = match_len;
    return;
  }
  if ( loc2.atEnd() ) {
    ans.resize(1);
    loc1.GetPath().CopyTail( loc1orig, ans[0].path );
    ans[0].right_end = 
      make_pair( 0, loc1.GetPath().NSegments()-1 - loc1.GetIndex() );
    ans[0].longest_perfect_match = match_len;
    return;
  }

  // Otherwise, we must have reached a gap.  Call the gap skipper.

  KmerPath initial_perfect_match;
  loc1.GetPath().CopySubpath( loc1orig, loc1, initial_perfect_match );

  vec<FarEnd> far_ends;

  GapSkipRight( loc1, loc2, far_ends, ngv );

  // Recursively call self on each far_end; build answers in ans.
  vec<MergedKmerPath> sub_answers;

  for( vec<FarEnd>::iterator far = far_ends.begin();
       far != far_ends.end(); far++ ) {
    if( far->DONE ) {   // Don't recurse: far->merged goes to end
      ans.push_back( MergedKmerPath(initial_perfect_match) );
      ans.back().path.Append( far->merged );
      ans.back().right_end = make_pair(far->stop1, far->stop2);
      ans.back().longest_perfect_match = match_len;
    } else {      
      sub_answers.clear();
      RecursiveAlignAndMergeRight( far->loc1, far->loc2, sub_answers, ngv );
      for( vec<MergedKmerPath>::iterator sub_ans = sub_answers.begin();
	   sub_ans != sub_answers.end(); sub_ans++ ) {
	// perfect match + merged gap section + recursive answer
	ans.push_back( initial_perfect_match );
	ans.back().path.Append( far->merged );
	ans.back().path.AppendNoFirstKmer( sub_ans->path );
	ans.back().right_end = sub_ans->right_end;
	ans.back().longest_perfect_match = 
	  max( match_len, sub_ans->longest_perfect_match );
      }
    }
  }
}


// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

#include "Vec.h"
#include "paths/KmerPath.h"
#include "paths/GapSkip.h"

#include "paths/ImproperMerge.h"
#include "paths/PathEmbedding.h"


// I doubt anyone but ImproperMergePaths will ever call these helpers,
// so I'm making them private to this file.

void RecursiveImproperMergeLeft( const KmerPathLoc& loc1orig,
				 const KmerPathLoc& loc2orig,
				 vec<ImproperMerger>& ans,
				 const NegativeGapValidator* ngv );
void RecursiveImproperMergeRight( const KmerPathLoc& loc1orig,
				  const KmerPathLoc& loc2orig,
				  vec<ImproperMerger>& ans,
				  const NegativeGapValidator* ngv );



void ImproperMergePaths( const KmerPath& p1, const KmerPath& p2, 
			 int ind1, int ind2, vec<ImproperMerger>& ans,
			 int min_perfect_match, 
			 const NegativeGapValidator* ngv ) {

  ans.clear();

  // The given intervals must be overlapping sequence.
  if ( !p1.isSeq(ind1) || !p2.isSeq(ind2) 
       || !p1.Segment(ind1).Overlaps(p2.Segment(ind2)) )
    return;

  // Create KmerPathLocs pointing to an aligning k-mer in each path.
  // If the alignment isn't proper, return immediately (with ans still empty)
  pair<KmerPathLoc,KmerPathLoc> lp = CreateAlignedLocs( p1, p2, ind1, ind2);

  // Look left:
  KmerPathLoc left1 = lp.first;
  KmerPathLoc left2 = lp.second;
  ScanLeftPerfectMatch(left1,left2);

  // Look right:
  KmerPathLoc right1 = lp.first;
  KmerPathLoc right2 = lp.second;
  ScanRightPerfectMatch(right1,right2);

  // Is the perfect match long enough?
  int perfect_match_len = right1 - left1 + 1;
  if( perfect_match_len < min_perfect_match )
    return;

  // Next we'll call the recusive aligners to the left and right.
  // If we didn't bump into a gap in one direction, this is
  // a bit inefficient, since we could already know now that
  // we'd reached the end.  In AlignAndMerge we do our best to
  // speed things up, and handle this without the extra function 
  // calls, but here we're lazy.

  // Hold on to the perfectly-matched middle:
  KmerPath middle;
  p1.CopySubpathNoFirstKmer( left1, right1, middle );

  // Call the recursive helpers:
  vec<ImproperMerger> left_ans, right_ans;
  RecursiveImproperMergeLeft( left1, left2, left_ans, ngv );
  RecursiveImproperMergeRight( right1, right2, right_ans, ngv );

  // Glue together the lefts, the middle, and the rights, in all ways:
  ans.reserve( left_ans.size() * right_ans.size() );

  for(uint i=0; i<left_ans.size(); i++)
    for(uint j=0; j<right_ans.size(); j++) {
      // left:
      ans.push_back( left_ans[i] );  // sets .path and .left_end*
      // middle (which was already NoFirstKmer'ed):
      ans.back().merged.Append(middle);
      // right:
      ans.back().merged.AppendNoFirstKmer( right_ans[j].merged );
      ans.back().right_end1 = right_ans[j].right_end1;
      ans.back().right_end2 = right_ans[j].right_end2;
    }

  return;
}




void RecursiveImproperMergeLeft( const KmerPathLoc& loc1orig,
				 const KmerPathLoc& loc2orig,
				 vec<ImproperMerger>& ans,
				 const NegativeGapValidator* ngv ) {

  KmerPathLoc loc1=loc1orig, loc2=loc2orig;

  // If this scan doesn't reach a gap, we're done!
  if ( ScanLeftPerfectMatch(loc1,loc2) == false  // kmer mismatch
       || loc1.atBegin() 
       || loc2.atBegin() ) {
    ans.resize(1);
    ans[0].left_end1 = loc1;
    ans[0].left_end2 = loc2;
    loc1.GetPath().CopySubpath( loc1, loc1orig, ans[0].merged );
    return;
  }

  // Otherwise, we must have reached a gap.  Call the gap skipper.

  KmerPath initial_perfect_match;
  loc1.GetPath().CopySubpath( loc1, loc1orig, initial_perfect_match );

  vec<FarEnd> far_ends;

  GapSkipLeft( loc1, loc2, far_ends, ngv );

  // If there's no way to skip the gap, we're done!
  if( far_ends.empty() ) {
    ans.resize(1);
    ans[0].left_end1 = loc1;
    ans[0].left_end2 = loc2;
    ans[0].merged = initial_perfect_match;
    return;
  }

  // Recursively call self on each far_end; build answers in ans.
  vec<ImproperMerger> sub_answers;

  for( vec<FarEnd>::iterator far = far_ends.begin();
       far != far_ends.end(); far++ ) {
    if( far->DONE ) {  // there was not another shared kmer -- done!
      ans.resize( ans.size()+1 );
      ans.back().merged = initial_perfect_match;
      ans.back().left_end1 = loc1;
      ans.back().left_end2 = loc2;
    } else {
      sub_answers.clear();
      RecursiveImproperMergeLeft( far->loc1, far->loc2, sub_answers, ngv );
      for( vec<ImproperMerger>::iterator sub_ans = sub_answers.begin();
	   sub_ans != sub_answers.end(); sub_ans++ ) {
	// recursive answer + merged gap section + perfect match
	ans.push_back( *sub_ans );
	ans.back().merged.AppendNoFirstKmer( far->merged );
	ans.back().merged.Append( initial_perfect_match );
      }
    }
  }
}


  
void RecursiveImproperMergeRight( const KmerPathLoc& loc1orig,
				  const KmerPathLoc& loc2orig,
				  vec<ImproperMerger>& ans,
				  const NegativeGapValidator* ngv ) {

  KmerPathLoc loc1=loc1orig, loc2=loc2orig;

  // If this scan doesn't reach a gap, we're done!
  if ( ScanRightPerfectMatch(loc1,loc2) == false  // kmer mismatch
       || loc1.atEnd() 
       || loc2.atEnd() ) {
    ans.resize(1);
    ans[0].right_end1 = loc1;
    ans[0].right_end2 = loc2;
    loc1.GetPath().CopySubpath( loc1orig, loc1, ans[0].merged );
    return;
  }

  // Otherwise, we must have reached a gap.  Call the gap skipper.

  KmerPath initial_perfect_match;
  loc1.GetPath().CopySubpath( loc1orig, loc1, initial_perfect_match );

  vec<FarEnd> far_ends;

  GapSkipRight( loc1, loc2, far_ends, ngv );

  // If there's no way to skip the gap, we're done!
  if( far_ends.empty() ) {
    ans.resize(1);
    ans[0].right_end1 = loc1;
    ans[0].right_end2 = loc2;
    ans[0].merged = initial_perfect_match;
    return;
  }

  // Recursively call self on each far_end; build answers in ans.
  vec<ImproperMerger> sub_answers;

  for( vec<FarEnd>::iterator far = far_ends.begin();
       far != far_ends.end(); far++ ) {
    if( far->DONE ) {  // there was not another shared kmer -- done!
      ans.resize( ans.size()+1 );
      ans.back().merged = initial_perfect_match;
      ans.back().right_end1 = loc1;
      ans.back().right_end2 = loc2;
    } else {
      sub_answers.clear();
      RecursiveImproperMergeRight( far->loc1, far->loc2, sub_answers, ngv );
      for( vec<ImproperMerger>::iterator sub_ans = sub_answers.begin();
	   sub_ans != sub_answers.end(); sub_ans++ ) {
	ans.resize( ans.size()+1 );
	// perfect match + merged gap section + recursive answer
	ans.back().merged = initial_perfect_match;
	ans.back().merged.Append( far->merged );
	ans.back().merged.AppendNoFirstKmer( sub_ans->merged );
	ans.back().right_end1 = sub_ans->right_end1;
	ans.back().right_end2 = sub_ans->right_end2;
      }
    }
  }
}


  
KmerPathLoc ImproperMerger::PushLoc1( const KmerPathLoc& loc1 ) const {

  // Get the matched sub-path of path 1:

  KmerPath sub_path;
  PathEmbedding into_read = SubpathEmbedding( loc1.GetPath(),
					      this->left_end1,
					      this->right_end1,
					      &sub_path );
  KmerPathLoc sub_loc = into_read.Preimage( loc1 );

  // Find the embedding(s) of the sub-path into the merger:
  vec<PathEmbedding> embeds;
  FindPathEmbeddings( sub_path, this->merged, embeds,
		      0, this->merged.NSegments()-1 );

  ForceAssert( embeds.size() > 0 );

  // It is possible that there are multiple embeddings.
  // Sometimes, with work, we could determine which one to use.
  // Sometimes this is impossible, it's inherently ambiguous.

  // For now, ignore all that and just pick one.
  KmerPathLoc ans = embeds[0].PushLoc( sub_loc );

  return ans;
}


KmerPathLoc ImproperMerger::PushLoc2( const KmerPathLoc& loc2 ) const {

  // Get the matched sub-path of path 2:

  KmerPath sub_path;
  PathEmbedding into_read = SubpathEmbedding( loc2.GetPath(),
					      this->left_end2,
					      this->right_end2,
					      &sub_path );
  KmerPathLoc sub_loc = into_read.Preimage( loc2 );

  // Find the embedding(s) of the sub-path into the merger:
  vec<PathEmbedding> embeds;
  FindPathEmbeddings( sub_path, this->merged, embeds,
		      0, this->merged.NSegments()-1 );

  ForceAssert( embeds.size() > 0 );

  // It is possible that there are multiple embeddings.
  // Sometimes, with work, we could determine which one to use.
  // Sometimes this is impossible, it's inherently ambiguous.

  // For now, ignore all that and just pick one.
  KmerPathLoc ans = embeds[0].PushLoc( sub_loc );

  return ans;
}



// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

#include "CoreTools.h"
#include "paths/GapSkip.h"

#include "paths/NegativeGapValidator.h"


// Things for gap-skipping:
  


// Just a wrapper for starting the potential recursion of GapOnGap
// and handling the flipping of locations and answers if necessary.
// Coming in, loc1 and loc2 point to the last matched kmers before a gap.
void GapSkipRight(KmerPathLoc loc1, KmerPathLoc loc2, 
		  vec<FarEnd>& far_ends, const NegativeGapValidator* ngv ) {

  // We need the gap to be on path1.
  Bool needs_flipping = ( ! loc1.GapToRight() );
  if( needs_flipping ) swap(loc1, loc2);
  // If this Assert ever happens, there's an algorithmic flaw in the caller:
  // this should only be called if one of the locs is next to a gap.
  ForceAssert( loc1.GapToRight() );

  vec< pair<int,int> > will_be_discarded;

  // Increment the locs to point to first non-matched kmer.
  // (can't just loc2.Increment*Gap(+1) -- would hop over gaps of length 0.
  loc1.IncrementInterval();  // Now points into the gap
  if( loc2.GapToRight() )
    loc2.IncrementInterval();  // Now points into the gap
  else
    loc2.IncrementHaltAtGap(+1);

  GapOnGapRight( loc1, loc2, 0, 0, far_ends, will_be_discarded, ngv );

  if( needs_flipping ) // un-flip the answers:
    for( vec<FarEnd>::iterator far = far_ends.begin(); 
	 far != far_ends.end(); far++ )
      far->flip();
}

void GapSkipLeft(KmerPathLoc loc1, KmerPathLoc loc2, 
		 vec<FarEnd>& far_ends, const NegativeGapValidator* ngv ) {

  // We need the gap to be on path1.
  Bool needs_flipping = ( ! loc1.GapToLeft() );
  if( needs_flipping ) swap(loc1, loc2);
  // If this Assert ever happens, there's an algorithmic flaw in the caller:
  // this should only be called if one of the locs is next to a gap.
  ForceAssert( loc1.GapToLeft() );

  vec< pair<int,int> > will_be_discarded;

  // Increment the locs to point to first non-matched kmer:
  // (can't just loc2.Increment*Gap(-1) -- would hop over gaps of length 0.
  loc1.DecrementInterval();  // Now points into the gap.
  if( loc2.GapToLeft() )
    loc2.DecrementInterval();  // Now points into the gap.
  else
    loc2.IncrementHaltAtGap(-1);

  GapOnGapLeft( loc1, loc2, 0, 0, far_ends, will_be_discarded, ngv );

  if( needs_flipping ) // un-flip the answers:
    for( vec<FarEnd>::iterator far = far_ends.begin(); 
	 far != far_ends.end(); far++ )
      far->flip();
}



// We have just processed a gap on loc2, and it can end
// somewhere in the middle of a gap on loc1, using up
// between min_gap_used and max_gap_used of the loc1 gap kmers.
// Continue aligning from there until:
//   case 1: the rest of read2 fits inside the gap on read1, or
//   case 2: we find a match on read2 to the kmer after the gap on read1, or
//   case 3: the gap on read1 might end in a gap on read2; recurse.
//
// Returns a FarEnd for each alignment it finds.
// Also for each, returns a pair of integers which should replace
// min_gap_used and max_gap_used in the caller, showing the least/most
// gap kmers which should be taken care of before the indicated
// alignment can begin.  (This is needed for setting surrounding gaps.)

// Coming in, loc2orig should point at the first kmer to be covered
// by un-accounted-for gap space on loc1, and loc1orig should point
// to the gap interval (only its index, not its loc, matters).

void GapOnGapRight( KmerPathLoc loc1orig,
		    KmerPathLoc loc2orig,
		    int min_gap_used, int max_gap_used,
		    vec<FarEnd>& far_ends,
		    vec< pair<int,int> >& real_gap_used,
		    const NegativeGapValidator* ngv ) {

  ForceAssert( loc1orig.isGap() );

  const KmerPath& path1 = loc1orig.GetPath();
  const KmerPath& path2 = loc2orig.GetPath();

  // If we have an NGV, use it on min_gap_used, if appropriate
  if( ngv != NULL && max_gap_used != 0  // =0 when starting from kmers
      && min_gap_used < ngv->GetK()-1 ) {

    const longlong kmer1 = path1.Stop(loc1orig.GetIndex()-1);
    const longlong kmer2 = loc2orig.GetKmer();
    min_gap_used = ngv->MinGap( kmer1, kmer2, min_gap_used );
    
    if( min_gap_used > max_gap_used ) return;  // no alignments possible!

    if( max_gap_used < ngv->GetK()-1 )
      max_gap_used = ngv->MaxGap( kmer1, kmer2, max_gap_used );
  }

  KmerPathLoc loc1 = loc1orig;
  loc1.IncrementInterval();  // so loc1 points to the first kmer after the gap

  // The fewest/most gap kmers available from the path1 gap:
  int min_gap_avail = max((longlong)0, loc1orig.Minimum() - max_gap_used);
  int max_gap_avail = loc1orig.Maximum() - min_gap_used;


  // Case 1: path2 may end inside the path1 gap.

  KmerPathLoc loc2end = path2.End();
  int ngv_buffer = ((ngv == NULL) ? 0 : ngv->MinGap( loc2end.GetKmer(),
						     loc1.GetKmer(), 0 ));

  if( DistMin(loc2orig, loc2end)+1 <= max_gap_avail - ngv_buffer ) {
    far_ends.push_back( FarEnd() );
    FarEnd& end_in_gap = far_ends.back();  // construct merged path in place
    end_in_gap.DONE = True;  // don't recurse on this

    // Tail of path2, with gaps adjusted:
    path2.CopySubpathAdjustGaps( loc2orig, loc2end, 
				 0, max_gap_avail - ngv_buffer,
				 end_in_gap.merged );
    // Remainder of the path1 gap:
    end_in_gap.merged.AddGap( max(ngv_buffer,
				  min_gap_avail-DistMax(loc2orig,loc2end)-1),
			      max_gap_avail - DistMin(loc2orig,loc2end)-1 );
    // Remainder of path1 after the gap:
    path1.CopyTail( loc1, end_in_gap.merged );
    // Set stop values:
    end_in_gap.stop1 = 0;
    end_in_gap.stop2 = path1.NSegments() - loc1orig.GetIndex();

    // push real gap_used interval
    real_gap_used.push_back(make_pair( min_gap_used,
				       min(max_gap_used, 
					   (int)loc1orig.Maximum() -
					   ngv_buffer -
					   DistMin(loc2orig, loc2end) - 1 )));
  }

  // If path2 *must* end inside the gap, we're done -- return now.
  if( DistMax(loc2orig,loc2end)+1 <= min_gap_avail )
    return;

  // Now scan through the places on read2 where the gap could end.

  KmerPathLoc loc2near = loc2orig, loc2far = loc2orig;
  loc2near.IncrementMaxGap( max(0, (int)loc1orig.Minimum() - max_gap_used) );
  loc2far.IncrementMinGap( loc1orig.Maximum() - min_gap_used );

  // these now point at possibilities for the first kmer past the loc1 gap
  // (no +1s: loc2 comes in pointing to the first kmer covered by the gap)
  longlong kmer1 = loc1.GetKmer();
  KmerPathLoc loc2 = loc2orig;

  for(int index = loc2near.GetIndex(); index <= loc2far.GetIndex(); index++) {
    loc2.SetIndex(index);

    // Case 2: we find a match for kmer1 on path2.

    if( loc2.isSeq() && loc2.GetSegment().Contains(kmer1) ) {
      // Deal with first/last interval cutoff:
      if( index==loc2near.GetIndex() && kmer1 < loc2near.GetKmer() ) continue;
      if( index==loc2far.GetIndex() && kmer1 > loc2far.GetKmer() ) continue;
      // Don't match to the kmer just after a gap -- will be caught in Case 3
      if( kmer1 == loc2.Start() && loc2.isGap(-1) 
	  && index != loc2near.GetIndex() ) continue;

      // Copy the part of path2 covered by the path1 gap.
      loc2.SetKmer( kmer1 );
      // loc2 now points at the successfully matched kmer.
      far_ends.push_back( FarEnd(loc1, loc2) );

      path2.CopySubpathAdjustGaps( loc2orig, loc2,
				   min_gap_avail+1, max_gap_avail+1,
				   far_ends.back().merged );
      // Those +1's since we're copying the kmer past the end of the gap too

      // Give feedback on how much of the gap we needed:
      int min_OK_used = loc1orig.Minimum() - DistMax(loc2orig,loc2);
      int max_OK_used = loc1orig.Maximum() - DistMin(loc2orig,loc2);
      real_gap_used.push_back( make_pair( max( min_gap_used, min_OK_used ),
					  min( max_gap_used, max_OK_used )));
    } // done with loc2.isSeq()


    // Case 3: the gap on path1 could end in a gap on path2.  Recurse.

    if( loc2.isGap() ) {
      vec<FarEnd> sub_far_ends;
      vec< pair<int,int> > sub_used;
      int sub_min = 0, sub_max = loc2.Maximum();
      // Deal with first/last interval cutoff:
      if( index == loc2near.GetIndex() ) sub_min = loc2near.GetLoc();
      if( index == loc2far.GetIndex() ) sub_max = loc2far.GetLoc();

      // How far into this gap can we cover?  WARNING: not necessarily
      // all of it, even if we can reach past its end!
      if( index != loc2far.GetIndex() ) {
	loc2.SetLoc(0);
	int max_gap_cover = max_gap_avail - DistMin(loc2orig, loc2);
	if( max_gap_cover < sub_max )
	  sub_max = max_gap_cover;
      }      

      // Recurse, with reads 1 and 2 flipped:
      GapOnGapRight( loc2, loc1, sub_min, sub_max,
		     sub_far_ends, sub_used, ngv );

      // Step through and process recursion results:
      for( unsigned int i=0; i<sub_far_ends.size(); i++ ) {

	far_ends.push_back( FarEnd() );
	FarEnd& far = far_ends.back();  // so I can build it in-place
	FarEnd& sub_far = sub_far_ends[i];
	// Copy most values, re-flipping reads 1 and 2 as needed:
	far.loc1 = sub_far.loc2;
	far.loc2 = sub_far.loc1;
	far.DONE = sub_far.DONE;
	far.stop1 = sub_far.stop2;
	far.stop2 = sub_far.stop1;

	// Construct the new merged:
	//  Segment of path2, with gaps adjusted (using recursion data!),
	//  up to but not including the gap on loc2 where the loc1 gap ends:
	KmerPathLoc loc2seq_end = loc2;
	loc2seq_end.DecrementInterval();  // points to end of seq before gap2
	path2.CopySubpathAdjustGaps( loc2orig, loc2seq_end,
				     max(0,min_gap_avail - sub_used[i].second),
				     max_gap_avail - sub_used[i].first,
				     far.merged );
	// The merged gap, where gap1 and gap2 end up overlapping:
	far.merged.AddGap(sub_used[i].first, sub_used[i].second);
	// The sub merged (no pun intended):
	far.merged.Append( sub_far.merged );

	// Give feedback on how much of the gap we needed:
	int min_OK_used = loc1orig.Minimum() 
	                  - DistMax(loc2orig,loc2seq_end)-1
	                  - sub_used[i].second;
	int max_OK_used = loc1orig.Maximum() 
	                  - DistMin(loc2orig,loc2seq_end)-1
	                  - sub_used[i].first;

	real_gap_used.push_back( make_pair( max( min_gap_used, min_OK_used ),
					    min( max_gap_used, max_OK_used )));
      }
    } // done with loc2.isGap()
  } // done stepping through indices from loc2near to loc2far.
}



void GapOnGapLeft( KmerPathLoc loc1orig,
		   KmerPathLoc loc2orig,
		   int min_gap_used, int max_gap_used,
		   vec<FarEnd>& far_ends,
		   vec< pair<int,int> >& real_gap_used,
		   const NegativeGapValidator* ngv ) {

  ForceAssert( loc1orig.isGap() );

  const KmerPath& path1 = loc1orig.GetPath();
  const KmerPath& path2 = loc2orig.GetPath();

  // If we have an NGV, use it on min_gap_used, if appropriate
  if( ngv != NULL && max_gap_used != 0  // =0 when starting from kmers
      && min_gap_used < ngv->GetK()-1 ) {

    const longlong kmer1 = path1.Start(loc1orig.GetIndex()+1);
    const longlong kmer2 = loc2orig.GetKmer();
    min_gap_used = ngv->MinGap( kmer2, kmer1, min_gap_used );
    
    if( min_gap_used > max_gap_used ) return;  // no alignments possible!

    if( max_gap_used < ngv->GetK()-1 )
      max_gap_used = ngv->MaxGap( kmer2, kmer1, max_gap_used );
  }

  KmerPathLoc loc1 = loc1orig;
  loc1.DecrementInterval();  // so loc1 points to the first kmer before the gap

  // The fewest/most gap kmers available from the path1 gap:
  int min_gap_avail = max((longlong)0, loc1orig.Minimum() - max_gap_used);
  int max_gap_avail = loc1orig.Maximum() - min_gap_used;


  // Case 1: path2 may end inside the path1 gap.

  KmerPathLoc loc2end = path2.Begin();
  int ngv_buffer = ((ngv == NULL) ? 0 : ngv->MinGap( loc1.GetKmer(),
						     loc2end.GetKmer(), 0 ));

  if( DistMin(loc2end, loc2orig)+1 <= max_gap_avail - ngv_buffer ) {
    far_ends.push_back( FarEnd() );
    FarEnd& end_in_gap = far_ends.back();  // construct merged path in place
    end_in_gap.DONE = True;  // don't recurse on this

    // part of path1 before the gap:
    path1.CopyHead( loc1, end_in_gap.merged );
    // Remainder of the path1 gap:
    end_in_gap.merged.AddGap( max(ngv_buffer,
				  min_gap_avail-DistMax(loc2end,loc2orig)-1),
			      max_gap_avail - DistMin(loc2end,loc2orig)-1 );
    // Tail of path2, with gaps adjusted:
    path2.CopySubpathAdjustGaps( loc2end, loc2orig, 
				 0, max_gap_avail - ngv_buffer ,
				 end_in_gap.merged );
    // Set stop values:
    end_in_gap.stop1 = 0;
    end_in_gap.stop2 = loc1orig.GetIndex() + 1;
    
    // push real gap_used interval
    real_gap_used.push_back(make_pair( min_gap_used,
				       min(max_gap_used, 
					   (int)loc1orig.Maximum() -
					   ngv_buffer -
					   DistMin(loc2end, loc2orig) - 1 )));
  }

  // If path2 *must* end inside the gap, we're done -- return now.
  if( DistMax(loc2end,loc2orig)+1 <= min_gap_avail )
    return;

  // Now scan through the places on read2 where the gap could end.

  KmerPathLoc loc2near = loc2orig, loc2far = loc2orig;
  loc2near.IncrementMaxGap( - max(0, (int)loc1orig.Minimum() - max_gap_used) );
  loc2far.IncrementMinGap( -(loc1orig.Maximum() - min_gap_used) );

  // these now point at possibilities for the first kmer past the loc1 gap
  // (no +1s: loc2 comes in pointing to the first kmer covered by the gap)
  longlong kmer1 = loc1.GetKmer();
  KmerPathLoc loc2 = loc2orig;

  for(int index = loc2near.GetIndex(); index >= loc2far.GetIndex(); index--) {
    loc2.SetIndex(index);

    // Case 2: we find a match for kmer1 on path2.

    if( loc2.isSeq() && loc2.GetSegment().Contains(kmer1) ) {

      // Deal with first/last interval cutoff:
      if( index==loc2near.GetIndex() && kmer1 > loc2near.GetKmer() ) continue;
      if( index==loc2far.GetIndex() && kmer1 < loc2far.GetKmer() ) continue;
      // Don't match to the kmer just before a gap -- will be caught in Case 3
      if( kmer1 == loc2.Stop() && loc2.isGap(+1) 
	  && index != loc2near.GetIndex() ) continue;

      // Copy the part of path2 covered by the path1 gap.
      loc2.SetKmer( kmer1 );
      // loc2 now points at the successfully matched kmer.
      far_ends.push_back( FarEnd(loc1, loc2) );

      path2.CopySubpathAdjustGaps( loc2, loc2orig,
				   min_gap_avail+1, max_gap_avail+1,
				   far_ends.back().merged );
      // Those +1's since we're copying the kmer past the end of the gap too

      // Give feedback on how much of the gap we needed:
      int min_OK_used = loc1orig.Minimum() - DistMax(loc2,loc2orig);
      int max_OK_used = loc1orig.Maximum() - DistMin(loc2,loc2orig);
      real_gap_used.push_back( make_pair( max( min_gap_used, min_OK_used ),
					  min( max_gap_used, max_OK_used )));
    } // done with loc2.isSeq()


    // Case 3: the gap on path1 could end in a gap on path2.  Recurse.

    if( loc2.isGap() ) {
      vec<FarEnd> sub_far_ends;
      vec< pair<int,int> > sub_used;
      int sub_min = 0, sub_max = loc2.Maximum();
      // Deal with first/last interval cutoff
      // These locs last moved *left*, so # gap-mers used is -GetLoc()-1.
      if( index == loc2near.GetIndex() ) sub_min = -loc2near.GetLoc()-1;
      if( index == loc2far.GetIndex() ) sub_max = -loc2far.GetLoc()-1;

      // How far into this gap can we cover?  WARNING: not necessarily
      // all of it, even if we can reach past its end!
      if( index != loc2far.GetIndex() ) {
	loc2.SetLoc(-1); // rightmost location in this gap
	int max_gap_cover = max_gap_avail - DistMin(loc2, loc2orig);
	if( max_gap_cover < sub_max )
	  sub_max = max_gap_cover;
      }      

      // Recurse, with reads 1 and 2 flipped:
      GapOnGapLeft( loc2, loc1, sub_min, sub_max,
		     sub_far_ends, sub_used, ngv );

      // Step through and process recursion results:
      for( unsigned int i=0; i<sub_far_ends.size(); i++ ) {

	far_ends.push_back( FarEnd() );
	FarEnd& far = far_ends.back();  // so I can build it in-place
	FarEnd& sub_far = sub_far_ends[i];
	// Copy most values, re-flipping reads 1 and 2 as needed:
	far.loc1 = sub_far.loc2;
	far.loc2 = sub_far.loc1;
	far.DONE = sub_far.DONE;
	far.stop1 = sub_far.stop2;
	far.stop2 = sub_far.stop1;

	// Construct the new merged:
	// The sub merged (no pun intended):
	far.merged.Append( sub_far.merged );
	// The merged gap, where gap1 and gap2 end up overlapping:
	far.merged.AddGap(sub_used[i].first, sub_used[i].second);
	//  Segment of path2, with gaps adjusted (using recursion data!),
	//  up to but not including the gap on loc2 where the loc1 gap ends:
	KmerPathLoc loc2seq_end = loc2;
	loc2seq_end.IncrementInterval();  // points to end of seq before gap2
	path2.CopySubpathAdjustGaps( loc2seq_end, loc2orig,
				     max(0,min_gap_avail - sub_used[i].second),
				     max_gap_avail - sub_used[i].first,
				     far.merged );

	// Give feedback on how much of the gap we needed:
	int min_OK_used = loc1orig.Minimum() 
	                  - DistMax(loc2seq_end,loc2orig)-1
	                  - sub_used[i].second;
	int max_OK_used = loc1orig.Maximum() 
	                  - DistMin(loc2seq_end,loc2orig)-1
	                  - sub_used[i].first;

	real_gap_used.push_back( make_pair( max( min_gap_used, min_OK_used ),
					    min( max_gap_used, max_OK_used )));
      }
    } // done with loc2.isGap()
  } // done stepping through indices from loc2near to loc2far.
}



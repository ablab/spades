// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

// A class to deal with K-mer read paths containing gaps of size <K.

// Note a bunch of +1s and -1s appear below because two kmers
// with a gap of 0 between them are already offset by 1.  (Thanks Jon.)

#include "paths/NegativeGapValidator.h"
#include "paths/AlignAndMerge.h"  // for MergedKmerPath struct, for FillCNGs

// Just say yes or no
bool NegativeGapValidator::Validate(const KmerPath& path, 
				    int seg_min, int seg_max) const {
  if( seg_max < 0 ) seg_max += path.NSegments(); // so eg -1 means last

  for(int seg = max(seg_min,1); 
      seg <= min(seg_max, path.NSegments()-2); seg++) {
    if( path.isSeq(seg) ) continue;
    if( path.Maximum(seg) >= K-1 ) continue;  // certainly OK with me

    // found a potentially actionable gap.  Check it out.
    int d = mp_kbb->MinOffset( path.Stop(seg-1), 
			       path.Start(seg+1), 
			       path.Minimum(seg)+1 ) - 1;

    if( verbosity > 1 )
      cout << path.Segment(seg) << " between kmers " 
	   << path.Stop(seg-1) << " and " << path.Start(seg+1)
	   << ", min offset " << d << ": ";

    if( d > path.Maximum(seg) ) { // No way this gap could exist!
      if( verbosity > 0 ) 
	cout << "rejecting path due to neg gap" << endl;
      return false;
    } else {
      if( verbosity > 1 )
	cout << "OK" << endl;
    }
  }

  // All gaps passed the test
  return true;
}



// Activist validator: restrict the limits of gaps if necessary
bool NegativeGapValidator::MakeValid(KmerPath& path,
				     int seg_min, int seg_max) const {
  if( seg_max < 0 ) seg_max += path.NSegments(); // so eg -1 means last

  for(int seg = max(seg_min,1); 
      seg <= min(seg_max, path.NSegments()-2); seg++) {
    if( path.isSeq(seg) ) continue;
    if( path.Minimum(seg) >= K-1 ) continue;  // can't touch these.

    // found a potentially actionable gap.  Check it out.

    int dmin = mp_kbb->MinOffset( path.Stop(seg-1), 
				  path.Start(seg+1), 
				  path.Minimum(seg)+1 ) - 1;
    int dmax = mp_kbb->MaxOffset( path.Stop(seg-1), 
				  path.Start(seg+1), 
				  path.Maximum(seg)+1 ) - 1;

    if( verbosity > 1 )
      cout << path.Segment(seg) << " between kmers " 
	   << path.Stop(seg-1) << " and " << path.Start(seg+1)
	   << ", real gap " << dmin << " to " << dmax << ": ";

    if( dmin > path.Maximum(seg) ) { // No way this gap could exist!
      if( verbosity > 0 ) 
	cout << "rejecting path due to negative gap" << endl;
      return false;
    }
    if( dmin > path.Minimum(seg) || dmax < path.Maximum(seg) ) { // correct it!
      path.SetGap(seg, dmin, dmax);
      if( verbosity > 0 )
	cout << "changed gap to " << path.Segment(seg) << endl;
    } else {
      if( verbosity > 1 )
	cout << "OK" << endl;
    }
  }

  // All gaps passed the test
  return true;
}

// If handed a non-NULL MergerKmerPath*, this updates the struct's
// left_end, right_end, and given fields.
bool NegativeGapValidator::FillConstantNegativeGaps( const KmerPath& path,
						     KmerPath& ans,
						     MergedKmerPath* mkp,
						     int seg_min, 
						     int seg_max ) const {
  ans.Clear();
  bool filled = false;

  if( seg_max < 0 ) seg_max += path.NSegments();

  KmerPath filling;

  // Data on how segments map:
  vec<int> seg_map(path.NSegments());

  for (int seg=0; seg < path.NSegments(); seg++) {
    if( seg_min <= seg && seg <= seg_max
	&& path.isGap(seg) 
	&& path.Stretch(seg)==0 
	&& path.Minimum(seg) < K
	&& mp_kbb->KmersBetween( path.Stop(seg-1), path.Start(seg+1),
				 path.Minimum(seg), filling ) ) {
      ans.Append(filling);
      filled = true;
    }
    else
      ans.AddSegment( path.Segment(seg) );

    seg_map[seg] = ans.NSegments()-1; // we never care about this value for gaps
  }
  
  if( mkp != NULL ) {
    mkp->left_end.first = seg_map[mkp->left_end.first];
    mkp->left_end.second = seg_map[mkp->left_end.second];
    mkp->given = seg_map[mkp->given];
    mkp->right_end.first = 
      ans.NSegments()-1 - seg_map[ path.NSegments()-1 - mkp->right_end.first ];
    mkp->right_end.second = 
      ans.NSegments()-1 - seg_map[ path.NSegments()-1 - mkp->right_end.second ];
  }

  return filled;
}


bool NegativeGapValidator::LegalConcat( const KmerPath& p, 
					const KmerPath& q ) const {

  if( p.IsEmpty() || q.IsEmpty() ) return true;
  if( p.NSegments()==1 && p.isGap(0) ) return true;
  if( q.NSegments()==1 && q.isGap(0) ) return true;

  int mingap=0, maxgap=0;

  KmerPathLoc pend = p.End();
  KmerPathLoc qbegin = q.Begin();

  if( pend.isGap() ) {
    mingap += pend.Minimum();
    maxgap += pend.Maximum();
    pend.DecrementInterval();
  }
  if( qbegin.isGap() ) {
    mingap += qbegin.Minimum();
    maxgap += qbegin.Maximum();
    qbegin.IncrementInterval();
  }
  
  return( maxgap >= K || 
	  MinGap( pend.GetKmer(), qbegin.GetKmer(), mingap ) <= maxgap );  

}





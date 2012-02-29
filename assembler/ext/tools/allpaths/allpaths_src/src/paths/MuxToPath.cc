/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////


#include "paths/MuxToPath.h"
#include "system/Assert.h"

void MuxToPath::ToPath( const vec<Mux>& muxes, 
			const KmerPath& path_so_far,
			KmerPath& ans ) const {
  if( muxes.empty() ) {
    ans = path_so_far;
    return;
  }

  ans.Clear();
  // Reserve the right number of segments for kmer_path.

  int numsegs = path_so_far.NSegments();
  vec<Mux>::const_iterator muxp = muxes.begin();
  if( path_so_far.IsEmpty() ) {
    numsegs += 
      muxp->GetPathId().GetPathPtr(*mp_pathsFw,*mp_pathsRc)->NSegments();
    muxp++;
  }
  
  for( ; muxp != muxes.end(); muxp++ )
    numsegs += muxp->GetSegment();

  ans.Reserve( numsegs );

  // Copy the snippets.  Can't call Extend, since we do them in
  // the efficient, backwards order.

  muxp = --muxes.end();
  const KmerPath* current_read_path_p;
  const KmerPath* next_read_path_p = 
    muxp->GetPathId().GetPathPtr(*mp_pathsFw,*mp_pathsRc);

  // First give ans its initial kmer, as we AppendNoFirstPath everything
  ans.AddSegment( next_read_path_p->Start(0), next_read_path_p->Start(0) );

  // Then copy the middle snippets
  while( muxp != muxes.begin() ) {
    current_read_path_p = next_read_path_p;
    KmerPathLoc joining_kmer( current_read_path_p, muxp->GetSegment() );
    muxp--;
    next_read_path_p = muxp->GetPathId().GetPathPtr(*mp_pathsFw,*mp_pathsRc);
    // ForceAssert that the extension makes sense.
    // If this fails, then the vec<Mux> argument is corrupt.
    ForceAssert(joining_kmer.GetSegment().Contains(next_read_path_p->Start(0)));
    joining_kmer.SetKmer( next_read_path_p->Start(0) );
    current_read_path_p->CopySubpathNoFirstKmer( current_read_path_p->Begin(),
						 joining_kmer,
						 ans );
  } // note: muxp = muxes.begin() is still used later

  // Then add the last snippet and path_so_far, as needed:
  if( path_so_far.IsEmpty() ) {
    // In this case, copy the entire read of the first mux:
    ans.AppendNoFirstKmer( *next_read_path_p );
  }
  else {
    // In this case, copy the snippet and path_so_far
    current_read_path_p = next_read_path_p;
    KmerPathLoc joining_kmer( current_read_path_p, muxp->GetSegment() );
    // ForceAssert that the extension makes sense.
    // If this fails, then the_mux does not represent a legal
    // extension of the argument path_so_far; ToPath called badly
    ForceAssert( joining_kmer.GetSegment().Contains(path_so_far.Start(0)) );
    joining_kmer.SetKmer( path_so_far.Start(0) );
    current_read_path_p->CopySubpathNoFirstKmer( current_read_path_p->Begin(),
						 joining_kmer,
						 ans );
    ans.AppendNoFirstKmer( path_so_far );
  }

}


void MuxToPath::ExtendByMux( const Mux& the_mux, 
			     const KmerPath& path_so_far,
			     KmerPath& ans ) const {

  const KmerPath& the_path = 
    *the_mux.GetPathId().GetPathPtr(*mp_pathsFw,*mp_pathsRc);

  if( path_so_far.IsEmpty() )
    ans = the_path;
  else {
    KmerPathLoc joining_kmer( the_path, the_mux.GetSegment() );
    // ForceAssert that the extension makes sense.
    // If this fails, then the_mux does not represent a legal
    // extension of the kmer_path; data corrupted upstream.
    ForceAssert( joining_kmer.GetSegment().Contains(path_so_far.Start(0)) );
    joining_kmer.SetKmer( path_so_far.Start(0) );

    ans.Clear();
    the_path.CopyHead( joining_kmer, ans );
    ans.AppendNoFirstKmer( path_so_far );
  }
}

// This returns true if there was a perfect match, in which case ans
// holds the extended path.  It returns false if the path did not
// match path_so_far on their overlap, in which case ans is empty.
// (So yes, the return value is redundant.)
bool MuxToPath::ExtendByKmersIfMatch( const OrientedKmerPathId& okpid,
				      const int numKmers,
				      const KmerPath& path_so_far,
				      KmerPath& ans ) const {

  const KmerPath& the_path = 
    *okpid.GetPathPtr(*mp_pathsFw,*mp_pathsRc);

  if( path_so_far.IsEmpty() ) {
    ans = the_path;
    return true;
  }
  else {
    ans.Clear();

    KmerPathLoc joining_kmer( the_path, 0, 0 );
    bool enough_kmers = joining_kmer.IncrementMinGap( numKmers );
    // IncrementMinGap, because the kmer counts in muxes are MinDist.
    if( ! enough_kmers )
      return false;

    // Check whether the joining kmer matches.
    // This should never happen if the okpid is a mux of
    // the final read of path_so_far, but could happen if
    // intermediate reads in the mux graph did not get used
    // -- for example, because they didn't match earlier parts 
    // of the path, or becuase they were declared bad.
    if( joining_kmer.isGap() ||
	joining_kmer.GetKmer() != path_so_far.Start(0) )
      return false;

    // Do the paths match perfectly to the right?
    KmerPathLoc loc1 = joining_kmer, loc2 = path_so_far.Begin();

    if( ScanRightPerfectMatchGaps( loc1, loc2 ) ) {
      // The paths match perfectly on their overlap.
      the_path.CopyHead( joining_kmer, ans );
      ans.AppendNoFirstKmer( path_so_far );
      // It's possible that the_path actually extends further right
      // than path_so_far.  But even if so, we want to trim the path
      // we're constructing at the right end of the opening read,
      // so don't use any kmers which stick out beyond it.
      return true;
    }
    else {
      // The paths disagree on their overlap.  This is possible
      // because of transitivity issues: if A->B->C are muxes,
      // there is no guarantee that A and C agree, off B.
      return false;
    }
  }
}


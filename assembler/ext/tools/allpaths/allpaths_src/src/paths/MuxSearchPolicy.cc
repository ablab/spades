/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
#include "paths/MuxSearchPolicy.h"

////////////////////////////////////////////////////////////////////////////
//
//  MSP_PairedPairs stuff:
//

// Internal helper function, thought at some point might want to
// extract this for broader use.
//
// Call this on two paths oriented the same direction.
// This just checks that the two paths have a perfect match and
// so do their partners.  It doesn't care at all about the separation
// and deviation of the pair, just the paths.
//
// If given a nonnegative segoff, it only checks the original
// reads with the first read extending segoff segments to the
// left of the second read.  (Muxes provide that information.)
// The partner reads may still align in any way.

bool MSP_PairedPairs::PairedPairsPerfectMatch( const OrientedKmerPathId& okpid1,
					       const OrientedKmerPathId& okpid2,
					       const int segoff ) const {

  if( okpid1.IsRc() != okpid2.IsRc() ) // point different directions
    return false;

  const int path_id1 = okpid1.GetId();
  const int path_id2 = okpid2.GetId();
  const int partner_id1 = mp_pairs_handler->GetPartnerId( path_id1 );
  const int partner_id2 = mp_pairs_handler->GetPartnerId( path_id2 );

  if( partner_id1<0 || partner_id2<0 )  // unpaired read
    return false;

  if( m_verbosity > 4 ) {
    cout << "  Reads " << path_id1 << " and " << path_id2 
	 << ( HavePerfectMatch( (*mp_pathsFw)[path_id1], 
				(*mp_pathsFw)[path_id2] ) ?
	      " DO " : " do NOT " )
	 << " match" << endl;
    cout << "  Partners " << partner_id1 << " and " << partner_id2 
	 << ( HavePerfectMatch( (*mp_pathsFw)[partner_id1], 
				(*mp_pathsFw)[partner_id2] ) ?
	      " DO " : " do NOT " )
	 << " match" << endl;
  }

  if( segoff >= 0 ) {
    // Check only one offset for path_id1 vs path_id2
    const KmerPath* path1 = okpid1.GetPathPtr( *mp_pathsFw, *mp_pathsRc );
    if( path1->NSegments() <= segoff ) return false; // but caller is dumb!

    const KmerPath* path2 = okpid2.GetPathPtr( *mp_pathsFw, *mp_pathsRc );

    KmerPathLoc loc1(path1, segoff), loc2(path2,0);

    return( path1->Segment(segoff).Contains(path2->Start(0))
	    && IsPerfectMatchToRight( loc1, loc2 )
	    && HavePerfectMatch( (*mp_pathsFw)[partner_id1], 
				 (*mp_pathsFw)[partner_id2] ) );
  }
  else {
    // Check all offsets for both paths
    return( HavePerfectMatch( (*mp_pathsFw)[path_id1], 
			      (*mp_pathsFw)[path_id2] ) &&
	    HavePerfectMatch( (*mp_pathsFw)[partner_id1], 
			      (*mp_pathsFw)[partner_id2] ) );
  }

}


MuxSearchPolicy::Opinion
MSP_PairedPairs::PushAndGiveOpinion( const Mux& final_mux, 
				     const MuxSearchState& search,
				     const vec<Bool>& is_good ) {
  // This MSP maintains no internal state.
  // So all we do is test whether this read is good.

  OrientedKmerPathId final_read = final_mux.GetPathId();
  const int partner_id = mp_pairs_handler->GetPartnerId( final_read.GetId() );
  if( partner_id < 0 )
    return( BAD );  // BAD: unpaired read

  // set up:
  int segs_back = final_mux.GetSegment();
  int kmers_back = final_mux.GetNumKmers();
  const int max_kmers_back = ((*mp_pathsFw)[ final_read.GetId() ]).KmerCount() - 1;

  // We need some criterion for declaring the beginning muxes all happy.
  // What should it be?  Certainly not clear to me.
  // For now, here's an expedient answer: any read which reaches back
  // to dist_to_here=0 (leftmost kmer of opener) is automatically happy.
  //
  // This makes some sense if the opener is as short as a paired pair end read.
  // With longer openers or pairs for openers, we'd need to rethink.
  if( search.dist_stack.empty() ||
      max_kmers_back >= search.dist_stack.back() + kmers_back )
    return( GOOD );  // GOOD: initial condition

  // Walk back through path_to_here, looking for happy same-orientation reads
  vec<Mux>::const_reverse_iterator muxp;
  vec<Bool>::const_reverse_iterator op;

  for( muxp = search.path_to_here.rbegin(), op = is_good.rbegin();
       muxp != search.path_to_here.rend() && kmers_back <= max_kmers_back; 
       segs_back += muxp->GetSegment(), kmers_back += muxp->GetNumKmers(), 
	 muxp++, op++ ) {
    if( *op != true ) continue;  // only look at happy reads
    
    if( PairedPairsPerfectMatch( final_read, muxp->GetPathId(), segs_back ) )
      return( GOOD );  // GOOD: Paired-pairs matched an earlier happy pair.
  }

  return( BAD );  // BAD: no paired-pair match.

  // As it currently stands, there is no communication between fw and rc
  // reads.  If one flavor has no happy reads for a while, it can never
  // regain them.  And there's no check that partners of happy are happy
  // if they also appear in the mux list.  Lots of room for fiddling.
}



MuxSearchPolicy::Opinion
MSP_TrackPath::PushAndGiveOpinion( const Mux& final_mux, 
				   const MuxSearchState& search,
				   const vec<Bool>& is_good ) {

  m_kmerpaths.push_back( KmerPath() );
  m_pathifier.ExtendByMux( final_mux, *(m_kmerpaths.end()-2), 
			   m_kmerpaths.back() );
  return GOOD;
}


// This is very inefficient.  If time is actually a problem,
// rewrite this: could track whether the current path
// is a match (if it's not for our parent, also not for us),
// and could check whether a mux still matches more efficiently
// than by building the whole extended path.
MuxSearchPolicy::Opinion
MSP_FateOfClosure::PushAndGiveOpinion( const Mux& final_mux, 
				       const MuxSearchState& search,
				       const vec<Bool>& is_good ) {
  MSP_TrackPath::PushAndGiveOpinion( final_mux, search, is_good );
  
  if( CurrentPath().IsEmpty() )
    return GOOD;

  // Does the current path match the closure so far?
  if( IsPerfectMatch( CurrentPath().End(), m_closure.End() ) ) {

    // Is this the closure?
    if( m_closure == CurrentPath() ) {
      cout << "Found the closure!" << endl;
      return GOOD;
    }

    // If not, are any of its muxes are also perfect matches?
    mp_muxGraph->GetMuxesOf( final_mux.GetPathId(), muxes );

    KmerPath extended_path;
    for( vec<Mux>::iterator mux = muxes.begin(); mux != muxes.end(); mux++ ) {

      m_pathifier.ExtendByMux( *mux, CurrentPath(), extended_path );

      if( IsPerfectMatch( extended_path.End(), m_closure.End() ) ) {
	// ...if so, all is still going well.  Keep quiet, unless verbose:
	if( m_verbose ) {
	  cout << "  So far so good: mux " << final_mux << "->" << *mux 
	       << " extending " << CurrentPath() << endl;
	}
	return GOOD;
      }
    }

    // Otherwise, the closure would need a missing mux.  Oh no!
    cout << "MSP_FateOfClosure: Mux search is diverging from closure!"
	 << "\nClosure: " << m_closure
	 << "\nCurrent: " << CurrentPath()
	 << "\nFinal read " << final_mux.GetPathId()
	 << " has " << muxes.size() 
	 << (muxes.size() == 1 ? " mux:" : " muxes:");
    for( vec<Mux>::iterator mux = muxes.begin(); mux != muxes.end(); mux++ )
      cout << "\n   " << *mux << "  " 
	   << *(mux->GetPathId().GetPathPtr(*mp_pathsFw,*mp_pathsRc));
    cout << "\nNo muxes match closure!" << endl;
  }

  return GOOD;
}




void MSP_ReachCloserOrAbort::Setup( int minExtLen, int maxExtLen,
				    const vec<OrientedKmerPathId>& closers,
				    const map< OrientedKmerPathId,
				               vec<SubsumptionRecord> >& 
				      closers_I_subsume ) {
  m_endpoints.clear();

  // As currently set up, we don't get anything here: the closers
  // themselves are "real" reads, with low id's, and we find them
  // through the unipath-extended super-reads subsuming them.
  int id;
  for( vec<OrientedKmerPathId>::const_iterator cl = closers.begin();
       cl != closers.end(); cl++ )
    if( UseThisRead(id = cl->GetId()) )
      m_endpoints.push_back( endpoint(id, minExtLen, maxExtLen) );

  // Here is where we pick up the subsumed closers.
  for( map< OrientedKmerPathId,vec<SubsumptionRecord> >::const_iterator
	 mov = closers_I_subsume.begin(); 
       mov != closers_I_subsume.end(); mov++ ) {
    for( vec<SubsumptionRecord>::const_iterator sub = mov->second.begin();
	 sub != mov->second.end(); sub++ ) {
      if( UseThisRead(id = sub->GetSuperPathId().GetId()) ) {
	int hang = sub->GetLeftOverhang();
	m_endpoints.push_back(endpoint( id, minExtLen+hang, maxExtLen+hang ));
      }
    } // end loop through vec<SubsumptionRecord>
  } // end loop through map

}


MuxSearchPolicy::Opinion
MSP_ReachCloserOrAbort::PushAndGiveOpinion( const Mux& final_mux, 
					    const MuxSearchState& search,
					    const vec<Bool>& is_good ) {
//  cout << "PAGO called on " << final_mux << endl;

  // Is this the right way to use an OffsetTracker??
  int here = final_mux.GetPathId().GetId();
  if( ! UseThisRead(here) )
    return GOOD;

  const ImplicitOffsetVec& can_get_here = mp_offsetTracker->GetOffsetsTo(here);

//  cout << "GetOffsetsTo returns " << can_get_here.size() << " offsets" << endl;

  int new_dist = 
    (search.dist_stack.empty() ? 0 : search.dist_stack.back() ) 
    + final_mux.GetNumKmers();

  for( vec<endpoint>::const_iterator endp = m_endpoints.begin();
       endp != m_endpoints.end(); endp++ )
    if(IsReachable( *endp, can_get_here, new_dist))
      return GOOD;

//  cout << "None of them are helpful" << endl;

  return PRUNE;

}

bool MSP_ReachCloserOrAbort::IsReachable( const endpoint& goal,
                                          const ImplicitOffsetVec& can_get_here,
					  int dist_so_far ) {
  const ImplicitOffset& io =
    *lower_bound( can_get_here.begin(), can_get_here.end(),
		  ImplicitOffset( goal.path_id, goal.min_dist - dist_so_far ) );

//   if( io.GetFrom() == goal.path_id && 
//       io.GetAmount() <= goal.max_dist - dist_so_far ) {
//     cout << "Endpoint [read " << goal.path_id 
// 	 << " at depth " << goal.min_dist << "-" << goal.max_dist
// 	 << "] matches offset [read " << io.GetFrom()
// 	 << " in " << io.GetAmount() << " kmers]" << endl;
//     return true;
//   }
//   else {
//     cout << "Endpoint [read " << goal.path_id 
// 	 << " at depth " << goal.min_dist << "-" << goal.max_dist
// 	 << "] doesn't match lower bound [read " << io.GetFrom()
// 	 << " in " << io.GetAmount() << " kmers]" << endl;
//     return false;
//   }

  return( io.GetFrom() == goal.path_id && 
	  io.GetAmount() <= goal.max_dist - dist_so_far );

}


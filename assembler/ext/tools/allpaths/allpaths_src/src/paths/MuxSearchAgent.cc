/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////


#include "paths/MuxSearchAgent.h"
#include <algorithm>  // for is_sorted

////////////////////////////////////////////////////////
//
// /////////////////////////////////////////////////////
// //
// //  Stuff for the base class:
// //



/////////////////////////////
//
//  Public member functions, called by KmerPathMuxSearcher::SearchDirector
//

// The MuxGraph encodes the minimal *leftwards* extensions of each read.
// Therefore we start searching on the right end of an insert and work
// our way to the left, opposite the sense of K.P.{Breadth,Depth}S.

void MuxSearchAgentBase::Setup( int openingReadId, 
				int closingReadId,
				int minLength,
				int maxLength,
				MuxSearchResult& result ) {
  search.clear();
  result.clear();
  dist_to_here = 0;

  mp_result = &result;

  minAcceptableExtLength = minLength;
  maxAcceptableExtLength = maxLength;

  rightReadId = openingReadId;
  leftReadId = closingReadId;

  // Seed search.to_explore with all the fills of the right read.
  SetupOpeners();

  // closers gets fills of closing read;
  // closers_I_subsume records any reads which subsume those fills.
  SetupClosers();

}



void MuxSearchAgentBase::SaveIfConverged() {
  vec<MuxWalkGraph::Label> graph_hits;
  mp_result->WalkGraph().FindNodes( final_mux.GetPathId(), graph_hits );

  if( ! graph_hits.empty() )
    // This doesn't mean we found a new closure for sure; still need
    // to check whether any hits yeild an acceptable length path.
    MergeIntoExistingClosures( graph_hits );
}



// This pushes final_mux onto the stack, but we don't push its muxes
// onto to_explore if final_mux goes beyond the end of the insert.
// The reason to push final_mux in this case is to check it for
// subsumed closures.

void MuxSearchAgentBase::PushMux() {

  search.path_to_here.push_back( final_mux );
  dist_to_here += final_mux.GetNumKmers();
  search.dist_stack.push_back( dist_to_here );
  // push its muxes onto to_explore, unless the path is already too long.
  // (But we do need to check too-long paths for mux-skipped closure.)
  search.to_explore.push_back( vec<Mux>() );
  if( dist_to_here <= maxAcceptableExtLength )
    mp_muxGraph->GetMuxesOf( final_mux.GetPathId(), search.to_explore.back() );
  
  if( m_verbosity > 3 ) {
    cout << "Pushing to explore mux ->" << final_mux.GetPathId() << endl;
    if( m_verbosity > 4 ) {
      if( dist_to_here <= maxAcceptableExtLength ) {
	cout << "GetEdgesFrom says " << final_mux.GetPathId()
	     << " has " << search.to_explore.back().size() 
	     << (search.to_explore.back().size() == 1 
		 ? " extension" : " extensions");
	for( uint i=0; i < search.to_explore.back().size(); i++ )
	  cout << " " << search.to_explore.back()[i];
	cout << endl;
      }
      else
	cout << "Not asking for extensions of " << final_mux.GetPathId()
	     << "; path is already too long" << endl;
    }
  }
}


// This may want modification to make it easier for derived classes to fiddle with.
void MuxSearchAgentBase::SaveIfClosed() {
  // Is it a closer?
  bool just_added_closer = false;
  if( binary_search( closers.begin(), closers.end(), final_mux.GetPathId() ) ) {
    const int closer_trim = mp_readfillDB->LeftTrimOfFill(final_mux.GetPathId().GetId());
    if( minAcceptableExtLength <= dist_to_here + closer_trim
	&& dist_to_here + closer_trim <= maxAcceptableExtLength ) {
      AddClosure( closer_trim );
      just_added_closer = true;
    }
  }

  if( ! just_added_closer
      && final_mux_subsumes_something 
      && minAcceptableExtLength <= dist_to_here ) {
    // We may have accidentally mux-skip a closer.
    // This can happen in the following (relatively unusual) case:
    //
    // old path:                   --------------------...
    // long mux:        -----------------
    // subsumed closer:    ---
    AddMuxSkippedClosures( subsume_iter->second );
  }
}



/////////////////////////////
//
//  Private member functions:
//

void MuxSearchAgentBase::SetupOpeners() {

  search.to_explore.resize( 1, vec<Mux>() );

  for( int f = mp_readfillDB->FirstFilling(rightReadId); 
       f <= mp_readfillDB->LastFilling(rightReadId); f++ ) {
    Mux this_opener( OrientedKmerPathId(f, true /* RC */),
		     -1,  // segments field -- irrelevant
		     -mp_readfillDB->RightTrimOfFill(f) // neg offset relative to orig read
		     );
    search.to_explore[0].push_back( this_opener );
  }

  if( m_verbosity>2 ) {
    cout << "Seeded search stack with " 
	 << mp_readfillDB->NumFillings(rightReadId)
	 << " starting fills" << endl;
  }
}


void MuxSearchAgentBase::SetupClosers() {

  closers.clear();
  closers_I_subsume.clear();

  for( int f = mp_readfillDB->FirstFilling(leftReadId); 
       f <= mp_readfillDB->LastFilling(leftReadId); f++ ) {
    // The genuine closer:
    OrientedKmerPathId this_closer(f, false /* FW */);
    closers.push_back( this_closer );
    // Reads which subsume it:
    vec<SubsumptionRecord> this_closer_subs;
    mp_subList->GetFullRecordsFor( this_closer, this_closer_subs );
    for( vec<SubsumptionRecord>::iterator sub_rec = this_closer_subs.begin();
	 sub_rec != this_closer_subs.end(); sub_rec++ ) {
      closers_I_subsume[sub_rec->GetSuperPathId()].push_back(*sub_rec);
    }

    if( m_verbosity > 5 ) {
      cout << "closer " << this_closer;
      if( ! this_closer_subs.empty() ) {
	cout << " subsumed by ";
	for( vec<SubsumptionRecord>::iterator sub_rec = this_closer_subs.begin();
	     sub_rec != this_closer_subs.end(); sub_rec++ )
	  cout << sub_rec->GetSuperPathId() << " ";
      }
      cout << endl;
    }
  }
  ForceAssert( is_sorted( closers.begin(), closers.end() ) );  

  if( m_verbosity>2 ) {
    cout << "Seeded with " << closers.size() << " closers" << endl;
  }
}




void MuxSearchAgentBase::MergeIntoExistingClosures( vec<MuxWalkGraph::Label>& graph_hits ) {

  if( m_verbosity > 6 )
    cout << "Read " << final_mux.GetPathId() << " has "
	 << graph_hits.size() << " hits in the walk graph" << endl;
  
  for( vec<MuxWalkGraph::Label>::iterator hitp = graph_hits.begin();
       hitp != graph_hits.end(); hitp++ ) {

    // do the lengths add up to something acceptable?
    int total_length = 
      dist_to_here + final_mux.GetNumKmers() + mp_result->WalkGraph().DistToEnd(*hitp);

    if( minAcceptableExtLength <= total_length &&
	total_length <= maxAcceptableExtLength ) {

      // Hooray! Merge this path into the walk graph.
      search.path_to_here.push_back( final_mux );
      mp_result->WalkGraph().MergeClosure( search.path_to_here, *hitp );
      search.path_to_here.pop_back();
      mp_result->num_closures_found++;
      if( m_verbosity > 1 ) {
	cout << "Merged into an existing closure!  At read "
	     << final_mux.GetPathId()
	     << ", depth " << search.dist_stack.back() + final_mux.GetNumKmers()
	     << " kmers" << endl;
	if( m_verbosity > 6 )
	  cout << "During exploration of state " 
	       << 1+mp_result->num_states_explored << endl;
      }
    }
  } // done checking final_mux for closure potential
}




void MuxSearchAgentBase::AddClosure( int closer_trim ) { 
  mp_result->WalkGraph().AddClosure( search.path_to_here, closer_trim );
  mp_result->num_closures_found++;
  if( m_verbosity > 1 ) {
    cout << "Found a closure!  At depth " 
	 << search.dist_stack.back() + closer_trim
	 << " kmers" << endl;
    if( m_verbosity > 6 )
      cout << "During exploration of state " 
	   << 1+mp_result->num_states_explored << endl;
  }
}




// Find closers hopped over by mux searching.
// This can happen in the following (relatively unusual) case:
//
// old path:                   --------------------...
// long mux:        -----------------
// subsumed closer:    ---
//
// Hmm: this has no notion of good/bad reads, and makes no reference to MayIClose().
// It's not clear what we want to do in this case.  Need to think about this more.

void MuxSearchAgentBase::AddMuxSkippedClosures( vec<SubsumptionRecord>& subsumptions ) {

  for( vec<SubsumptionRecord>::iterator sub_rec = subsumptions.begin();
       sub_rec != subsumptions.end(); sub_rec++ ) {
    // check that the subsumed closing read:
    // (1) would yeild an acceptable-length closure, and
    // (2) doesn't appear in exactly that position in path_to_here already!
    const int inferred_length = search.dist_stack.back() - sub_rec->GetLeftOverhang()
      + mp_readfillDB->LeftTrimOfFill(sub_rec->GetSubPathId().GetId());
    if( inferred_length < minAcceptableExtLength ||
	inferred_length > maxAcceptableExtLength )
      continue;

    int overhang = sub_rec->GetLeftOverhang();
    bool closer_present = false;
    for( vec<Mux>::reverse_iterator lookback = search.path_to_here.rbegin();
	 lookback != search.path_to_here.rend() && overhang >= 0; 
	 overhang -= lookback->GetNumKmers(), lookback++ ) {
      if( overhang==0 && lookback->GetPathId() == sub_rec->GetSubPathId() ) {
	closer_present = true;
	break;
      }
    }
    if( closer_present ) continue;

    // Hooray!  Found a mux-skipped closer!
    mp_result->WalkGraph().AddClosure( search.path_to_here, 
				  inferred_length - search.dist_stack.back() );
    mp_result->num_closures_found++;
    if( m_verbosity > 1 )
      cout << "Found a closure with subsumed end read!  "
	   << search.path_to_here.back().GetPathId() << " subsumes "
	   << sub_rec->GetSubPathId() << " at depth "
	   << inferred_length << " kmers" << endl;
  }
}


////////////////////////////////////////////////////////////////////////////
//
//  Display helpers used at high verbosity
//

void MuxSearchAgentBase::ReportSearchState() const {
  // Print the size of each level's exploration stack
  cout << "d=" << dist_to_here;
  if( !search.path_to_here.empty() ) {
    cout << ", ends with ";
    if( search.path_to_here.size() > 1 )
      cout << (search.path_to_here.end()-2)->GetPathId() << "->";
    cout << search.path_to_here.back().GetPathId();
  }
  cout << ", stack=";
  PrintExploreSizes( search.to_explore );
  cout << "." << endl;
}



void MuxSearchAgentBase::PrintExploreSizes( const vec< vec<Mux> >& to_explore ) const {
  int zeros=0, noptions;
  for(uint i=0; i<to_explore.size(); i++) {
    noptions = to_explore[i].size();
    if ( noptions == 0 )
      zeros++;
    else {
      if( zeros < 5 )
	for(;zeros>0;zeros--) cout << "0";
      else
	cout << "(0^" << zeros << ")";
      zeros=0;
      if( noptions > 9 )
	cout << " " << noptions << " ";
      else
	cout << noptions;
    }
  }
  if( zeros < 6 )
    for(;zeros>0;zeros--) cout << "0";
  else
    cout << "(0^" << zeros << ")";
}







////////////////////////////////////////////////////////
//
// /////////////////////////////////////////////////////
// //
// //  Stuff for derived classes:
// //


bool MuxSearchAgentGoodReads::MayIPush() {

  const Mux& final_mux = GetFinalMux();
  const int new_dist = GetDist() + final_mux.GetNumKmers();
  const int last_good_dist = m_lastGoodDist.empty() ? 0 : m_lastGoodDist.back();

  // If this mux sticks out so far that we could never reach
  // back to the last good read, we should prune immediately.
  if( last_good_dist < new_dist - m_max_read_length + 1 ) {
    if( m_verbosity > 2 )
      cout << "NOT pushing mux ->" << final_mux.GetPathId()
	   << ": good read coverage says to prune" << endl;
    return( false );
  }

  // Otherwise, we plan to push.  Tell/query the policies.

  bool this_is_good = true;
  bool prune = false;

  MuxSearchPolicy::Opinion opin;
  for( uint i=0; i < m_policies.size(); i++ ) {
    opin = m_policies[i]->PushAndGiveOpinion( final_mux,
					      GetSearch(),
					      m_isGood );
    this_is_good &= (opin == MuxSearchPolicy::GOOD);
    prune |= (opin== MuxSearchPolicy::PRUNE);
  }

  if( prune ) {
    // Some policy said we should prune this branch immediately.
    if( m_verbosity > 3 )
      cout << "  Some policy said to PRUNE this search branch!" << endl;
    PopPolicies();
    return( false );
  }
  // Otherwise, we are going to push the mux for sure.

  m_mayClose = this_is_good;

  // Also check for min_perfect_match
  // A read without a long enough overlap is deemed bad,
  // but if it's the closer, we allow its use anyway.
  this_is_good &=
    new_dist - (*mp_readLengths)[final_mux.GetPathId().GetId()] 
    <= last_good_dist - m_min_perfect_match;

  m_isGood.push_back( this_is_good );
  
  // m_lastGoodDist keeps track of where solid coverage by good reads ends.
  // If this read is good, does it reach back far enough?

  const int this_read_len = (*mp_readLengths)[final_mux.GetPathId().GetId()];
  const bool this_reaches = ( last_good_dist >= new_dist - this_read_len + 1 );
  
  if( this_is_good && this_reaches )
    m_lastGoodDist.push_back( new_dist );
  else
    m_lastGoodDist.push_back( last_good_dist );

  if( m_verbosity > 3 ) {
    cout << "  Good reads tracking says: " 
	 << final_mux << " is "
	 << (this_is_good ? "GOOD" : "BAD");
    if( this_is_good && ! this_reaches )
      cout << ", BUT not long enough for good coverage";
    cout << endl;
  }

  return( true );
}

// The hash functions used below are based on iterating h -> h*c + d
// for a fixed constant c and each new datum d.  Not deep, but good enough.
// Pick large odd c with large order mod 2^64 (almost any will do).

#define MUX_HASH_CONST1  123456789
#define MUX_HASH_CONST2  987654321

longlong MuxSearchAgentGoodReads::GoodReadsHash( const MuxSearchState& search,
						 const int dist_to_here,
						 const int min_dist) {
  longlong good_reads_hash = 0;

  int j=0, k=0;

  for( int i=search.path_to_here.size()-1; 
       i>=0 && search.dist_stack[i] >= min_dist; i-- ) {
    k += search.path_to_here[i].GetNumKmers();
    if( m_isGood[i] )
      good_reads_hash = ( good_reads_hash * MUX_HASH_CONST1 )
	+ ( longlong( search.path_to_here[i].GetPathId().GetIdRc() ) << 32 )
	+ ( k << 8 ) + j++;
  }

  return good_reads_hash;
}


pair<longlong,longlong> MuxSearchAgentGoodReads::PolicyHashOfPosition() {
  
  pair<longlong,longlong> here;
  
  here.first = BasicHashOfPosition();
  here.second = 0;

  const MuxSearchState& search = GetSearch();
  const int dist_to_here = GetDist();
  const int min_dist = dist_to_here - m_max_read_length;

  // Now give all policies a chance to contribute to the hash.
  //
  // Policies which don't track any internal state will return 0,
  // and multiplication by any odd number is invertible mod 2^n,
  // so this may be a wasted op but it causes no hash collisions.
  // 
  // A return value of 1 will be replaced by the good_reads_hash.
  //
  // Any other return value will be incorporated into the hash.

  longlong hashval;
  for( uint i=0; i < m_policies.size(); i++ ) {
    hashval = m_policies[i]->HashOfPosition( search, dist_to_here,
					     GetFinalMux(), m_isGood );
    if( hashval == 1 )
      hashval = GoodReadsHash( search, dist_to_here, min_dist );

    here.second = (here.second * MUX_HASH_CONST2) + hashval;
  }

  if( m_verbosity > 6 )
    cout << "state hash " 
	 << hex << here.first 
	 << ":" << here.second 
	 << dec << endl;

  return here;
}


// Pass vec<Bool> of which reads are good to the walk graph.
// Also, some extra reporting at high verbosity.
void MuxSearchAgentGoodReads::AddClosure( int closer_trim ) { 
  const MuxSearchState& search = GetSearch();

  // Somtimes a read is considered bad for walking purposes, but
  // good for closure purposes.  If that's the case, we should label
  // it good in the result graph, so nothing downstream is confused.
  bool last_is_good = m_isGood.back();
  m_isGood.back() = true;
  mp_result->WalkGraph().AddClosure( search.path_to_here, closer_trim, &m_isGood );
  m_isGood.back() = last_is_good;

  mp_result->num_closures_found++;

  if( m_verbosity > 1 ) {
    cout << "Found a closure!  At depth " 
	 << search.dist_stack.back() + closer_trim
	 << " kmers" << endl;
    if( m_verbosity > 4 ) {
      // Print the series of good muxes
      cout << "    Good path:\n    opener: " 
	   << search.path_to_here[0].GetPathId()
	   << ( m_isGood[0] ? " (good)" : " (BAD?!)" ) << endl;
      for(uint i=1; i < m_isGood.size(); i++)
	if( m_isGood[i] || i==m_isGood.size()-1 ) {
	  int skip = search.dist_stack[i] - m_lastGoodDist[i-1];
	  cout << "    -> " << search.path_to_here[i].GetPathId()
	       << " at d=" << search.dist_stack[i]
	       << " (skip " << skip
	       << ", read len " << (*mp_readLengths)[i]
	       << ( (skip < (*mp_readLengths)[i]) ? ")" : " (?!) )" )
	       << endl;
	}
      if( m_isGood.back() == false )
	cout << "    Note: final read is bad for walking, but a good closer" << endl;
    }
  }
}


// MergeIntoExistingClosures is harder in the good reads case:
// before we can merge, we need to know whether this read would
// be considered good or bad (which is recorded in the result).
// So we wrap this in a full push/pop.  This is slightly more
// work than is needed -- we have no use for the muxes of the
// pushed final_mux -- but it's easy, robust, and infrequent.

void MuxSearchAgentGoodReads::MergeIntoExistingClosures( vec<MuxWalkGraph::Label>& graph_hits ) {

  const Mux& final_mux = GetFinalMux();
  const int dist_to_here = GetDist();
  const MuxSearchState& search = GetSearch() ;

  if( m_verbosity > 6 )
    cout << "Read " << final_mux.GetPathId() << " has "
	 << graph_hits.size() << " hits in the walk graph" << endl;

  // find opinions and push this mux
  // If MayIPush() says no, there wouldn't have been any closures anyway.
  if( this->MayIPush() ) {
    this->PushMux();

    for( vec<MuxWalkGraph::Label>::iterator hitp = graph_hits.begin();
	 hitp != graph_hits.end(); hitp++ ) {

      // Does this hit agree on goodness of the current read?
      if( m_isGood.back() != mp_result->WalkGraph().IsGood(*hitp) )
	continue;

      // do the lengths add up to something acceptable?
      int total_length = 
	dist_to_here + final_mux.GetNumKmers() + mp_result->WalkGraph().DistToEnd(*hitp);

      if( total_length < GetMinAcceptableExtLength() || 
	  total_length > GetMaxAcceptableExtLength() )
	continue;

      // No need to push anything here; already done.
      mp_result->WalkGraph().MergeClosure( search.path_to_here, *hitp, &m_isGood );
      mp_result->num_closures_found++;
      if( m_verbosity > 1 )
	cout << "Merged into an existing closure!  At read "
	     << final_mux.GetPathId()
	     << ", depth " << search.dist_stack.back() + final_mux.GetNumKmers()
	     << " kmers" << endl;
    } // end loop over graph_hits

    this->Pop();
  }
}

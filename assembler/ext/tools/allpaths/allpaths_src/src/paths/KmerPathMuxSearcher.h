// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

#ifndef PATHS_KMERPATHMUXSEARCHER_H
#define PATHS_KMERPATHMUXSEARCHER_H


#include "paths/MuxGraph.h"
#include "paths/ReadFillDatabase.h"
#include "paths/SubsumptionList.h"

#include "paths/MuxSearchResult.h"
#include "paths/MuxSearchAgent.h"
#include "paths/MuxSearchPolicy.h"

/// KmerPathMuxSearcher encapsulates the depth-first search algorithm for
/// walking an insert by hopping from each read to its minimal extensions.
/// It can be passed any number of MuxSearchPolicies; if given any policies,
/// the insert is walked requiring coverage by reads which the policies all
/// consider good.
///
/// For an overview of Mux-based insert walking, see paths/doc/mux.pdf

class KmerPathMuxSearcher {

public:
  KmerPathMuxSearcher( const MuxGraph* p_mug, 
		       const ReadFillDatabase* p_fillDB,
		       const SubsumptionList* p_subList,
		       const vec<int>* p_readLengthsInKmers,
		       int verbosity=0 ) 
    : mp_muxGraph    ( p_mug ),
      mp_readfillDB  ( p_fillDB ),
      mp_subList     ( p_subList ),
      mp_readLengths ( p_readLengthsInKmers ),
      m_min_perfect_match( 1 ),
      m_verbosity   ( verbosity ),
      m_search_limit( 0 ),
      m_answer_size_limit( 0 )
  {
    m_max_read_length = *max_element( mp_readLengths->begin(), 
				      mp_readLengths->end() );
    if( m_verbosity > 2 ) PRINT(m_max_read_length);
  }

  // We own all of our policies.  Delete them in destructor.
  ~KmerPathMuxSearcher() { ClearPolicies(); }

  // You may not copy a KmerPathMuxSearcher -- who would own the policies?
private:
  KmerPathMuxSearcher( const KmerPathMuxSearcher& );
  KmerPathMuxSearcher& operator= ( const KmerPathMuxSearcher& );


public:
  // All policies are derived from MuxSearchPolicy, and come from here:
  void AddPolicy( MuxSearchPolicy* MSP ) { m_policies.push_back( MSP ); }
  void AddPolicies( vec<MuxSearchPolicy*> MSPs ) {
    m_policies.insert( m_policies.end(), MSPs.begin(), MSPs.end() );
  }

  // Caller is responsible for deleting this policy.
  void RemovePolicy( MuxSearchPolicy* MSP ) { m_policies.EraseValue( MSP ); }

  // We own all of our policies, so delete them when we forget them.
  void ClearPolicies() {
    for( ; ! m_policies.empty(); m_policies.pop_back() )
      delete m_policies.back();
  }


  /// NOTE convention: mux walking goes right-to-left.  That is,
  /// the opening read is on the *right*, closing read the *left*.
  /// This is opposite the Breadth/Depth directions.
  ///
  /// ALSO NOTE this unusual (but mux-appropriate) convention:
  /// the length arguments to FindClosures are the min/max dist
  /// between the left-most points in the opening and closing reads.
  /// That is, they are pair separation plus length of left read.

  void FindClosures( const int closingReadId,
                     const int openingReadId, 
                     const int minAcceptableExtLength,
                     const int maxAcceptableExtLength,
                     MuxSearchResult &result ) const;

private:

  // Templatized, so the source code appears below in this file.
  template<typename MuxSearchAgent>
  void SearchDirector( MuxSearchAgent sa,
		       const int openingReadId, 
		       const int closingReadId,
		       const int minAcceptableExtLength,
		       const int maxAcceptableExtLength,
		       MuxSearchResult &result ) const;

public:

  void SetVerbosity( int verb ) { m_verbosity = verb; }

  // A search limit of zero means unlimited
  // As of this writing, one minute gets to about 10 million.
  void SetSearchLimit( longlong lim ) { m_search_limit = lim; }

  void SetAnswerSizeLimit( longlong lim ) { m_answer_size_limit = lim; }


  void SetMinPerfectMatch( int mpm ) { m_min_perfect_match = mpm; }

private:

  bool HitSearchLimit( MuxSearchResult& result ) const {
    if( ++result.num_states_explored == m_search_limit ) {
      result.hit_search_limit = true;
      if( m_verbosity>0 )
	cout << "Aborting -- hit search limit of " << m_search_limit << endl;
      return true;
    }
    else
      return false;
  }

  bool HitAnswerSizeLimit( MuxSearchResult& result ) const {
    if( m_answer_size_limit>0 && result.WalkGraph().size() >= m_answer_size_limit ) {
      result.hit_search_limit = true;
      if( m_verbosity>0 )
	cout << "Aborting -- answer grew to " << result.WalkGraph().size() << " nodes"
	     << endl;
      return true;
    }
    else
      return false;
  }


private:

  const MuxGraph* mp_muxGraph;
  const ReadFillDatabase* mp_readfillDB;
  const SubsumptionList* mp_subList;
  const vec<int>* mp_readLengths;   // read lengths, in kmers
  int m_max_read_length;      // the maximum of mp_readLengths
  int m_min_perfect_match;

  vec<MuxSearchPolicy*> m_policies;

  int m_verbosity;
  longlong m_search_limit;
  longlong m_answer_size_limit;

};



////////////////////////////////////////////////////////////////////////////
//
//  The main search function
//

// This is templatized over a MuxSearchAgent class, which implements
// the actual primitive operations of the search.  These are slightly
// different for a mux-search with and without additional policies.
// Templates allow us to reuse the search code with no run-time overhead.

template<typename MuxSearchAgent>
void KmerPathMuxSearcher::SearchDirector( MuxSearchAgent msa,
					  int openingReadId, 
					  int closingReadId,
					  int minAcceptableExtLength,
					  int maxAcceptableExtLength,
					  MuxSearchResult& result ) const {

  msa.Setup( openingReadId, closingReadId,
	     minAcceptableExtLength, maxAcceptableExtLength,
	     result );

  while( ! msa.Done() ) {

    msa.TopOfLoop(); // ForceAssert loop invariant; report at high verbosity
    
    // If the current final mux has no more extensions to explore,
    // then we are done with it: pop one level on the mux stack.
    if( ! msa.HasNextMux() ) {
      msa.Pop();
      continue;
    }

    msa.GetNextMux();

    if( msa.NextMuxGoesTooFar() ) // returns false if it might subsume a closer
      continue;

    if( msa.StateWorthRecording() && ! msa.IsStateNewRecordIfSo() ) {
      // If we know we've been here before:
      if( m_verbosity > 5 )
	cout << "Found a searched state -- not exploring again" << endl;
      msa.SaveIfConverged();
      if( HitAnswerSizeLimit(result) )
	break;
    }
    else {
      // If this will be our first visit to this state:
      if( HitSearchLimit( result ) )
	break;
      if( msa.MayIPush() ) {
	msa.PushMux();
	if( msa.MayIClose() ) {
	  msa.SaveIfClosed();
	  if( HitAnswerSizeLimit(result) )
	    break;
	}
	// Note that we keep exploring deeper even if this was a closure:
	// we might be in a repeat region where we'll close again later.
      }
    }
  }
}


#endif

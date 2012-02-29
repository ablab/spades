// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

#ifndef PATHS_MUXSEARCHAGENT_H
#define PATHS_MUXSEARCHAGENT_H

#include <algorithm>

#include "paths/Mux.h"
#include "paths/MuxGraph.h"
#include "paths/MuxWalkGraph.h"
#include "paths/ReadFillDatabase.h"
#include "paths/SubsumptionList.h"

#include "paths/MuxSearchResult.h"
#include "paths/MuxSearchState.h"
#include "paths/MuxSearchPolicy.h"

/// The search agent directed by the SearchDirector of a
/// KmerPathMuxSearcher.
///
/// For an overview of Mux-based insert walking, see paths/doc/mux.pdf
///
/// The basic agent, MuxSearchAgentBase, takes care of the stacks
/// of muxes, identifies and stores closures, and stops the search
/// when it reaches the max insert length.  Searching with this
/// agent is possible but slow, as search branches never coalesce.
///
/// Derived classes should prune the search tree when reaching a
/// state isomorphic to one they have previously processed, and may 
/// optinally prune on any other criteria (eg require coverage by 
/// "happy" reads, whatever "happy" means).
//
//  For things inside the innermost loop: inline as much as possible.
//  There are no virtual methods, so no run-time polymorphism cost.

class MuxSearchAgentBase {

  /////////////////////////////
  //
  //  Constructor:
  //

public:
  MuxSearchAgentBase( const MuxGraph* muxGraph, 
		      const ReadFillDatabase* readfillDB, 
		      const SubsumptionList* subList, 
		      int verbosity )
    : mp_muxGraph   ( muxGraph ),
      mp_readfillDB ( readfillDB ),
      mp_subList    ( subList ),
      m_verbosity   ( verbosity )
  { }

  // Run-time polymorphism used only very rarely -- can't afford it
  // for things called every time through the main loop.
  // The only reason we need this is to silence warning messages.
  virtual ~MuxSearchAgentBase() {}

  /////////////////////////////
  //
  //  Public member functions, called by KmerPathMuxSearcher::SearchDirector
  //

public:
  void Setup( int openingReadId, 
	      int closingReadId,
	      int minAcceptableExtLength,
	      int maxAcceptableExtLength,
	      MuxSearchResult& result );
  
  bool Done() { return ( search.path_to_here.empty() && 
			 search.to_explore[0].empty() ); }

  void TopOfLoop() {
    ForceAssert( search.to_explore.size() == search.path_to_here.size() + 1 );
    ForceAssert( ( dist_to_here==0 && search.dist_stack.empty() ) 
		 || ( dist_to_here == search.dist_stack.back() ) );

    if( m_verbosity>3 ) ReportSearchState();
  }

  bool HasNextMux() {
    extensions_remaining = &(search.to_explore.back());
    return ( ! extensions_remaining->empty() );
  }

  void Pop() {
    ForceAssert( ! search.path_to_here.empty() ); // neither HasNextMux() nor Done()
    search.to_explore.pop_back();
    dist_to_here -= search.path_to_here.back().GetNumKmers();
    search.dist_stack.pop_back();
    search.path_to_here.pop_back();
  }

  void GetNextMux() {
    final_mux = extensions_remaining->back();
    extensions_remaining->pop_back();
  }

  bool NextMuxGoesTooFar() {
    // Even if so, we need to keep it if it subsumes a closer.
    subsume_iter = closers_I_subsume.find(final_mux.GetPathId());
    final_mux_subsumes_something = subsume_iter != closers_I_subsume.end();

    if( ! final_mux_subsumes_something &&
	dist_to_here + final_mux.GetNumKmers() > maxAcceptableExtLength ) {
      if( m_verbosity > 5 )
	cout << "Not exploring " << final_mux.GetPathId()
	     << ", would make path too long" << endl;
      return true;
    }
    else
      return false;
  }

  // Is it worth the time and space to calculate the invariant of the 
  // state and store it in the set of seen states?  The actual calculation
  // and storage are the job of the derived classes, and I suppose
  // some derived class might want to override this, but it seems
  // pretty general, so I'll include it in the base.
  bool StateWorthRecording() {
    return ( mp_muxGraph->NumMuxesOf( final_mux.GetPathId() ) > 1 );
  }

  bool IsStateNewRecordIfSo() { return true; }

  void SaveIfConverged();

  // Guarantee: if MayIPush() returns true, the search will call PushMux().
  // So derived classes can be efficient by doing some of their
  // pushing within MayIPush() if they're going to return true.
  bool MayIPush() { return true; }

  void PushMux();

  // Guarantee: this is always called immediately after PushMux().
  // So derived classes may cache the answer, if it's known during
  // MayIPush() or PushMux().
  bool MayIClose() { return true; }

  void SaveIfClosed();


  /////////////////////////////
  //
  //  Private member functions:
  //

private:
  void SetupClosers();
  void SetupOpeners();

  void ReportSearchState() const;
  void PrintExploreSizes( const vec< vec<Mux> >& to_explore ) const;

  // rarely called; virtual OK
  virtual void MergeIntoExistingClosures( vec<MuxWalkGraph::Label>& graph_hits );
  virtual void AddClosure( int closer_trim );
  void AddMuxSkippedClosures( vec<SubsumptionRecord>& subsumptions );


  /////////////////////////////
  //
  //  Data:
  //

private:
  const MuxGraph* mp_muxGraph;
  const ReadFillDatabase* mp_readfillDB;
  const SubsumptionList* mp_subList;

protected:
  int m_verbosity;

private:  
  int leftReadId, rightReadId;
  int minAcceptableExtLength, maxAcceptableExtLength;
  
  // The genuine closers (fillings of the left read).
  vec<OrientedKmerPathId> closers;
  // The other possible closers: reads which subsume the above ones.
  // Look up a super-read, the map gives you the closers it subsumes.
  map< OrientedKmerPathId,vec<SubsumptionRecord> > closers_I_subsume;

protected:
  // Of course derived classes need to be able to fiddle with the result
  MuxSearchResult* mp_result;

private:
  // The various stacks that are used in the search
  // are bundled together in a MuxSearchState.
  MuxSearchState search;
  int dist_to_here;

  // Loop variables (set once per loop, for expected re-use later in that pass)
  Mux final_mux;
  vec<Mux>* extensions_remaining;
  map< OrientedKmerPathId,vec<SubsumptionRecord> >::iterator subsume_iter;
  bool final_mux_subsumes_something;

  /////////////////////////////
  //
  // Protected read-only accessors to search state internals:
  //
  // The search state is private, with protected read-only accessors below
  // for the use of derived classes.  I can't think of why a derived class 
  // would need write access, but perhaps that's just a failure of my 
  // imagination -- if so, modify this organization in the future.
protected:
  const MuxSearchState& GetSearch() { return search; }
  const int& GetDist() { return dist_to_here; }
  const Mux& GetFinalMux() { return final_mux; }
  const int& GetMinAcceptableExtLength() { return minAcceptableExtLength; }
  const int& GetMaxAcceptableExtLength() { return maxAcceptableExtLength; }
  const vec<OrientedKmerPathId>& GetClosers() { return closers; }
  const map< OrientedKmerPathId,vec<SubsumptionRecord> >& GetCloserSubs()
    { return closers_I_subsume; }

  // This helper isn't used at all in the base class.
  // But it's used in every derived class, so I'm putting it here
  // to avoid code duplication.
protected:
  longlong BasicHashOfPosition() {
    return ( GetFinalMux().GetPathId().GetIdRc() 
	     ^ (longlong(GetDist()+GetFinalMux().GetNumKmers()) << 32) );
  }


};




///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//
//
//   Derived MuxSearchAgent classes.
//
//




/////////////////////////////
///
///  This search agent prunes any pair<okpid, dist> it has already seen
///

class MuxSearchAgentSimple : public MuxSearchAgentBase {
public:
  MuxSearchAgentSimple( const MuxGraph* muxGraph, 
			const ReadFillDatabase* readfillDB, 
			const SubsumptionList* subList, 
			int verbosity )
    : MuxSearchAgentBase( muxGraph, readfillDB, subList, verbosity ) { }

private:
  set<longlong> seen;

public:
  void Setup( int opening, int closing, int min, int max, MuxSearchResult& result )
  {
    MuxSearchAgentBase::Setup( opening, closing, min, max, result );
    seen.clear();
  }

  bool IsStateNewRecordIfSo() {
    return( (seen.insert(BasicHashOfPosition())).second );
  }

  longlong BasicHashOfPosition() {
    longlong here = MuxSearchAgentBase::BasicHashOfPosition();
    if( m_verbosity > 6 )
      cout << "state hash " << hex << here << dec << endl;
    return here;
  }

};


/////////////////////////////
///
///  This search is constructed around a vector of policy classes,
///  each of which judges reads to be good or bad.  A read is good if
///  all of the policies agree that it is, and a search must be
///  covered by good reads.  If that becomes impossible (becuase the
///  last good read was longer than max-path-length ago), prune.

class MuxSearchAgentGoodReads : public MuxSearchAgentBase {
public:
  MuxSearchAgentGoodReads( const MuxGraph* muxGraph, 
			   const ReadFillDatabase* readfillDB, 
			   const SubsumptionList* subList, 
			   const vec<int>* readLengths,
			   const int max_read_length,
			   const int min_perfect_match,
			   int verbosity,
			   const vec<MuxSearchPolicy*>& policies )
    : MuxSearchAgentBase ( muxGraph, readfillDB, subList, verbosity ),
      mp_readLengths     ( readLengths ),
      m_max_read_length  ( max_read_length ),
      m_policies         ( policies ),
      m_min_perfect_match( min_perfect_match )
  {}


private:

  const vec<int>* mp_readLengths; 
  const int m_max_read_length;
  const vec<MuxSearchPolicy*>& m_policies;
  int m_min_perfect_match;

  set< pair<longlong,longlong> > seen;


  // Internal state vectors which push and pop with the search:

  vec<Bool> m_isGood;
  vec<int> m_lastGoodDist;  // cached dist_to_here of most recent good read

  bool m_mayClose; // may differ from m_isGood.back() for min_perfect_match reasons


  /////////////////////////////
  //
  //  Override base class methods:
  //

public:
  bool IsStateNewRecordIfSo() {
    return( (seen.insert(PolicyHashOfPosition())).second );
  }

  void Setup( int opening, int closing, int min, int max, MuxSearchResult& result )
  {
    MuxSearchAgentBase::Setup( opening, closing, min, max, result );
    seen.clear();
    SetupPolicies();
  }

  void TopOfLoop() {
    MuxSearchAgentBase::TopOfLoop();
    ForceAssert( m_isGood.size() == GetSearch().path_to_here.size() );
    ForceAssert( m_isGood.size() == m_lastGoodDist.size() );
  }

  void Pop() {
    MuxSearchAgentBase::Pop();
    m_isGood.pop_back();
    m_lastGoodDist.pop_back();
    PopPolicies();
  }

  bool MayIPush();

  bool MayIClose() {
    return m_mayClose;
  }

private:
  // These pass on the good/bad read designations to the result MuxWalkGraph
  void AddClosure(int);
  void MergeIntoExistingClosures( vec<MuxWalkGraph::Label>& graph_hits );

private:
  // helper functions:
  void SetupPolicies() {
    for( vec<MuxSearchPolicy*>::const_iterator pol = m_policies.begin();
	 pol != m_policies.end(); pol++ )
      (*pol)->Setup( GetMinAcceptableExtLength(), 
		     GetMaxAcceptableExtLength(), 
		     GetClosers(), GetCloserSubs() );
  }

  void PopPolicies() {
    for_each( m_policies.begin(), m_policies.end(), 
	      mem_fun( &MuxSearchPolicy::Pop ) );
  }

  longlong GoodReadsHash( const MuxSearchState& search,
			  const int dist_to_here,
			  const int min_dist );

  pair<longlong,longlong> PolicyHashOfPosition();

};


#endif

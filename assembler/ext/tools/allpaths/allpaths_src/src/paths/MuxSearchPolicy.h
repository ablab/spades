/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////


#ifndef PATHS_MUXSEARCHPOLICY
#define PATHS_MUXSEARCHPOLICY

#include "PairsHandler.h"

#include "paths/Mux.h"
#include "paths/MuxSearchState.h"
#include "paths/KmerPath.h"
#include "paths/MuxToPath.h"
#include "paths/MuxGraph.h"
#include "paths/OffsetTracker.h"
#include "paths/SubsumptionList.h"


#if __GNUC__ > 2
#include <ext/hash_set>
using __gnu_cxx::hash_set;
#else
#include <hash_set>
#endif


////////////////////////////////////////////////////////////////////////////
///
///  MuxSearchPolicy is a virtual class from which all policies derive.
///
///  For an overview of Mux-based insert walking, see paths/doc/mux.pdf

class MuxSearchPolicy {
public:

  enum Opinion { GOOD, BAD, PRUNE };

  virtual ~MuxSearchPolicy() {}

  // The policy's opportunity to record any information about
  // the insert to be walked.  Most policies ignore this.
  virtual void Setup( int minExtLength, int maxExtLength,
		      const vec<OrientedKmerPathId>& closers,
		      const map< OrientedKmerPathId,vec<SubsumptionRecord> >& 
		        closers_I_subsume ) {}

  // Return value indicates whether this policy thinks the given mux
  // is good.  In either case, the policy should update any internal
  // state to push the mux.
  //
  // The return value PRUNE indicates that this read is bad, and
  // moreover this policy will never decide that any read further
  // along this branch of the search is good, so prune the search now.
  // This has no effect on the search result, but may save time.

  virtual Opinion PushAndGiveOpinion( const Mux& final_mux, 
				      const MuxSearchState& search,
				      const vec<Bool>& is_good ) = 0;

  virtual void Pop() {}


  /// The owning MuxSearchAgent will consider two nodes the same if
  /// they end with the same mux at the same depth.  If a policy may
  /// distinguish between two such nodes, it needs to offer its own
  /// hash of the state to further tell things apart.
  ///
  /// One common amount of extra state to record is the good_reads_hash,
  /// which encodes all the muxes declared GOOD within max_read_length.
  /// The special return value 1 will be converted to the good_reads_hash.
  /// Any value other than 1 will be incorporated into the hash directly.

  virtual longlong HashOfPosition( const MuxSearchState& search,
				   const int& dist_to_here,
				   const Mux& final_mux,
				   const vec<Bool>& is_good ) { return 0; }
};

////////////////////////////////////////////////////////////////////////////
//
//  Derived Mux Search Policies, all named MSP_Whatever
//


/// The trivial policy (accept everything).  For testing purposes.
class MSP_Trivial : public MuxSearchPolicy {
public:
  Opinion PushAndGiveOpinion( const Mux& final_mux, 
			      const MuxSearchState& search,
			      const vec<Bool>& is_good )
  { return GOOD; }
};

/// The trivial policy (accept everything), but use the good_reads_hash.
/// This can be derived from if you want to use that hash.
class MSP_GoodReadsHash : public MSP_Trivial {
public:
  longlong HashOfPosition( const MuxSearchState& search,
			   const int& dist_to_here,
			   const Mux& final_mux,
			   const vec<Bool>& is_good ) 
  { return 1; }
};


/// MSP_PairedPairs is a MuxSearchPolicy that says reads are good
/// if they paired-pair overlap with an earlier good read;
/// paired-pair overlapping means both ends reads of the insert align.
/// Intended for libraries with very tight length distribution and
/// high depth of coverage.
//
//  It maintains no internal state.
//
//  It uses the good_reads_hash.

class MSP_PairedPairs : public MSP_GoodReadsHash {
public:

  MSP_PairedPairs( const vecKmerPath* pathsFw,
		   const vecKmerPath* pathsRc,
		   const phandler* pairs_handler,
		   int verbosity=0 )
    : mp_pathsFw( pathsFw ),
      mp_pathsRc( pathsRc ),
      mp_pairs_handler( pairs_handler ),
      m_verbosity( verbosity )
  { }

  // This MSP maintains no state, so PAGO just tests whether a read is good
  Opinion PushAndGiveOpinion( const Mux& final_mux, 
			      const MuxSearchState& search,
			      const vec<Bool>& is_good );

private:
  // helper function for PushAndGiveOpinion:
  bool PairedPairsPerfectMatch( const OrientedKmerPathId& okpid1,
				const OrientedKmerPathId& okpid2,
				const int segoff = -1) const;

private:
  const vecKmerPath* mp_pathsFw;
  const vecKmerPath* mp_pathsRc;
  const phandler* mp_pairs_handler;
  int m_verbosity;
};




/// MSP_TheseReadsOnly is a search policy that only likes a fixed
/// subset of reads.  (If reads were gap-filled, the vec<int> passed
/// should hold the fill ids, not the original read ids.)

class MSP_TheseReadsOnly : public MuxSearchPolicy {
public:

  MSP_TheseReadsOnly( const vec<int>& theseReads )
    : m_good_reads( theseReads.begin(), theseReads.end(), 
		    2*theseReads.size() )
                    // Keep hash table load factor under 50%.
  {}

  // This MSP maintains no state, so PAGO just tests whether a read is good
  Opinion PushAndGiveOpinion( const Mux& final_mux, 
			      const MuxSearchState& search,
			      const vec<Bool>& is_good ) {
    if( m_good_reads.count(final_mux.GetPathId().GetId()) )
      return GOOD;
    else
      return BAD;
  }

private:
  hash_set<int> m_good_reads;

};



/// MSP_TheseOrientedReadsOnly is like MSP_TheseReadsOnly, but sees fw/rc.

class MSP_TheseOrientedReadsOnly : public MuxSearchPolicy {
public:

  MSP_TheseOrientedReadsOnly( const vec<int>& opener_direction_reads,
			      const vec<int>& closer_direction_reads )
    : m_good_reads_fw( closer_direction_reads.begin(), 
		       closer_direction_reads.end(), 
		       2*closer_direction_reads.size() ),
      m_good_reads_rc( opener_direction_reads.begin(), 
		       opener_direction_reads.end(), 
		       2*opener_direction_reads.size() )
  {}

  // This MSP maintains no state, so PAGO just tests whether a read is good
  Opinion PushAndGiveOpinion( const Mux& final_mux, 
			      const MuxSearchState& search,
			      const vec<Bool>& is_good ) {
    bool found = final_mux.GetPathId().IsFw()
      ? m_good_reads_fw.count( final_mux.GetPathId().GetId() )
      : m_good_reads_rc.count( final_mux.GetPathId().GetId() );
    return( found ? GOOD : BAD );
  }

private:
  hash_set<int> m_good_reads_fw;
  hash_set<int> m_good_reads_rc;

};


/// MSP_TrackPath tracks the current path.
/// All paths are accepted; this is only useful if derived from
class MSP_TrackPath : public MuxSearchPolicy {
public:
  MSP_TrackPath( const vecKmerPath* pathsFw, 
		 const vecKmerPath* pathsRc )
    : m_pathifier( pathsFw, pathsRc ),
      m_kmerpaths( 1 )  // seed with an empty kmer path
  { }

  Opinion PushAndGiveOpinion( const Mux& final_mux, 
			      const MuxSearchState& search,
			      const vec<Bool>& is_good );
  
  void Pop() { m_kmerpaths.pop_back(); }

  const KmerPath& CurrentPath() { return m_kmerpaths.back(); }

protected:
  MuxToPath m_pathifier;
  vec<KmerPath> m_kmerpaths;
};


/// MSP_PrintEachPath is for debugging; all paths are accepted.
class MSP_PrintEachPath : public MSP_TrackPath {
public:

  MSP_PrintEachPath( const vecKmerPath* pathsFw, const vecKmerPath* pathsRc )
    : MSP_TrackPath( pathsFw, pathsRc ) { }

  Opinion PushAndGiveOpinion( const Mux& final_mux, 
			      const MuxSearchState& search,
			      const vec<Bool>& is_good ) {
    MSP_TrackPath::PushAndGiveOpinion( final_mux, search, is_good );
    cout << CurrentPath() << endl;
    return( GOOD );
  }
};

/// MSP_FateOfClosure is for debugging: constructor takes a
/// KmerPath which someone thinks is a closure (eg a depth
/// search); policy points out why we didn't find it.
class MSP_FateOfClosure : public MSP_TrackPath {
public:
  MSP_FateOfClosure( const KmerPath& closure, 
		     const vecKmerPath* pathsFw, 
		     const vecKmerPath* pathsRc,
		     const MuxGraph* theMuxGraph )
    : MSP_TrackPath( pathsFw, pathsRc ),
      m_closure( closure ),
      mp_pathsFw( pathsFw ),
      mp_pathsRc( pathsRc ),
      mp_muxGraph( theMuxGraph ),
      m_verbose( false )
  {
    cout << "MSP_FateOfClosure: Tracking fate of putative closure " 
	 << m_closure << endl;
  }

  MSP_FateOfClosure( const String& closure_string, 
		     const vecKmerPath* pathsFw, 
		     const vecKmerPath* pathsRc,
		     const MuxGraph* theMuxGraph )
    : MSP_TrackPath( pathsFw, pathsRc ),
      m_closure( KmerPath(closure_string) ),
      mp_pathsFw( pathsFw ),
      mp_pathsRc( pathsRc ),
      mp_muxGraph( theMuxGraph ),
      m_verbose( false )
  {
    cout << "MSP_FateOfClosure: Tracking fate of putative closure " 
	 << m_closure << endl;
    if( closure_string[closure_string.size()-1] == 'v' )
      SetVerbose( true );
  }

  void SetVerbose( bool verbose ) { m_verbose = verbose; }

  Opinion PushAndGiveOpinion( const Mux& final_mux, 
			      const MuxSearchState& search,
			      const vec<Bool>& is_good ); 

private:
  const KmerPath m_closure;
  const vecKmerPath* mp_pathsFw;
  const vecKmerPath* mp_pathsRc;
  const MuxGraph* mp_muxGraph;

  bool m_verbose;
  
  vec<Mux> muxes;
};



/// MSP_ReachCloserOrAbort is a policy build around an OffsetTracker.
/// It causes us to abort any search path as soon as it is
/// not possible to reach a closer.

class MSP_ReachCloserOrAbort : public MuxSearchPolicy {
public:

  MSP_ReachCloserOrAbort( const OffsetTracker* offsetTracker )
    : mp_offsetTracker( offsetTracker ) { }

  void Setup( int minExtLength, int maxExtLength,
	      const vec<OrientedKmerPathId>& closers,
	      const map< OrientedKmerPathId,vec<SubsumptionRecord> >& 
	        closers_I_subsume );


  Opinion PushAndGiveOpinion( const Mux& final_mux, 
			      const MuxSearchState& search,
			      const vec<Bool>& is_good );

private:

  const OffsetTracker* mp_offsetTracker;

  // preprocessed list of where we would be happy to get
  // (relative offsets for subsumed closers already dealt with)
  struct endpoint {
    endpoint( int pid, int min_d, int max_d )
      : path_id(pid), min_dist(min_d), max_dist(max_d) { }
    int path_id;
    int min_dist;
    int max_dist;
  };
  vec<endpoint> m_endpoints;

  bool IsReachable( const endpoint& goal,
		    const ImplicitOffsetVec& can_get_here,
		    int dist_so_far );

  // Should we look into offset information for this read?
  // We could instead require a constructor argument telling it
  // how many reads to ignore, but this seems both more flexible
  // and less demanding on the user.
  bool UseThisRead( const int& read_id ) {
    return( ! mp_offsetTracker->GetOffsetsTo(read_id).empty() );
  }

};






#endif

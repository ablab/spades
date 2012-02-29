/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2010) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef KMER_BASE_BROKER
#define KMER_BASE_BROKER

#include "Basevector.h"
#include "Bitvector.h"
#include "CoreTools.h"
#include "paths/KmerPath.h"
#include "paths/SuperBaseVector.h"

#include <algorithm>  // for set_union
#include <map>

/**
   Class: KmerBaseBrokerTemplate

   A class to answer questions about kmers and paths that require 
   knowing the underlying sequence (as opposed to just the <kmer ids>
   that make up a <kmer path>).  For example:
   
   - what is the base sequence of a given kmer, given its <kmer number>?
   - what is the base sequence of a <KmerPathInterval>?
   - what is the base sequence of a <KmerPath>?
   
   Quick-and-dirty -- keeps its own copy of all the big files
   (reads.{fastb,paths{_rc},pathsdb}) in memory.
   Want this more efficient in both time and space later.
   NOTE: there is now a constructor which takes the address of a preloaded
   paths{_rc} and pathsdb, thereby avoiding duplication of these structures.
   
   Note: regarding references to "high-quality k-mers", for ALLPATHS purposes
   all kmers are high-quality.
*/
template <class TAG>
class KmerBaseBrokerTemplate {
public:

  KmerBaseBrokerTemplate( ) 
    : self_owned(False) 
  { }

  /// Constructor that loaded paths, paths_rc, pathsdb and reads direct from file
  KmerBaseBrokerTemplate(String RunDir, int k, const String readsBase = "reads", const String pathsBase = "paths", bool Verbose = false) 
    : K( k ),
      bases( RunDir + "/" + readsBase + ".fastb" ),
      hqkmersfile( RunDir + "/" + readsBase + ".hqkmers.k" + ToString(K) ),
      verbose( Verbose )
  { 
    pathsp = new vecKmerPath( RunDir + "/" + readsBase + "." + pathsBase + ".k" + ToString(K) );
    paths_rcp = new vecKmerPath( RunDir + "/" + readsBase + "." + pathsBase + "_rc.k" + ToString(K) );
    vec<TAG>* nonconst_pathsdbp = new vec<TAG>;
    BREADX2( RunDir + "/" + readsBase + "." + pathsBase + "db" + ( TAG::IS_BIG ? "_big" : "") + ".k" + ToString(K), *nonconst_pathsdbp ); 
    pathsdbp = nonconst_pathsdbp;
    self_owned = True; 
  }
  
  void Initialize(String RunDir, int k, const String readsBase = "reads", const String pathsBase = "paths", bool Verbose = false)
  { 
    if ( self_owned ) {
      delete pathsp;
      delete paths_rcp;
      delete pathsdbp;
    }
    K = k;
    bases.ReadAll( RunDir + "/" + readsBase + ".fastb" );
    hqkmersfile = RunDir + "/" + readsBase + ".hqkmers.k" + ToString(K);
    verbose = Verbose;
    pathsp = new vecKmerPath( RunDir + "/" + readsBase + "." + pathsBase + ".k" + ToString(K) );
    paths_rcp = new vecKmerPath( RunDir + "/" + readsBase + "." + pathsBase + "_rc.k" + ToString(K) );
    vec<TAG>* nonconst_pathsdbp = new vec<TAG>;
    BREADX2( RunDir + "/" + readsBase + "." + pathsBase + "db" + (TAG::IS_BIG ? "_big" : "") + ".k" + ToString(K), *nonconst_pathsdbp ); 
    pathsdbp = nonconst_pathsdbp;
    self_owned = True; 
  }



  /// Constructor that takes preloaded paths, paths_rc, pathsdb
  /// but loads reads directly from file
  KmerBaseBrokerTemplate( String RunDir, int k, 
                  const vecKmerPath& paths, const vecKmerPath& paths_rc,
                  const vec<TAG>& pathsdb, 
                  const String readsBase = "reads", bool Verbose = false )
    : K( k ),
      bases( RunDir + "/" + readsBase + ".fastb" ),
      hqkmersfile( RunDir + "/"+ readsBase + ".hqkmers.k" + ToString(K) ),
      verbose( Verbose )
  { 
    hqkmersfile = "";
    pathsp = &paths;
    paths_rcp = &paths_rc;
    pathsdbp = &pathsdb;
    self_owned = False;
  }


  /// Constructor that takes preloaded paths, paths_rc, pathsdb and reads
  /// Ignores hqkmers
  KmerBaseBrokerTemplate( int k, 
                  const vecKmerPath& paths, const vecKmerPath& paths_rc,
                  const vec<TAG>& pathsdb, 
                  const vecbasevector& reads, bool Verbose = false )
    : K( k ),
      bases( reads ),
      verbose( Verbose )
  { 
    hqkmersfile = "";
    pathsp = &paths;
    paths_rcp = &paths_rc;
    pathsdbp = &pathsdb;
    self_owned = False;
    // ForceAssertEq( paths.size( ), paths_rc.size( ) );
    ForceAssertEq( static_cast<size_t>(paths.size( )), reads.size( ) );
  }

  KmerBaseBrokerTemplate( int k, 
                  const vecKmerPath& paths, const vecKmerPath& paths_rc,
                  const vec<TAG>& pathsdb, 
                  const String& reads_file,
                  bool Verbose = false )
    : K( k ),
      verbose( Verbose )
  { 
    bases.ReadAll(reads_file);
    hqkmersfile = "";
    pathsp = &paths;
    paths_rcp = &paths_rc;
    pathsdbp = &pathsdb;
    self_owned = False;
    // ForceAssertEq( paths.size( ), paths_rc.size( ) );
    ForceAssertEq( static_cast<size_t>(paths.size( )), bases.size( ) );
  }

  /// Constructor that takes preloaded paths, pathsdb and reads, ignores hqkmers
  /// Special case where there are no rc paths - unipaths for example
  KmerBaseBrokerTemplate( int k, 
                  const vecKmerPath& paths, 
                  const vec<TAG>& pathsdb, 
                  const vecbasevector& reads, bool Verbose = false )
    : K( k ),
      bases( reads ),
      verbose( Verbose )
  { 
    hqkmersfile = "";
    pathsp = &paths;
    paths_rcp = 0;
    pathsdbp = &pathsdb;
    self_owned = False;
  }


  // The following does not correctly set hqkmersfile.  This might cause problems.

  void Initialize( int k, 
                  const vecbasevector& Bases,
                  const vecKmerPath& paths, const vecKmerPath& paths_rc,
                  const vec<TAG>& pathsdb, 
                  bool Verbose = false )
  { K = k;
    bases = Bases;
    hqkmersfile = "";
    verbose = Verbose;
    pathsp = &paths;
    paths_rcp = &paths_rc;
    pathsdbp = &pathsdb;
    self_owned = False;
  }
  
  ~KmerBaseBrokerTemplate( )
  {
    if (self_owned) {    
      delete pathsp;
      delete paths_rcp;
      delete pathsdbp;
    }
  }

  int GetK() const { return K; }

  /// Method: Bases(kmer_id_t)
  /// Convert a <kmer number> to its sequence.
  /// See also: <ClearBasesCache()>.
  const basevector& Bases(kmer_id_t k) const;

  /// Method: Bases(KmerPathInterval)
  /// Convert an interval of kmers to its sequence
  basevector Bases(KmerPathInterval rpi) const;

  // Method: ClearBasesCache
  // <Bases(kmer_id_t)> caches and returns references to the cache.
  // This lets you clear the cache, which also invalidates those refs.
  void ClearBasesCache() { bases_cache.clear(); }

  /// Method: ToSequence
  /// Convert a KmerPath to its sequence-with-gaps.
  ///
  /// This will assert if passed a corrupted KmerPath, ie one
  /// not representable in base space because of bad kmer proximities.
  ///
  /// In the returned object, gap lengths are in bases, not kmers,
  /// so their lengths are possibly negative (as low as -K+1).
  /// Negative gaps whose length is determined by the bases are just
  /// returned as continuous sequence.
  SuperBaseVector ToSequence( const KmerPath& path, String name="" ) const;

  /// ToSequenceSafe: same as ToSequence, but just check to see if ToSequence
  /// would assert.  Return True if OK.

  Bool ToSequenceSafe( const KmerPath& path, String name="" ) const;

  /// Method: Seq
  /// Return the sequence of a gap-free <KmerPath>.
  basevector Seq( const KmerPath& path ) const;




  // These utilities are primarily useful for negative gap validation:
  // Return the minimum d>=0 such that k2 can start d bases past k1 start.
  // Returns 0 only if k1=k2; returns 1 if k2 can follow k1; etc.
  // Maximum possible return value is K.
  // Begins searching at d_min, if given.
  int MinOffset(kmer_id_t k1, kmer_id_t k2, int d_min=0) const;

  // Same for maximum offset.  Of course this only makes sense if
  // d_max < K.  The goal is to completely fix the size of small gaps.
  int MaxOffset(kmer_id_t k1, kmer_id_t k2, int d_max) const;

  // Is the given offset possible?
  bool PossibleOffset( kmer_id_t k1, kmer_id_t k2, int d ) const;

  /// If you have two kmers separated by a gap of fixed size <K,
  /// the intervening kmers are determined (but might not have numbers!).
  /// This returns the kmers that go in the gap, if it is easy to figure
  /// them out from (=those two kmers at that separation appear in) the reads.
  //  Philosophically doesn't belong here because it doesn't involve actual
  //  bases, but a KBB owns the ungapped paths and pathsdb.
  bool KmersBetween( kmer_id_t k1, kmer_id_t k2, int gapsize, KmerPath& ans ) const;

  // Certainly one should never copy one of these!
  KmerBaseBrokerTemplate(const KmerBaseBrokerTemplate<TAG>&) {}

private:

  int K;

  const vecKmerPath *pathsp, *paths_rcp;
  const vec<TAG>* pathsdbp;
  Bool self_owned;

  vecbasevector bases;
  mutable map<kmer_id_t, basevector> bases_cache;

  String hqkmersfile;  // remember this if needed for later kmerHQ loading
  mutable vecbitvector kmerHQ;  // ToReadBase SparseReads into this as needed

  bool verbose;


  struct cache_key {
    kmer_id_t m_k1, m_k2;
    int m_gapsize;
    cache_key(kmer_id_t k1, kmer_id_t k2, int g)
      : m_k1(k1), m_k2(k2), m_gapsize(g) { };
    friend 
    bool operator<(const cache_key& lhs, const cache_key& rhs) {
      return(lhs.m_k1 < rhs.m_k1 ||
	     (lhs.m_k1 == rhs.m_k1 &&
	      (lhs.m_gapsize < rhs.m_gapsize ||
	       (lhs.m_gapsize == rhs.m_gapsize &&
		(lhs.m_k2 < rhs.m_k2)))));
    }
  };

  mutable map<cache_key,KmerPath> between_cache;

  mutable map<cache_key,int> min_offset_cache;
  mutable map<cache_key,int> max_offset_cache;

};

typedef KmerBaseBrokerTemplate<tagged_rpint> KmerBaseBroker;
typedef KmerBaseBrokerTemplate<big_tagged_rpint> KmerBaseBrokerBig;

#endif

// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

#ifndef PATHS_KMERPATHDATABASE_H
#define PATHS_KMERPATHDATABASE_H

#include "paths/KmerPath.h"

/**
   Class: KmerPathDatabaseTemplate

   A <path iterval database>.
   
   An encapsulated vec<TAG>, which allows per-database
   caching for lookups.
*/
template <class TAG>
class KmerPathDatabaseTemplate
{
 public:
  KmerPathDatabaseTemplate();

  KmerPathDatabaseTemplate( const String& filename );
  KmerPathDatabaseTemplate( const KmerPathDatabaseTemplate<TAG>& other );
  KmerPathDatabaseTemplate( const vec<TAG>* p_vecOfTaggedRpints );
  KmerPathDatabaseTemplate( const vecKmerPath& fwdPaths,
                    const vecKmerPath& revPaths );
  KmerPathDatabaseTemplate( const vec<KmerPath>& fwdPaths,
                    const vec<KmerPath>& revPaths );
  
  KmerPathDatabaseTemplate& operator= ( const KmerPathDatabaseTemplate<TAG>& other );
  
  ~KmerPathDatabaseTemplate();


  // MethodDecl: Contains
  // Find all occurrences of the given kmer.  Return an array of indices
  // of the path intervals in the database that contain this kmer.
  void Contains( kmer_id_t kmerId, 
                 vec<path_interval_id_t>& answer, bool append = false ) const;




  void GetHighFrequencyKmers( vec<kmer_id_t>& hiFreqKmers, const copy_num_t minFreq ) const;

  const TAG& operator[] ( const size_t index ) const
  {
    return (*mp_taggedRpints)[ index ];
  }

  const vec<TAG>& ConstRef() const { return *mp_taggedRpints; }

  typename vec<TAG>::const_iterator Begin() const { return mp_taggedRpints->begin(); }
  typename vec<TAG>::const_iterator End()   const { return mp_taggedRpints->end(); }

  longlong size() const { return mp_taggedRpints->size(); }

  // The mp_rpints vector is written out to or read in from the given filename.
  void Write( const String& filename ) const;
  void Read( const String& filename );


 private:
  template <class KmerPathVector>
  void ConstructFromKmerPaths( const KmerPathVector& fwPaths, const KmerPathVector& rcPaths );
  
  void CopyFrom( const KmerPathDatabaseTemplate<TAG>& other );

  const vec<TAG>* mp_taggedRpints;
  bool m_ownRpints;
  mutable unsigned int m_to;
  mutable longlong m_hits;
  mutable longlong m_misses;
  mutable longlong m_missDistance;
};

typedef KmerPathDatabaseTemplate<tagged_rpint> KmerPathDatabase;
typedef KmerPathDatabaseTemplate<big_tagged_rpint> KmerPathDatabaseBig;

#endif
